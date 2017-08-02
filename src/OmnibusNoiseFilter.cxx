#include "WireCellSigProc/OmnibusNoiseFilter.h"

#include "WireCellSigProc/Diagnostics.h"

#include "WireCellUtil/Response.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellUtil/NamedFactory.h"

#include "FrameUtils.h"          // fixme: needs to move to somewhere more useful.

WIRECELL_FACTORY(OmnibusNoiseFilter, WireCell::SigProc::OmnibusNoiseFilter,
                 WireCell::IFrameFilter, WireCell::IConfigurable);

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusNoiseFilter::OmnibusNoiseFilter()
    : m_intag("orig")
    , m_outtag("raw")
{
}
OmnibusNoiseFilter::~OmnibusNoiseFilter()
{
}

void OmnibusNoiseFilter::configure(const WireCell::Configuration& cfg)
{
    //std::cerr << "OmnibusNoiseFilter: configuring with:\n" << cfg << std::endl;
    auto jmm = cfg["maskmap"];
    for (auto name : jmm.getMemberNames()) {
        m_maskmap[name] = jmm[name].asString();
	//	std::cerr << name << " " << m_maskmap[name] << std::endl;
    }

    for (auto jf : cfg["channel_filters"]) {
        auto filt = Factory::find_tn<IChannelFilter>(jf.asString());
        std::cerr << "OmnibusNoiseFilter: adding channel filter: " << m_perchan.size()
                  << " \"" << jf.asString() << "\"\n";
        m_perchan.push_back(filt);
    }
    for (auto jf : cfg["channel_status_filters"]) {
        auto filt = Factory::find_tn<IChannelFilter>(jf.asString());
        std::cerr << "OmnibusNoiseFilter: adding channel status filter: " << m_perchan_status.size()
                  << " \"" << jf.asString() << "\"\n";
        m_perchan_status.push_back(filt);
    }
    for (auto jf : cfg["grouped_filters"]) {
        auto filt = Factory::find_tn<IChannelFilter>(jf.asString());
        std::cerr << "OmnibusNoiseFilter: adding grouped filter: " << m_grouped.size()
                  << " \"" << jf.asString() << "\"\n";
        m_grouped.push_back(filt);
    }


    auto jcndb = cfg["noisedb"];
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(jcndb.asString());
    std::cerr << "OmnibusNoiseFilter: using channel noise DB object: " 
                  << " \"" << jcndb.asString() << "\"\n";

    m_intag = get(cfg, "intraces", m_intag);
    m_outtag = get(cfg, "outtraces", m_outtag);
}

WireCell::Configuration OmnibusNoiseFilter::default_configuration() const
{
    Configuration cfg;
    cfg["maskmap"]["chirp"] = "bad";
    cfg["maskmap"]["noisy"] = "bad";
    
    cfg["channel_filters"][0] = "mbOneChannelNoise";
    cfg["channel_status_filters"][0] = "mbOneChannelStatus";
    cfg["grouped_filters"][0] = "mbCoherentNoiseSub";

    // user must supply.  "OmniChannelNoiseDB" is a likely choice.
    // Avoid SimpleChannelNoiseDB.
    cfg["noisedb"] = ""; 

    // The tags for input and output traces
    cfg["intraces"] = m_intag;
    cfg["outtraces"] = m_outtag;
    return cfg;
}


bool OmnibusNoiseFilter::operator()(const input_pointer& inframe, output_pointer& outframe)
{
    if (!inframe) {             // eos
        outframe = nullptr;
        return true;
    }

    auto traces = wct::sigproc::tagged_traces(inframe, m_intag);
    if (traces.empty()) {
        std::cerr << "OmnibusNoiseFilter: warning: no traces for tag \"" << m_intag << "\"\n";
	return true;
    }

    // Warning: this implicitly assumes a dense frame (ie, all tbin=0 and all waveforms same size).
    const int nsamples = traces.at(0)->charge().size();

    // For now, just collect any and all masks and interpret them as "bad".
    Waveform::ChannelMaskMap input_cmm = inframe->masks();
    Waveform::ChannelMaskMap cmm;
    Waveform::merge(cmm,input_cmm,m_maskmap);
    
    // Get the ones from database and then merge
    std::vector<int> bad_channels = m_noisedb->bad_channels();
    {
	Waveform::BinRange bad_bins;
	bad_bins.first = 0;
	bad_bins.second = nsamples;
	Waveform::ChannelMasks temp;
	for (size_t i = 0; i< bad_channels.size();i++){
	    temp[bad_channels.at(i)].push_back(bad_bins);
	    //std::cout << temp.size() << " " << temp[bad_channels.at(i)].size() << std::endl;
	}
	Waveform::ChannelMaskMap temp_map;
	temp_map["bad"] = temp;
	Waveform::merge(cmm,temp_map,m_maskmap);
    }

    // Collect our working area indexed by channel.
    std::map<int, IChannelFilter::signal_t> bychan;
    for (auto trace : traces) {
    	int ch = trace->channel();

	// make local copy which will be filled only by good channel.
	bychan[ch].resize(nsamples,0);

	// good
	if (find(bad_channels.begin(), bad_channels.end(),ch) == bad_channels.end()) {

	    auto const& charge = trace->charge();
	    const int ncharges = charge.size();	    

	    { // sanity check.  This "should" never print.
		if (ncharges != nsamples) { // this block "should" never be called
		    std::cerr << "OmnibusNoiseFilter: found different length waveforms: "
			      << nsamples << " != " << ncharges << " in channel " << ch
			      << ". Will normalize to the first one."
			      << std::endl;
		}
	    }

	    // Do assignment with care not to overflow input nor output
	    bychan[ch].assign(charge.begin(), charge.begin() + std::min(nsamples, ncharges));
	}

    	IChannelFilter::signal_t& signal = bychan[ch]; // ref

        int filt_count=0;
        for (auto filter : m_perchan) {
            auto masks = filter->apply(ch, signal);

	    // fixme: probably should assure these masks do not lead to out-of-bounds...

            Waveform::merge(cmm,masks,m_maskmap);
            ++filt_count;
        }
    }

    int group_counter = 0;
    for (auto group : m_noisedb->coherent_channels()) {
        ++group_counter;

        int flag = 1;

        IChannelFilter::channel_signals_t chgrp;
        for (auto ch : group) {	    // fix me: check if we don't actually have this channel
	    // std::cout << group_counter << " " << ch << " " << std::endl;
            if (bychan.find(ch)==bychan.end()) {
                std::cerr << "OmnibusNoiseFilter: warning: unknown channel " << ch << "\n";
                flag = 0;
            }
            else{
                chgrp[ch] = bychan[ch]; // copy...
            }
        }
      
        if (flag == 0) continue;
      
        for (auto filter : m_grouped) {
            auto masks = filter->apply(chgrp);

            Waveform::merge(cmm,masks,m_maskmap);
        }

        for (auto cs : chgrp) {
            bychan[cs.first] = cs.second; // copy
        }
    }

    // run status
    for (auto trace : traces) {
    	int ch = trace->channel();
    	IChannelFilter::signal_t& signal = bychan[ch]; // ref
	for (auto filter : m_perchan_status) {
    	    auto masks = filter->apply(ch, signal);

	    Waveform::merge(cmm,masks,m_maskmap);
	}
    }
    
    // pack up output
    ITrace::vector itraces;
    for (auto cs : bychan) {    // fixme: that tbin though
        itraces.push_back(std::make_shared<SimpleTrace>(cs.first, 0, cs.second));
    }
    auto sframe = new SimpleFrame(inframe->ident(), inframe->time(), itraces, inframe->tick(), cmm);
    IFrame::trace_list_t indices(itraces.size());
    for (size_t ind=0; ind<itraces.size(); ++ind) {
        indices[ind] = ind;
    }
    sframe->tag_traces(m_outtag, indices);
    sframe->tag_frame("noisefilter");
    outframe = IFrame::pointer(sframe);
    return true;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
