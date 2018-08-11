#include "WireCellSigProc/OmnibusNoiseFilter.h"

#include "WireCellSigProc/Diagnostics.h"

#include "WireCellUtil/Response.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellUtil/NamedFactory.h"
// #include "WireCellUtil/ExecMon.h" // debugging

#include "FrameUtils.h"          // fixme: needs to move to somewhere more useful.

#include <unordered_map>

WIRECELL_FACTORY(OmnibusNoiseFilter, WireCell::SigProc::OmnibusNoiseFilter,
                 WireCell::IFrameFilter, WireCell::IConfigurable);

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusNoiseFilter::OmnibusNoiseFilter(std::string intag, std::string outtag)
    : m_nsamples(0)
    , m_intag(intag)            // orig
    , m_outtag(outtag)          // raw
{
}
OmnibusNoiseFilter::~OmnibusNoiseFilter()
{
}

void OmnibusNoiseFilter::configure(const WireCell::Configuration& cfg)
{
    //std::cerr << "OmnibusNoiseFilter: configuring with:\n" << cfg << std::endl;
    m_nsamples = (size_t)get(cfg, "nsamples", (int)m_nsamples);

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
    cfg["nsamples"] = (int)m_nsamples;
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

    // ExecMon em("starting NoiseFilter");

    auto traces = wct::sigproc::tagged_traces(inframe, m_intag);
    if (traces.empty()) {
        std::cerr << "OmnibusNoiseFilter: warning: no traces for tag \"" << m_intag << "\"\n";
	return true;
    }
    // em("got tagged traces");

    if (m_nsamples) {
        std::cerr << "OmnibusNoiseFilter: will resize working waveforms from "
                  << traces.at(0)->charge().size()
                  << " to " << m_nsamples << std::endl;
    }
    else {
        // Warning: this implicitly assumes a dense frame (ie, all tbin=0 and all waveforms same size).
        // It also won't stop triggering a warning inside OneChannelNoise if there is a mismatch.
        m_nsamples = traces.at(0)->charge().size();
    }

    // For now, just collect any and all masks and interpret them as "bad".
    Waveform::ChannelMaskMap input_cmm = inframe->masks();
    Waveform::ChannelMaskMap cmm;
    Waveform::merge(cmm,input_cmm,m_maskmap);
    
    // Get the ones from database and then merge
    std::vector<int> bad_channels = m_noisedb->bad_channels();
    {
	Waveform::BinRange bad_bins;
	bad_bins.first = 0;
	bad_bins.second = (int) m_nsamples;
	Waveform::ChannelMasks temp;
	for (size_t i = 0; i< bad_channels.size();i++){
	    temp[bad_channels.at(i)].push_back(bad_bins);
	    //std::cout << temp.size() << " " << temp[bad_channels.at(i)].size() << std::endl;
	}
	Waveform::ChannelMaskMap temp_map;
	temp_map["bad"] = temp;
	Waveform::merge(cmm,temp_map,m_maskmap);
    }

    // em("Starting loop on traces");

    int nchanged_samples = 0;

    // Collect our working area indexed by channel.
    std::unordered_map<int, SimpleTrace*> bychan;
    for (auto trace : traces) {
    	int ch = trace->channel();

	// make working area directly in simple trace to avoid memory fragmentation
	SimpleTrace* signal = new SimpleTrace(ch, 0, m_nsamples);
	bychan[ch] = signal;

	// if good
	if (find(bad_channels.begin(), bad_channels.end(),ch) == bad_channels.end()) {

	    auto const& charge = trace->charge();
	    const size_t ncharges = charge.size();	    

            signal->charge().assign(charge.begin(), charge.begin() + std::min(m_nsamples, ncharges));
	    signal->charge().resize(m_nsamples, 0.0);
            
            if (ncharges != m_nsamples) {
                nchanged_samples += std::abs((int)m_nsamples - (int)ncharges);
            }

	}

        int filt_count=0;
        for (auto filter : m_perchan) {
            auto masks = filter->apply(ch, signal->charge());

	    // fixme: probably should assure these masks do not lead to out-of-bounds...

            Waveform::merge(cmm,masks,m_maskmap);
            ++filt_count;
        }
    }
    traces.clear();		// done with our copy of vector of shared pointers

    if (nchanged_samples) {
        std::cerr << "OmnibusNoiseFilter: warning, truncated or extended " << nchanged_samples << " samples\n";
    }

    // em("starting coherent loop");

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
                chgrp[ch] = bychan[ch]->charge(); // copy...
            }
        }
      
        if (flag == 0) continue;
      
        for (auto filter : m_grouped) {
            auto masks = filter->apply(chgrp);

            Waveform::merge(cmm,masks,m_maskmap);
        }

        for (auto cs : chgrp) {
	    // cs.second; // copy
            bychan[cs.first]->charge().assign(cs.second.begin(), cs.second.end());
        }
    }

    // em("starting run status");

    // run status
    for (auto& it : bychan) {
	const int ch = it.first;
    	IChannelFilter::signal_t& signal = it.second->charge();
	for (auto filter : m_perchan_status) {
    	    auto masks = filter->apply(ch, signal);

	    Waveform::merge(cmm,masks,m_maskmap);
	}
    }
    
    // em("starting packing output");

    ITrace::vector itraces;
    for (auto& cs : bychan) {    // fixme: that tbin though
        itraces.push_back(ITrace::pointer(cs.second));
    }

    // em("made ouput");
    bychan.clear();
    // em("cleared map");

    auto sframe = new SimpleFrame(inframe->ident(), inframe->time(), itraces, inframe->tick(), cmm);
    IFrame::trace_list_t indices(itraces.size());
    for (size_t ind=0; ind<itraces.size(); ++ind) {
        indices[ind] = ind;
    }
    sframe->tag_traces(m_outtag, indices);
    sframe->tag_frame("noisefilter");
    outframe = IFrame::pointer(sframe);

    // em("finished NoiseFilter");
    // std::cerr << em.summary() << std::endl;
    return true;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
