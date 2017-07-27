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

    // For now, just collect any and all masks and interpret them as "bad"
    Waveform::ChannelMaskMap input_cmm = inframe->masks();
    Waveform::ChannelMaskMap cmm;
    //Waveform::ChannelMasks bad_regions;
    //cmm["bad"] = bad_regions;


    Waveform::merge(cmm,input_cmm,m_maskmap);
    
    //    for (auto const& it: input_cmm) {
    //	bad_regions = Waveform::merge(bad_regions, it.second);
    // }
    
    
    // Get the ones from database and then merge
    int nsamples = m_noisedb->number_samples();
    std::vector<int> bad_channels = m_noisedb->bad_channels();
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
    
    //    bad_regions = Waveform::merge(bad_regions, temp);
    // for (int i = 0; i< bad_channels.size();i++){
    //   std::cout << bad_regions[bad_channels.at(i)].size() << std::endl;
    // }

    Waveform::merge(cmm,temp_map,m_maskmap);

    std::map<int, IChannelFilter::signal_t> bychan;

    auto traces = wct::sigproc::tagged_traces(inframe, m_intag);
    if (traces.empty()) {
        std::cerr << "No traces for tag \"" << m_intag << "\"\n";
        THROW(ValueError() << errmsg{"OmnibusNoiseFilter: no traces found"});
    }
    for (auto trace : traces) {
    	int ch = trace->channel();

	if (find(bad_channels.begin(),bad_channels.end(),ch)!=bad_channels.end()){
            bychan[ch].resize(nsamples,0);
            //std::cout << "Xin3 " << bychan[ch].at(10) << std::endl;
	}
        else{
            bychan[ch] = trace->charge(); // copy
	}

    	IChannelFilter::signal_t& signal = bychan[ch]; // ref

        int filt_count=0;
        for (auto filter : m_perchan) {
            auto masks = filter->apply(ch, signal);

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
