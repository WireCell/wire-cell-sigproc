#include "WireCellSigProc/OmnibusNoiseFilter.h"

#include "WireCellSigProc/Diagnostics.h"

#include "WireCellUtil/Response.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(OmnibusNoiseFilter, WireCell::SigProc::OmnibusNoiseFilter,
                 WireCell::IFrameFilter, WireCell::IConfigurable);

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusNoiseFilter::OmnibusNoiseFilter()
{
    configure(default_configuration());
}
OmnibusNoiseFilter::~OmnibusNoiseFilter()
{
}

void OmnibusNoiseFilter::configure(const WireCell::Configuration& config)
{
    auto jmm = config["maskmap"];

    for (auto name : jmm.getMemberNames()) {
        m_maskmap[name] = jmm[name].asString();
	//	std::cerr << name << " " << m_maskmap[name] << std::endl;
    }

    
}
WireCell::Configuration OmnibusNoiseFilter::default_configuration() const
{
    Configuration cfg1;
    cfg1["chirp"] = "bad";
    cfg1["noisy"] = "bad";
    
    Configuration cfg;
    cfg["maskmap"] = cfg1;
    
    return cfg;
}


bool OmnibusNoiseFilter::operator()(const input_pointer& in, output_pointer& out)
{
    // For now, just collect any and all masks and interpret them as "bad"
    Waveform::ChannelMaskMap input_cmm = in->masks();
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

    auto traces = in->traces();
    for (auto trace : *traces.get()) {
    	int ch = trace->channel();

	if (find(bad_channels.begin(),bad_channels.end(),ch)!=bad_channels.end()){
	  bychan[ch].resize(nsamples,0);
	  //std::cout << "Xin3 " << bychan[ch].at(10) << std::endl;
	}else{
	  bychan[ch] = trace->charge(); // copy
	}

    	IChannelFilter::signal_t& signal = bychan[ch]; // ref

	//	if (ch>=7200 && ch<7200+48){
	  //std::cout << "Xin1: " << ch << " " << signal.at(10) << std::endl;
	  for (auto filter : m_perchan) {
    	    auto masks = filter->apply(ch, signal);

	    Waveform::merge(cmm,masks,m_maskmap);
	    
	    //    	    for (auto const& it: masks) {
	    //  bad_regions = Waveform::merge(bad_regions, it.second);
	    // }
	  }
	  //std::cout << "Xin2: " << signal.at(10) << std::endl;
	  //	}
    }

    
    // int counter = 0;
    for (auto group : m_noisedb->coherent_channels()) {
      // std::cout << counter << " " << group.size() << std::endl;
      // counter ++;

      int flag = 1;

      IChannelFilter::channel_signals_t chgrp;
      for (auto ch : group) {	    // fix me: check if we don't actually have this channel
    	if (bychan.find(ch)==bychan.end()) {
	  flag = 0;
	}else{
	  chgrp[ch] = bychan[ch]; // copy...
	}
      }
      
      if (flag == 0) continue;
      
      // std::cout << "Xin1: " << chgrp.size() << std::endl;
      for (auto filter : m_grouped) {
    	auto masks = filter->apply(chgrp);

	Waveform::merge(cmm,masks,m_maskmap);
    	//for (auto const& it: masks) {
	// bad_regions = Waveform::merge(bad_regions, it.second);
	//	}
      }

      for (auto cs : chgrp) {
    	//std::cout << bychan[cs.first].at(0) << " ";
    	bychan[cs.first] = cs.second; // copy
    	//std::cout << bychan[cs.first].at(0) << " " << cs.second.at(0) << std::endl;
      }
    }

    // run status
    for (auto trace : *traces.get()) {
    	int ch = trace->channel();
    	IChannelFilter::signal_t& signal = bychan[ch]; // ref
	for (auto filter : m_perchan_status) {
    	    auto masks = filter->apply(ch, signal);

	    Waveform::merge(cmm,masks,m_maskmap);
	}
    }
    
    
    // pack up output
   
    
    
    //std::cout << "Xin: " << bad_regions.size() << std::endl;
    
    ITrace::vector itraces;
    for (auto cs : bychan) {
	// fixme: that tbin though
	SimpleTrace *trace = new SimpleTrace(cs.first, 0, cs.second);
	itraces.push_back(ITrace::pointer(trace));
    }
    SimpleFrame* sframe = new SimpleFrame(in->ident(), in->time(), itraces, in->tick(), cmm);
    out = IFrame::pointer(sframe);

    
    return true;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
