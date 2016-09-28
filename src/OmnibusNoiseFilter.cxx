#include "WireCellSigProc/OmnibusNoiseFilter.h"

#include "WireCellSigProc/Diagnostics.h"
#include "WireCellSigProc/Response.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

using namespace WireCell;

using namespace WireCellSigProc;

OmnibusNoiseFilter::OmnibusNoiseFilter()
{
}
OmnibusNoiseFilter::~OmnibusNoiseFilter()
{
}

void OmnibusNoiseFilter::configure(const WireCell::Configuration& config)
{
}
WireCell::Configuration OmnibusNoiseFilter::default_configuration() const
{
}


bool OmnibusNoiseFilter::operator()(const input_pointer& in, output_pointer& out)
{
    // For now, just collect any and all masks and interpret them as "bad"
    Waveform::ChannelMaskMap input_cmm = in->masks();
    Waveform::ChannelMasks bad_regions;
    for (auto const& it: input_cmm) {
	bad_regions = Waveform::merge(bad_regions, it.second);
    }

    // Get the ones from database and then merge
    int nsamples = m_noisedb->number_samples();
    std::vector<int> bad_channels = m_noisedb->bad_channels();
    Waveform::BinRange bad_bins;
    bad_bins.first = 0;
    bad_bins.second = nsamples;
    Waveform::ChannelMasks temp;
    for (int i = 0; i< bad_channels.size();i++){
      temp[bad_channels.at(i)].push_back(bad_bins);
      //std::cout << temp.size() << " " << temp[bad_channels.at(i)].size() << std::endl;
    }
    bad_regions = Waveform::merge(bad_regions, temp);
    // for (int i = 0; i< bad_channels.size();i++){
    //   std::cout << bad_regions[bad_channels.at(i)].size() << std::endl;
    // }

    std::map<int, IChannelFilter::signal_t> bychan;

    auto traces = in->traces();
    for (auto trace : *traces.get()) {
	int ch = trace->channel();

	bychan[ch] = trace->charge(); // copy
	IChannelFilter::signal_t& signal = bychan[ch]; // ref

	for (auto filter : m_perchan) {
	    auto masks = filter->apply(ch, signal);
	    for (auto const& it: masks) {
		bad_regions = Waveform::merge(bad_regions, it.second);
	    }
	}
    }


    // for (auto group : m_noisedb->coherent_channels()) {
    // 	IChannelFilter::channel_signals_t chgrp;
    // 	for (auto ch : group) {	    // fix me: check if we don't actually have this channel
    // 	    chgrp[ch] = bychan[ch]; // copy...
    // 	}
    // 	for (auto filter : m_perchan) {
    // 	    auto masks = filter->apply(chgrp);
    //      for (auto const& it: masks) {
    //          bad_regions = Waveform::merge(bad_regions, it.second);
    //      }
    // 	}
    // 	for (auto cs : chgrp) {
    // 	    bychan[cs.first] = cs.second; // copy
    // 	}
    // }

    {
	// pack up output
	Waveform::ChannelMaskMap cmm;
	cmm["bad"] = bad_regions;

	//std::cout << "Xin: " << bad_regions.size() << std::endl;

	ITrace::vector itraces;
	for (auto cs : bychan) {
	    // fixme: that tbin though
	    SimpleTrace *trace = new SimpleTrace(cs.first, 0, cs.second);
	    itraces.push_back(ITrace::pointer(trace));
	}
	SimpleFrame* sframe = new SimpleFrame(in->ident(), in->time(), itraces, in->tick(), cmm);
	out = IFrame::pointer(sframe);
    }
}


