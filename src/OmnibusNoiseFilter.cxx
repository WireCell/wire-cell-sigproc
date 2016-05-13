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
    std::map<int, IChannelFilter::signal_t> bychan;

    auto traces = in->traces();
    for (auto trace : *traces.get()) {
	int ch = trace->channel();

	bychan[ch] = trace->charge(); // copy
	IChannelFilter::signal_t& signal = bychan[ch]; // ref

	for (auto filter : m_perchan) {
	    filter->apply(ch, signal);
	}
    }
    for (auto group : m_noisedb->coherent_channels()) {
	IChannelFilter::channel_signals_t chgrp;
	for (auto ch : group) {	    // fix me: check if we don't actually have this channel
	    chgrp[ch] = bychan[ch]; // copy...
	}
	for (auto filter : m_perchan) {
	    filter->apply(chgrp);
	}
	for (auto cs : chgrp) {
	    bychan[cs.first] = cs.second; // copy
	}
    }

    {
	// pack up output
	ITrace::vector itraces;
	for (auto cs : bychan) {
	    // fixme: that tbin though
	    SimpleTrace *trace = new SimpleTrace(cs.first, 0, cs.second);
	    itraces.push_back(ITrace::pointer(trace));
	}
	SimpleFrame* sframe = new SimpleFrame(in->ident(), in->time(), itraces, in->tick());
	out = IFrame::pointer(sframe);
    }
}


