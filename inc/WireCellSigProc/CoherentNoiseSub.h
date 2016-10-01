#ifndef WIRECELLSIGPROC_COHERENTNOISESUB
#define WIRECELLSIGPROC_COHERENTNOISESUB

#include "WireCellIface/IChannelFilter.h"
#include "WireCellSigProc/Derivations.h"
#include "WireCellSigProc/Operations.h"

namespace WireCellSigProc {

    class CoherentNoiseSub : public WireCell::IChannelFilter { // no iconfigurable
    public:

	CoherentNoiseSub();
	virtual ~CoherentNoiseSub();

	//// IChannelFilter interface

	/** Filter in place the signal `sig` from given `channel`. */
	virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

	/** Filter in place a group of signals together. */
	virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

    };

}

#endif
