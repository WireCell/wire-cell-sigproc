#ifndef WIRECELLSIGPROC_COHERENTNOISESUB
#define WIRECELLSIGPROC_COHERENTNOISESUB

#include "WireCellIface/IChannelFilter.h"

namespace WireCellSigProc {

    class CoherentNoiseSub : public WireCell::IChannelFilter { // no iconfigurable
    public:

	CoherentNoiseSub();
	virtual ~CoherentNoiseSub();

	//// IChannelFilter interface

	/** Filter in place the signal `sig` from given `channel`. */
	virtual void apply(int channel, signal_t& sig) const {
	    return ;		// not implemented / no-op
	}

	/** Filter in place a group of signals together. */
	virtual void apply(channel_signals_t& chansig) const;

    };

}

#endif
