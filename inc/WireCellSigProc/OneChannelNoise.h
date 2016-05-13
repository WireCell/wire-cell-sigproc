#ifndef WIRECELLSIGPROC_ONECHANNELNOISE
#define WIRECELLSIGPROC_ONECHANNELNOISE

#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelFilter.h"

#include "WireCellSigProc/Diagnostics.h"


namespace WireCellSigProc {

    class OneChannelNoise : public WireCell::IChannelFilter, public WireCell::IConfigurable {
    public:

	OneChannelNoise();
	virtual ~OneChannelNoise();

	//// IChannelFilter interface

	/** Filter in place the signal `sig` from given `channel`. */
	virtual void apply(int channel, signal_t& sig) const;

	/** Filter in place a group of signals together. */
	virtual void apply(channel_signals_t& chansig) const { return; }

	/// IConfigurable configuration interface
	virtual void configure(const WireCell::Configuration& config);
	virtual WireCell::Configuration default_configuration() const;

	/// Direct injection of needed service interfaces.
	/** Set the sampling used when digitizing the waveform. */
	void set_channel_noisedb(WireCell::IChannelNoiseDatabase::pointer ndb) {
	    m_noisedb = ndb;
	}

    private:

	Diagnostics::Chirp m_check_chirp; // fixme, these should be done via service interfaces
	Diagnostics::Partial m_check_partial; // at least need to expose them to configuration
	WireCell::IChannelNoiseDatabase::pointer m_noisedb;
    };

}

#endif
