/** Inherently MicroBooNE-specific functions and classes
 */

#ifndef WIRECELLSIGPROC_MICROBOONE
#define WIRECELLSIGPROC_MICROBOONE

#include "WireCellUtil/Waveform.h"
#include "WireCellIface/IChannelFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellSigProc/Diagnostics.h"


namespace WireCell {
    namespace SigProc {
	namespace Microboone {

	    bool Chirp_raise_baseline(WireCell::Waveform::realseq_t& sig, int bin1, int bin2);
	    bool SignalFilter(WireCell::Waveform::realseq_t& sig);
	    float CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig);
	    bool RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig);

	    bool RemoveFilterFlags(WireCell::Waveform::realseq_t& sig);
	    bool NoisyFilterAlg(WireCell::Waveform::realseq_t& spec, float min_rms, float max_rms);

	    bool SignalProtection(WireCell::Waveform::realseq_t& sig, const WireCell::Waveform::compseq_t& respec, int res_offset, int pad_f, int pad_b);
	    bool Subtract_WScaling(WireCell::IChannelFilter::channel_signals_t& chansig, const WireCell::Waveform::realseq_t& medians);


	    /** Microboone style coherent noise subtraction.
	     *
	     * Fixme: in principle, this class could be general purpose
	     * for other detectors.  However, it uses the functions above
	     * which hard code microboone-isms.  If those
	     * microboone-specific parameters can be pulled out to a
	     * higher layer then this class can become generic and move
	     * outside of this file.
	     */
	    class CoherentNoiseSub : public WireCell::IChannelFilter { // no iconfigurable
	    public:

		CoherentNoiseSub();
		virtual ~CoherentNoiseSub();

		//// IChannelFilter interface

		/** Filter in place the signal `sig` from given `channel`. */
		virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

		/** Filter in place a group of signals together. */
		virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

		/// Direct injection of needed service interfaces.
		/** Set the sampling used when digitizing the waveform. */
		void set_channel_noisedb(WireCell::IChannelNoiseDatabase::pointer ndb) {
		    m_noisedb = ndb;
		}

            private:
		WireCell::IChannelNoiseDatabase::pointer m_noisedb;

	    };


	    /** Microboone style single channel noise subtraction.
	     *
	     * Fixme: in principle, this class could be general purpose
	     * for other detectors.  However, it uses the functions above
	     * which hard code microboone-isms.  If those
	     * microboone-specific parameters can be pulled out to a
	     * higher layer then this class can become generic and move
	     * outside of this file.
	     */

	    class OneChannelNoise : public WireCell::IChannelFilter, public WireCell::IConfigurable {
	    public:

		OneChannelNoise();
		virtual ~OneChannelNoise();

		//// IChannelFilter interface

		/** Filter in place the signal `sig` from given `channel`. */
		virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

		/** Filter in place a group of signals together. */
		virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

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

	    class ADCBitShift : public WireCell::IChannelFilter, public WireCell::IConfigurable {
	    public:
		ADCBitShift();
		virtual ~ADCBitShift();

			/** Filter in place the signal `sig` from given `channel`. */
		virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

		/** Filter in place a group of signals together. */
		virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

		virtual void configure(const WireCell::Configuration& config);
		virtual WireCell::Configuration default_configuration() const;
		
	    private:
		// Number of bits 12
		// How many ADC: 500
		// Threshold for ADC Bit shift: 7.5 sigma
		// Threshold for correction: 80%
	    };
	    
	}

    }

}

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
