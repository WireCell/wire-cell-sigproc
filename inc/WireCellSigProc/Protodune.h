/** Inherently protoDUNE-specific functions and classes
 *  Modified from Microboone.h
 */

#ifndef WIRECELLSIGPROC_PROTODUNE
#define WIRECELLSIGPROC_PROTODUNE

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Bits.h"
#include "WireCellIface/IChannelFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IAnodePlane.h"

#include "WireCellSigProc/Diagnostics.h"


namespace WireCell {
    namespace SigProc {
	namespace Protodune {

		bool LinearInterpSticky(WireCell::Waveform::realseq_t& signal, std::vector<std::pair<int,int> >& st_ranges);
		bool FftInterpSticky(WireCell::Waveform::realseq_t& signal, std::vector<std::pair<int,int> >& st_ranges);
		bool FftShiftSticky(WireCell::Waveform::realseq_t& signal, double toffset, std::vector<std::pair<int,int> >& st_ranges);
		bool FftScaling(WireCell::Waveform::realseq_t& signal, int nsamples);

        // hold common config stuff
        class ConfigFilterBase : public WireCell::IConfigurable {
        public:

		ConfigFilterBase(const std::string& anode = "AnodePlane",
                         const std::string& noisedb = "OmniChannelNoiseDB");
        virtual ~ConfigFilterBase();

		/// IConfigurable configuration interface
		virtual void configure(const WireCell::Configuration& config);
		virtual WireCell::Configuration default_configuration() const;

        // FIXME: this method needs to die.
		void set_channel_noisedb(WireCell::IChannelNoiseDatabase::pointer ndb) {
		    m_noisedb = ndb;
		}
        protected:
		std::string m_anode_tn, m_noisedb_tn;
		IAnodePlane::pointer m_anode;
		IChannelNoiseDatabase::pointer m_noisedb;
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

	    class StickyCodeMitig : public WireCell::IChannelFilter, public ConfigFilterBase {
	    public:

		StickyCodeMitig(const std::string& anode_tn = "AnodePlane",
                        const std::string& noisedb = "OmniChannelNoiseDB");
		virtual ~StickyCodeMitig();

		//// IChannelFilter interface

		/** Filter in place the signal `sig` from given `channel`. */
		virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

		/** Filter in place a group of signals together. */
		virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

		
	    private:

		// Diagnostics::Chirp m_check_chirp; // fixme, these should be done via service interfaces
		// Diagnostics::Partial m_check_partial; // at least need to expose them to configuration
                
	    };

	    class FembClockReSmp : public WireCell::IChannelFilter, public ConfigFilterBase {
	    public:

		FembClockReSmp(const std::string& anode_tn = "AnodePlane",
                       const std::string& noisedb = "OmniChannelNoiseDB");
		virtual ~FembClockReSmp();

		//// IChannelFilter interface

		/** Filter in place the signal `sig` from given `channel`. */
		virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

		/** Filter in place a group of signals together. */
		virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

	    };

   
	}

    }

}

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
