/** Remove all possible noise from a microboone-like detector.
 *
 * This filter is a kitchen sink class and is a candidate for
 * factoring.
 */

#ifndef WIRECELLSIGPROC_OMNIBUSNOISEFILTER
#define WIRECELLSIGPROC_OMNIBUSNOISEFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellUtil/Waveform.h"

#include <vector>
#include <map>

namespace WireCellSigProc {
    class OmnibusNoiseFilter : public WireCell::IFrameFilter, public WireCell::IConfigurable {
    public:
	/// Create an OmnibusNoiseFilter.
	OmnibusNoiseFilter(WireCell::Waveform::Period readout_window = WireCell::Waveform::Period(0,10.0), 
			   int nsamples = 9600); // fixme: these numbers need to be configurable
	virtual ~OmnibusNoiseFilter();

	/// WireCell node method.
	virtual bool operator()(const input_pointer& in, output_pointer& out);

	/// Explicitly deliver a channel quality database before
	/// running this node.  This may be done implicitly via configure().
	void set_channel_quality(WireCell::IChannelNoiseDatabase::pointer noise_db) {
	    m_noisedb = noise_db;
	}

	virtual void configure(const WireCell::Configuration& config);
	virtual WireCell::Configuration default_configuration() const;

    private:
	
	WireCell::Waveform::Period m_period;
	int m_nsamples;
	WireCell::IChannelNoiseDatabase::pointer m_noisedb;
    };

}

#endif
