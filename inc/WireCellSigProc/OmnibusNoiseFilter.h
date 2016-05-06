/** Remove all possible noise from a microboone-like detector.
 *
 * This filter is a kitchen sink class and is a candidate for
 * factoring.
 */

#ifndef WIRECELLSIGPROC_OMNIBUSNOISEFILTER
#define WIRECELLSIGPROC_OMNIBUSNOISEFILTER

#include "WireCellIface/IFrameFilter.h"

#include <vector>
#include <map>

namespace WireCellSigProc {

    // fixme: make honest interface
    class IChannelQualityDatabase {
    public:
	virtual ~IChannelQualityDatabase() {};

	/// Return baseline correction (additive offset)
	virtual double baseline(int channel) const = 0;

	/// Return gain correction (multiplicative scaling);
	virtual double gain(int channel) const = 0;
	
	/// Return correction for electronics response (ie, RC+RC
	/// shaping) This is multiplied to signal in frequency space.
	virtual const WireCell::Waveform::complex_t& electronics(int channel) const  = 0;

	/// Return an overall response correction (eg, to fix
	/// incorrectly set gain/shape parameters).  This is
	/// multiplied to signal in frequency space.
	virtual const WireCell::Waveform::complex_t& adhoc(int channel) const  = 0;

	/// Return multiplicative filter applied in frequency space.
	virtual const WireCell::Waveform::complex_t& noise_filter(int channel) const = 0;

    };

    class OmnibusNoiseFilter : public IFrameFilter {
    public:
	/// Create an OmnibusNoiseFilter.
	OmnibusNoiseFilter(WireCell::Waveform::Period readout_window = WireCell::Waveform::Period(0,10.0), 
			   int nsamples = 9600); // fixme: these numbers need to be configurable
	virtual ~OmnibusNoiseFilter();

	/// WireCell node method.
	virtual bool operator()(const input_pointer& in, output_pointer& out);

	/// Deliver a channel quality database  before running this node.
	void set_channel_quality(std::shared_ptr<IChannelQualityDatabase> cqdb) {
	    m_cqdb = cqdb;
	}


    private:
	
	WireCell::Waveform::Period m_period;
	int m_nsamples;
	std::shared_ptr<IChannelQualityDataBase> m_cqdb;
    };

}

#endif
