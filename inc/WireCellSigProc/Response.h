#ifndef WIRECELLSIGPROC_RESPONSE
#define WIRECELLSIGPROC_RESPONSE

#include "WireCellUtil/Waveform.h"

namespace WireCellSigProc {

    namespace Response {

	/// The cold electronics response function.
	double coldelec(double time_us, double gain_mVfC=7.8, double shaping_us=1.0);

	/// A functional object caching gain and shape.
	class ColdElec {
	    double _g, _s;
	public:
	    ColdElec(double gain_mVfC=7.8, double shaping_us=1.0);
	    // return the response at given time
	    double operator()(double time_us) const;

	    /// Lay down the function into a binned waveform.
	    WireCell::Waveform::signal_t generate(double tick_us=0.5, double begin_us=0.0, double end_us=10.0);
	};
    }
}

#endif
