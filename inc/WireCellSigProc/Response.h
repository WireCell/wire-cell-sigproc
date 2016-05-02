#ifndef WIRECELLSIGPROC_RESPONSE
#define WIRECELLSIGPROC_RESPONSE

#include "WireCellUtil/Waveform.h"

namespace WireCellSigProc {

    namespace Response {

	/// The cold electronics response function.
	double coldelec(double time_us, double gain_mVfC=7.8, double shaping_us=1.0);

	class Generator {
	public:
	    virtual ~Generator();
	    virtual double operator()(double time_us) const = 0;

	    /// Lay down the function into a binned waveform.
	    WireCell::Waveform::timeseq_t generate(double tick_us=0.5, double begin_us=0.0, double end_us=10.0);
	};

	/// A functional object caching gain and shape.
	class ColdElec : public Generator {
	    const double _g, _s;
	public:
	    ColdElec(double gain_mVfC=7.8, double shaping_us=1.0);
	    virtual ~ColdElec();

	    // return the response at given time
	    virtual double operator()(double time_us) const;

	};

	/// A functional object giving the response as a function of time to a simple RC circuit.
	class SimpleRC : public Generator {
	    const double _width, _offset;
	public:
	    //
	    SimpleRC(double width_us, double offset_us=0.0);
	    virtual ~SimpleRC();

	    // return the response at a given time
	    virtual double operator()(double time_us) const;

	};

    }
}

#endif
