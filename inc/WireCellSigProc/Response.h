#ifndef WIRECELLSIGPROC_RESPONSE
#define WIRECELLSIGPROC_RESPONSE

#include "WireCellUtil/Waveform.h"

namespace WireCellSigProc {

    namespace Response {

	/// The cold electronics response function.
	double coldelec(double time, double gain=7.8, double shaping=1.0);

	class Generator {
	public:
	    virtual ~Generator();
	    virtual double operator()(double time) const = 0;

	    /// Lay down the function into a binned waveform.
	    WireCell::Waveform::realseq_t generate(const WireCell::Waveform::Domain& domain, int nsamples);
	};

	/// A functional object caching gain and shape.
	class ColdElec : public Generator {
	    const double _g, _s;
	public:
	    // Create cold electronics response function.  Gain is an
	    // arbitrary scale, typically in mV/fC and shaping time in
	    // microsecond.  Shaping time in units consistent with
	    // calling the function.
	    ColdElec(double gain=7.8, double shaping=1.0);
	    virtual ~ColdElec();

	    // Return the response at given time.  Time in units consistent with shaping.
	    virtual double operator()(double time) const;

	};

	/// A functional object giving the response as a function of
	/// time to a simple RC circuit.
	class SimpleRC : public Generator {
	    const double _width, _offset;
	public:
	    // Create (current) response function for a simple RC
	    // circuit where a unit of charge is placed on the cap at
	    // time offset and circuit has RC time constant of given
	    // width.  Times are in units consistent with value used
	    // to call the function.
	    SimpleRC(double width, double offset=0.0);
	    virtual ~SimpleRC();

	    // Return the response at a given time.  Time in units
	    // consistent with width and offset.  Warning: to get the
	    // delta function, one must call *exactly* at the offset
	    // time.
	    virtual double operator()(double time) const;

	};

    }
}

#endif
