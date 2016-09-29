#include "WireCellSigProc/Response.h"
#include <cmath>


/*
  Cold Electronics response function.

  How was this function derived?

    1. Hucheng provided a transfer function of our electronics, which
       is the Laplace transformation of the shaping function.

    2. We did a anti-Laplace inverse of the shaping function

    3. The one you saw in the code is basically the result of the inversion

  - time_us is time in microsecond

  - gain_par is proportional to the gain, basically at 7.8 mV/fC, the
    peak of the shaping function should be at 7.8 mV/fC. In the code,
    you can find what value that I set to reach 14 mV/fC.

  - shaping_us is the shaping time (us)

  - the hard-coded numbers are the result of the inverting the
    Lapalace transformation in Mathematica.

 */


using namespace WireCellSigProc;


Response::Generator::~Generator()
{
}

WireCell::Waveform::realseq_t Response::Generator::generate(const WireCell::Waveform::Domain& domain, int nsamples)
{
    WireCell::Waveform::realseq_t ret(nsamples);
    const double tick = (domain.second-domain.first)/nsamples;
    for (int ind=0; ind < nsamples; ++ind) {
	double t = domain.first + ind*tick;
	ret[ind] = (*this)(t);
    }
    return ret;
}





double Response::coldelec(double time, double gain, double shaping)
{
    // leave this up to caller.
    // if (time_us <=0 || time_us >= 10) { // range of validity
    // 	return 0.0;
    // }

    const double reltime = time/shaping;

    // a scaling is needed to make the anti-Lapalace peak match the expected gain
    gain *= 10*1.012;

    return 4.31054*exp(-2.94809*reltime)*gain
	-2.6202*exp(-2.82833*reltime)*cos(1.19361*reltime)*gain
	-2.6202*exp(-2.82833*reltime)*cos(1.19361*reltime)*cos(2.38722*reltime)*gain
	+0.464924*exp(-2.40318*reltime)*cos(2.5928*reltime)*gain
	+0.464924*exp(-2.40318*reltime)*cos(2.5928*reltime)*cos(5.18561*reltime)*gain
	+0.762456*exp(-2.82833*reltime)*sin(1.19361*reltime)*gain
	-0.762456*exp(-2.82833*reltime)*cos(2.38722*reltime)*sin(1.19361*reltime)*gain
	+0.762456*exp(-2.82833*reltime)*cos(1.19361*reltime)*sin(2.38722*reltime)*gain
 	-2.620200*exp(-2.82833*reltime)*sin(1.19361*reltime)*sin(2.38722*reltime)*gain 
	-0.327684*exp(-2.40318*reltime)*sin(2.5928*reltime)*gain + 
	+0.327684*exp(-2.40318*reltime)*cos(5.18561*reltime)*sin(2.5928*reltime)*gain
	-0.327684*exp(-2.40318*reltime)*cos(2.5928*reltime)*sin(5.18561*reltime)*gain
	+0.464924*exp(-2.40318*reltime)*sin(2.5928*reltime)*sin(5.18561*reltime)*gain;
}

Response::ColdElec::ColdElec(double gain, double shaping)
    : _g(gain)
    , _s(shaping)
{
}
Response::ColdElec::~ColdElec()
{
}

double Response::ColdElec::operator()(double time) const
{
    return coldelec(time, _g, _s);
}


Response::SimpleRC::SimpleRC(double width, double tick, double offset)
  : _width(width), _tick(tick), _offset(offset)

{
}
Response::SimpleRC::~SimpleRC()
{
}
double Response::SimpleRC::operator()(double time) const
{
    double ret = -_tick/_width * exp(-(time-_offset)/_width);
    if (time == _offset) {	// this is a sketchy comparison
	ret += 1.0;		// delta function
    }
    return ret;
}


