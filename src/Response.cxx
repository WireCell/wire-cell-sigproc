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

WireCell::Waveform::timeseq_t Response::Generator::generate(double tick_us, double begin_us, double end_us)
{
    const int nticks = (end_us-begin_us)/tick_us;

    WireCell::Waveform::timeseq_t ret(nticks);
    for (int ind=0; ind < nticks; ++ind) {
	double t = begin_us + ind*tick_us;
	ret(ind) = (*this)(t);
    }
    return ret;
}





double Response::coldelec(double time_us, double gain_mVfC, double shaping_us)
{
    if (time_us <=0 || time_us >= 10) { // range of validity
	return 0.0;
    }

    // a scaling is needed to make the anti-Lapalace peak match the expected gain
    const double gain = gain_mVfC * 10*1.012;
    const double shaping  = shaping_us;

    return 4.31054*exp(-2.94809*time_us/shaping)*gain
	-2.6202*exp(-2.82833*time_us/shaping)*cos(1.19361*time_us/shaping)*gain
	-2.6202*exp(-2.82833*time_us/shaping)*cos(1.19361*time_us/shaping)*cos(2.38722*time_us/shaping)*gain
	+0.464924*exp(-2.40318*time_us/shaping)*cos(2.5928*time_us/shaping)*gain
	+0.464924*exp(-2.40318*time_us/shaping)*cos(2.5928*time_us/shaping)*cos(5.18561*time_us/shaping)*gain
	+0.762456*exp(-2.82833*time_us/shaping)*sin(1.19361*time_us/shaping)*gain
	-0.762456*exp(-2.82833*time_us/shaping)*cos(2.38722*time_us/shaping)*sin(1.19361*time_us/shaping)*gain
	+0.762456*exp(-2.82833*time_us/shaping)*cos(1.19361*time_us/shaping)*sin(2.38722*time_us/shaping)*gain
 	-2.620200*exp(-2.82833*time_us/shaping)*sin(1.19361*time_us/shaping)*sin(2.38722*time_us/shaping)*gain 
	-0.327684*exp(-2.40318*time_us/shaping)*sin(2.5928*time_us/shaping)*gain + 
	+0.327684*exp(-2.40318*time_us/shaping)*cos(5.18561*time_us/shaping)*sin(2.5928*time_us/shaping)*gain
	-0.327684*exp(-2.40318*time_us/shaping)*cos(2.5928*time_us/shaping)*sin(5.18561*time_us/shaping)*gain
	+0.464924*exp(-2.40318*time_us/shaping)*sin(2.5928*time_us/shaping)*sin(5.18561*time_us/shaping)*gain;
}

Response::ColdElec::ColdElec(double gain_mVfC, double shaping_us)
    : _g(gain_mVfC)
    , _s(shaping_us)
{
}
Response::ColdElec::~ColdElec()
{
}

double Response::ColdElec::operator()(double time_us) const
{
    return coldelec(time_us, _g, _s);
}


Response::SimpleRC::SimpleRC(double width_us, double offset_us)
    : _width(width_us), _offset(offset_us)

{
}
Response::SimpleRC::~SimpleRC()
{
}
double Response::SimpleRC::operator()(double time_us) const
{
    double ret = -1.0/_width * exp(-(time_us-_offset)/_width);
    if (time_us == _offset) {	// sketchy
	ret += 1.0;		// delta function
    }
    return ret;
}


