/** Things operating in frequency space
 */
#ifndef WIRECELLSIGPROC_FREQUENCY
#define WIRECELLSIGPROC_FREQUENCY

#include "WireCellUtil/Waveform.h"

namespace WireCellSigProc {

    // /// A pair of set gain and shaping times and their corresponding
    // /// correct values.
    // struct GainShapingCorrection {
    // 	double set_gain, set_shaping;
    // 	double cor_gain, cor_shaping;

    // 	// Units here are mV/fC and microsecond.
    // 	GainShapingCorrection(double sg=7.8, double ss=1, double cg=14.0, double cs=2.0)
    // 	    : set_gain(sg), set_shaping(ss), cor_gain(cg), cor_shaping(cs) { }
    // };

    // /// An attenuation to apply (multiply) to frequencies inside some band.
    // struct BandFilter {
    // 	WireCell::Waveform::Domain band; // frequency band to over which to apply attenuation
    // 	double attenuation;    // multiplied to each DFT frequency bin

    // 	BandFilter(WireCell::Waveform::Domain band, double atten = 0.0)
    // 	    : band(band), attenuation(atten) { }
    // };

    // namespace Frequency {

    // 	WireCell::Waveform::realseq_t apply(WireCell::Waveform::realseq_t wave, WireCell::Waveform::compseq_t filter);
	
    // 	/// Correct the given Fourier spectrum for gain/shaping.
    // 	void correct_gainshape(WireCell::Waveform::complex_t& spec, const GainShapingCorrection& gsc);

    // 	/// Apply the given band filter to the Fourier spectrum which spans give domain.
    // 	void filter_band(WireCell::Waveform::complex_t& spec, const WireCell::Waveform::Domain& domain, const BandFilter& bf);

    // }
}

#endif

