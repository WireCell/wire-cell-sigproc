#include "WireCellSigProc/OmnibusNoiseFilter.h"

#include "WireCellSigProc/Diagnostics.h"
#include "WireCellSigProc/Response.h"

using namespace WireCell;

using namespace WireCellSigProc;

OmnibusNoiseFilter::OmnibusNoiseFilter(WireCell::Waveform::Period readout_window,
				       int nsamples)
    : m_period(readout_window)
    , m_nsamples(nsamples)
{
}
OmnibusNoiseFilter::~OmnibusNoiseFilter()
{
}

void OmnibusNoiseFilter::configure(const WireCell::Configuration& config)
{
}
WireCell::Configuration OmnibusNoiseFilter::default_configuration() const
{
}


bool OmnibusNoiseFilter::operator()(const input_pointer& in, output_pointer& out)
{
    // NO HARD CODED MAGIC NUMBERS HERE!
    //...
    // (just hard coded magic functionality!)

    Diagnostics::Chirp check_chirp; // fixme, there are magic numbers hidden here
    Diagnostics::Partial check_partial; // fixme, here too.

    auto traces = in->traces();
    for (auto trace : *traces.get()) {
	int ch = trace->channel();

	// fixme: some channels are just bad can should be skipped.

	// get signal with nominal baseline correction
	float baseline = m_noisedb->nominal_baseline(ch);
	auto signal = trace->charge(); // copy
	Waveform::increase(signal, baseline);

	// get signal with nominal gain correction 
	float gc = m_noisedb->gain_correction(ch);
	auto signal_gc = signal; // copy, need to keep original signal
	Waveform::scale(signal_gc, gc);

	// determine if chirping
	Waveform::BinRange chirped_bins;
	bool is_chirp = check_chirp(signal_gc, chirped_bins.first, chirped_bins.second);
	
	auto spectrum = Waveform::dft(signal);
	bool is_partial = check_partial(spectrum); // Xin's "IS_RC()"

	if (!is_partial) {
	    Waveform::scale(spectrum, m_noisedb->rcrc(ch));
	}

	Waveform::scale(spectrum, m_noisedb->config(ch));

	Waveform::scale(spectrum, m_noisedb->noise(ch));

	signal = Waveform::idft(spectrum);

	//....

    }
}


