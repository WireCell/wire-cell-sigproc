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

bool OmnibusNoiseFilter::operator()(const input_pointer& in, output_pointer& out)
{
    Diagnostics::Chirp check_chirp; // fixme, there are magic numbers hidden here
    Diagnostics::Partial check_partial; // fixme, here too.

    auto rc1sig = Response::ColdElec( 7.8, 1.0).generate(m_period, m_samples)
    auto rc2sig = Response::ColdElec(14.0, 2.0).generate(m_period, m_samples)


    auto traces = in->traces();
    for (auto trace : *traces.get()) {
	int ch = trace->channel();

	// fixme: some channels are just bad can should be skipped.

	// get signal with nominal baseline correction
	float baseline = m_cqdb->baseline(ch);
	auto signal = trace->charge(); // copy
	Waveform::shift(signal, baseline);

	// get signal with nominal gain correction 
	float gc = m_cqdb->gain(ch);
	auto signal_gc = signal; // copy, need to keep original signal
	Waveform::scale(signal_gc, gc);

	// determine if chirping
	Waveform::BinRange chirped_bins;
	bool is_chirp = check_chirp(signal_gc, chirped_bins.first, chirped_bins.second);
	
	auto spectrum = Waveform::dft(signal);
	bool is_partial = check_partial(spectrum); // Xin's "IS_RC()"

	if (!is_partial) {
	    Waveform::scale(spectrum, m_cqdb->electronics(ch));
	}

	Waveform::scale(spectrum, m_cqdb->adhoc(ch));

	Waveform::scale(spectrum, m_cqdb->noise_filter(ch));

	signal = Waveform::idft(spectrum);

	//....

    }
}


