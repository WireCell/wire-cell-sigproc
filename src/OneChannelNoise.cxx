#include "WireCellSigProc/OneChannelNoise.h"

using namespace WireCell;

using namespace WireCellSigProc;

OneChannelNoise::OneChannelNoise()
    : m_check_chirp() // fixme, there are magic numbers hidden here
    , m_check_partial() // fixme, here too.
{
}
OneChannelNoise::~OneChannelNoise()
{
}

void OneChannelNoise::configure(const WireCell::Configuration& config)
{
    // fixme!
}
WireCell::Configuration OneChannelNoise::default_configuration() const
{
}


Waveform::ChannelMaskMap OneChannelNoise::apply(int ch, signal_t& signal) const
{
    Waveform::ChannelMaskMap ret;

    // NO HARD CODED MAGIC NUMBERS HERE!
    //...
    // (just hard coded magic functionality!)

    // fixme: some channels are just bad can should be skipped.

    // get signal with nominal baseline correction
    float baseline = m_noisedb->nominal_baseline(ch);
    Waveform::increase(signal, baseline *(-1));

    // // // get signal with nominal gain correction 
    // float gc = m_noisedb->gain_correction(ch);
    // auto signal_gc = signal; // copy, need to keep original signal
    // Waveform::scale(signal_gc, gc);

    // // // determine if chirping
    // Waveform::BinRange chirped_bins;
    // bool is_chirp = m_check_chirp(signal_gc, chirped_bins.first, chirped_bins.second);
    // if (is_chirp) {
    // 	ret["chirp"][ch].push_back(chirped_bins);
    // }
	
    // auto spectrum = Waveform::dft(signal);
    // // bool is_partial = m_check_partial(spectrum); // Xin's "IS_RC()"

    // // if (!is_partial) {
    // // 	Waveform::scale(spectrum, m_noisedb->rcrc(ch));
    // // }

    // // Waveform::scale(spectrum, m_noisedb->config(ch));

    // // Waveform::scale(spectrum, m_noisedb->noise(ch));

    // // remove the DC component 
    // spectrum.front() = 0;
    // signal = Waveform::idft(spectrum);

    // fixme: still need to add final rebaselining and "still noisy finding"

    return ret;
}


Waveform::ChannelMaskMap OneChannelNoise::apply(channel_signals_t& chansig) const
{
    return Waveform::ChannelMaskMap();
}
