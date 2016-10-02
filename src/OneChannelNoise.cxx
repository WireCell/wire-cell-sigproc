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

    // get signal with nominal gain correction 
    float gc = m_noisedb->gain_correction(ch);
    // if (ch < 2400)
    //   std::cout << gc << " " << ch << std::endl;
    auto signal_gc = signal; // copy, need to keep original signal
    Waveform::scale(signal_gc, gc);

    // determine if chirping
    Waveform::BinRange chirped_bins;
    bool is_chirp = m_check_chirp(signal_gc, chirped_bins.first, chirped_bins.second);
    if (is_chirp) {
      ret["chirp"][ch].push_back(chirped_bins);
      // for (int i=chirped_bins.first;i!=chirped_bins.second;i++){
      // 	signal.at(i) = 0;
      // }
    }

    auto spectrum = Waveform::dft(signal);
    bool is_partial = m_check_partial(spectrum); // Xin's "IS_RC()"
    // if (is_partial){
    //   std::cout << ch << std::endl;
    // }
    
    if (!is_partial) {
      //std::cout << "Xin: " << spectrum.front().real() << " " ;
      Waveform::shrink(spectrum, m_noisedb->rcrc(ch));
      //std::cout << spectrum.front().real() << std::endl;
    }

    // if (ch==2000) std::cout << "2000" << " " << m_noisedb->config(ch).at(1) << " " << m_noisedb->gain_correction(ch) << std::endl;
    // if (ch==2016) std::cout << "2016" << " " << m_noisedb->config(ch).at(1) << " " << m_noisedb->gain_correction(ch) << std::endl;

    Waveform::scale(spectrum, m_noisedb->config(ch));
    Waveform::scale(spectrum, m_noisedb->noise(ch));

    // remove the DC component 
    spectrum.front() = 0;
    signal = Waveform::idft(spectrum);

    //Now calculate the baseline ...
    baseline = Waveform::median(signal);
    //correct baseline
    Waveform::increase(signal, baseline *(-1));

    // Now do adaptive baseline for the chirping channels
    if (is_chirp){
      Operations::Chirp_raise_baseline(signal,chirped_bins.first, chirped_bins.second);
      Operations::SignalFilter(signal);
      Operations::RawAdapativeBaselineAlg(signal);
    }
    // Now do the adaptive baseline for the bad RC channels
    if (is_partial){
      Operations::SignalFilter(signal);
      Operations::RawAdapativeBaselineAlg(signal);
    }

    // Identify the Noisy channels ... 
    Operations::SignalFilter(signal);
    bool is_noisy = Operations::NoisyFilterAlg(signal,ch);
    Operations::RemoveFilterFlags(signal);

    // if (is_noisy){
    //   std::cout << "Xin: " << signal.at(1) << std::endl;
    // }

    // std::cout << ch << " " << is_chirp << " " << is_partial << " " << is_noisy << std::endl;

    if (is_noisy){
      chirped_bins.first = 0;
      chirped_bins.second = signal.size();
      ret["noisy"][ch].push_back(chirped_bins);
    }

    return ret;
}


Waveform::ChannelMaskMap OneChannelNoise::apply(channel_signals_t& chansig) const
{
    return Waveform::ChannelMaskMap();
}
