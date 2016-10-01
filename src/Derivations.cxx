#include "WireCellSigProc/Derivations.h"

#include <iostream>
using namespace WireCellSigProc;

WireCell::Waveform::realseq_t Derivations::CalcMedian(const WireCell::IChannelFilter::channel_signals_t& chansig){
  float max_rms = 0;

  //std::cout << "Xin3: " << chansig.size() << std::endl;
  for (auto it: chansig){
    int ch = it.first;
    WireCell::IChannelFilter::signal_t& signal = it.second;
    //std::cout << ch << " " << signal.size() << std::endl;
  }
  


  WireCell::Waveform::realseq_t medians;
  return medians;
}
