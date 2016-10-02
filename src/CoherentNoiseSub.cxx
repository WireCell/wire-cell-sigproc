#include "WireCellSigProc/CoherentNoiseSub.h"
#include <iostream>

using namespace WireCellSigProc;

CoherentNoiseSub::CoherentNoiseSub()
{
}
CoherentNoiseSub::~CoherentNoiseSub()
{
}

WireCell::Waveform::ChannelMaskMap
CoherentNoiseSub::apply(channel_signals_t& chansig) const
{
  // std::cout << "Xin2: " << std::endl;
  // find the median among all 
  WireCell::Waveform::realseq_t medians = Derivations::CalcMedian(chansig);
  //std::cout << medians.size() << " " << medians.at(0) << " " << medians.at(1) << std::endl;
  

  // do the signal protection and adaptive baseline
  Operations::SignalProtection(medians);

  // calculate the scaling coefficient and subtract
  Operations::Subtract_WScaling(chansig, medians);
  
  // for (auto it: chansig){
  //   std::cout << "Xin3 " << it.second.at(0) << std::endl;
  // }
  
  return WireCell::Waveform::ChannelMaskMap();		// not implemented
}
WireCell::Waveform::ChannelMaskMap
CoherentNoiseSub::apply(int channel, signal_t& sig) const
{
    return WireCell::Waveform::ChannelMaskMap();		// not implemented
}

