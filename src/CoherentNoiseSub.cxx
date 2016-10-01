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

  // calculate the RMS 
  
  // do the signal protection and adaptive baseline

  // calculate the scaling coefficient

  // Subtract 
  

    return WireCell::Waveform::ChannelMaskMap();		// not implemented
}
WireCell::Waveform::ChannelMaskMap
CoherentNoiseSub::apply(int channel, signal_t& sig) const
{
    return WireCell::Waveform::ChannelMaskMap();		// not implemented
}

