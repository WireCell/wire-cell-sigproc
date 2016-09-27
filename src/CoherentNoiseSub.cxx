#include "WireCellSigProc/CoherentNoiseSub.h"

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
    return WireCell::Waveform::ChannelMaskMap();		// not implemented
}
WireCell::Waveform::ChannelMaskMap
CoherentNoiseSub::apply(int channel, signal_t& sig) const
{
    return WireCell::Waveform::ChannelMaskMap();		// not implemented
}

