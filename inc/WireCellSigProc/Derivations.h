#ifndef WIRECELLSIGPROC_DERIVATIONS
#define WIRECELLSIGPROC_DERIVATIONS

#include "WireCellUtil/Waveform.h"
#include "WireCellIface/IChannelFilter.h"

namespace WireCellSigProc{
  namespace Derivations{
    WireCell::Waveform::realseq_t CalcMedian(const WireCell::IChannelFilter::channel_signals_t& chansig);
  }
}

#endif
