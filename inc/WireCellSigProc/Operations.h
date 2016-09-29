#ifndef WIRECELLSIGPROC_OPERATIONS
#define WIRECELLSIGPROC_OPERATIONS

#include "WireCellUtil/Waveform.h"

namespace WireCellSigProc {
  namespace Operations {
    bool Chirp_raise_baseline(WireCell::Waveform::realseq_t& sig, int bin1, int bin2);
    bool SignalFilter(WireCell::Waveform::realseq_t& sig);
    float CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig);
    bool RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig);

    bool RemoveFilterFlags(WireCell::Waveform::realseq_t& sig);
    bool NoisyFilterAlg(WireCell::Waveform::realseq_t& spec, int ch);
  }
}


#endif
