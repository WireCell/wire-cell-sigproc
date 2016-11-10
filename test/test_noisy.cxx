#include "WireCellSigProc/Diagnostics.h"
#include "WireCellSigProc/Microboone.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Testing.h"

#include <iostream>
#include <string>


// provides vectors "horig" and "hfilt"
// generate like:
// $ dump-root-hist-to-waveform sigproc/test/example-chirp.root horig hfilt > sigproc/test/example-chirp.h
#include <vector>
#include "example-noisy.h"

using namespace std;

using namespace WireCell;
using namespace WireCellSigProc;

int main(int argc, char* argv[])
{
    int ch = 0;
    Microboone::SignalFilter(horig);
    bool is_noisy = Microboone::NoisyFilterAlg(horig,ch);
    assert(is_noisy);
}
