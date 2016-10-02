#include "WireCellSigProc/Diagnostics.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Testing.h"

#include <iostream>
#include <string>


// provides vectors "horig" and "hfilt"
// generate like:
// $ dump-root-hist-to-waveform sigproc/test/example-chirp.root horig hfilt > sigproc/test/example-chirp.h
#include <vector>
#include "example-chirping.h"

using namespace std;

using namespace WireCell;
using namespace WireCellSigProc;

int main(int argc, char* argv[])
{
  
   int beg=-1, end=-1;
    Diagnostics::Chirp chirp;
    bool found = chirp(horig, beg, end);
    Assert(found);
    
    // the function should find something starting at the beginning.
    Assert(beg == 0);
    Assert(end == 3720);
    return 0;
}
