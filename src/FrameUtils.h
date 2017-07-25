#ifndef WIRECELL_SIGPROC_PRIVATE_FRAMEUTILS
#define WIRECELL_SIGPROC_PRIVATE_FRAMEUTILS

#include "WireCellIface/IFrame.h"
#include "WireCellIface/ITrace.h"
#include "WireCellUtil/Binning.h"
#include "WireCellUtil/Array.h"

#include <vector>

namespace wct {                 // eventually make this namespace
    namespace sigproc {         // pattern global.

        // Print some info to cerr about frame
        void dump_frame(WireCell::IFrame::pointer frame);


        // Raster a collection of traces into a 2D array block.  Each
        // row corresponds to one channel as indicated by the channels
        // vector.  The trace's tbin is used as an offset from column
        // 0.  If the array block is undersized, missed samples will
        // be quietly ignored.  Trace charge is added to any existing
        // values in the array block.
        void raster(WireCell::Array::array_xxf& block,
                    WireCell::ITrace::vector traces,
                    const std::vector<int>& channels);
         
        // Return a baseline calculated on a collection of traces as
        // the most probable binned sample value.  By default use a 12
        // bit ADC Binning.  The index of the binning is returned.
        int maxcount_baseline(const WireCell::ITrace::vector& traces, 
                              const WireCell::Binning& binning = WireCell::Binning(4096,0,4096));

    }
}

#endif

