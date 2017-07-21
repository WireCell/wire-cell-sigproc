#include "WireCellUtil/Waveform.h"
#include "FrameUtils.h"

#include <iostream>

using namespace WireCell;

void wct::sigproc::dump_frame(WireCell::IFrame::pointer frame)
{
    auto traces = frame->traces();
    const size_t ntraces = traces->size();
    std::vector<double> means, rmses, lengths, tbins;
    for (auto trace : *traces) {
        auto const& charge = trace->charge();

        auto mr = Waveform::mean_rms(charge);
        means.push_back(mr.first);
        rmses.push_back(mr.second*mr.second);
        const int nsamps = charge.size();
        lengths.push_back(nsamps);
        tbins.push_back(trace->tbin());
        
        if (std::isnan(mr.second)) {
            std::cerr << "Frame: channel " << trace->channel() << " rms is NaN\n";
        }

        for (int ind=0; ind<nsamps; ++ind) {
            float val = charge[ind];
            if (std::isnan(val)) {
                std::cerr << "Frame: channel " << trace->channel() << " sample " << ind << " is NaN\n";
            }
            if (std::isinf(val)) {
                std::cerr << "Frame: channel " << trace->channel() << " sample " << ind << " is INF\n";
            }
        }
    }
    double meanmean = Waveform::sum(means)/ntraces;
    double totrms = sqrt(Waveform::sum(rmses));
    double meanlen = Waveform::sum(lengths)/ntraces;
    double meantbin = Waveform::sum(tbins)/ntraces;
    std::cerr << "Frame: " << ntraces << " traces,"
              << " <mean>=" << meanmean
              << " TotRMS=" << totrms
              << " <len>=" << meanlen
              << " <tbin>=" << meantbin
              << std::endl;
    for (auto it : frame->masks()) {
        std::cerr << "\t" << it.first << " : " << it.second.size() << std::endl;
    }
}
