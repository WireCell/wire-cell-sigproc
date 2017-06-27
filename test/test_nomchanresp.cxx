#include "WireCellSigProc/NominalChannelResponse.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Testing.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TAxis.h"

#include <vector>
#include <iostream>

using namespace WireCell;
using WireCell::units::mV;
using WireCell::units::fC;
using WireCell::units::us;
using WireCell::SigProc::NominalChannelResponse;

using namespace std;

int main(int argc, char* argv[])
{
    NominalChannelResponse ncr;
    auto cfg = ncr.default_configuration();
    Binning binning(100,0,10*us);
    cfg["nbins"] = binning.nbins();
    cfg["tmin"] = binning.min();
    cfg["tmax"] = binning.max();

    const double GU = mV/fC;

    const std::vector<double> gains{7.8*GU, 14*GU};
    const std::vector<double> shapes{0.5*us,1.0*us,2.0*us,3.0*us};

    std::vector<TGraph*> graphs;

    for (auto gain : gains) {
        cfg["gain"]  = gain;
        for (auto shape : shapes) {
            cfg["shaping"] = shape;
            ncr.configure(cfg);
            auto wave = ncr.channel_response(0);
            const int nbins = wave.size();
            cerr << nbins << " " << binning.nbins() << endl;
            Assert(nbins == binning.nbins());
            TGraph* g = new TGraph(binning.nbins());
            for (int ind=0; ind<binning.nbins(); ++ind) {
                g->SetPoint(ind, binning.center(ind)/us, wave[ind]/GU);
            }
            graphs.push_back(g);
        }
    }

    const int colors[] = {1,2,4,6};
    TCanvas c("c","c",500,500);
    auto frame = c.DrawFrame(0.0, 0.0, 10.0, 15.0);
    frame->SetTitle("Nominal Channel Response various gains/shaping");
    frame->GetXaxis()->SetTitle("time [us]");
    frame->GetYaxis()->SetTitle("gain [mV/fC]");
    for (size_t ind=0; ind<graphs.size(); ++ind) {
        TGraph* g = graphs[ind];
        g->SetLineColor(colors[ind%4]);
        g->SetMarkerColor(colors[ind%4]);
        g->Draw("C*");
    }
    c.Print(Form("%s.pdf", argv[0]), "pdf");

    return 0;
}
