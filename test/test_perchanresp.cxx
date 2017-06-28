#include "WireCellUtil/Units.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Exceptions.h"

/// needed to pretend like we are doing WCT internals
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IChannelResponse.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"


#include "TCanvas.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TAxis.h"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

using namespace WireCell;
using WireCell::units::mV;
using WireCell::units::fC;
using WireCell::units::us;

using namespace std;

int main(int argc, char* argv[])
{
    /// Real component code would get this info from its
    /// configuration.
    const std::string cr_tn = "PerChannelResponse";
    const std::string ap_tn = "AnodePlane";
    const std::string pcr_filename = "calib_resp_v1.json.bz2";
    const std::string fr_filename = "ub-10-wnormed.json.bz2";
    const std::string wires_filename = "microboone-celltree-wires-v2.1.json.bz2";


    /// User code should never do this.
    try {                           
        PluginManager& pm = PluginManager::instance();
        pm.add("WireCellSigProc");
        pm.add("WireCellGen");

        auto icrcfg = Factory::lookup<IConfigurable>(cr_tn);
        auto cfg = icrcfg->default_configuration();
        cfg["filename"] = pcr_filename;
        icrcfg->configure(cfg);

        auto iapcfg = Factory::lookup<IConfigurable>(ap_tn);
        cfg = iapcfg->default_configuration();
        cfg["fields"] = fr_filename;
        cfg["wires"] = wires_filename;
        iapcfg->configure(cfg);
    }
    catch (Exception& e) {
        cerr << "caught Exception: " << errstr(e) << endl;
        return 1;
    }


    /// Finally, we now pretend to be real component code.
    auto cr = Factory::find<IChannelResponse>(cr_tn);
    auto ap = Factory::find<IAnodePlane>(ap_tn);

    std::vector<int> planechans[3]; // fixme: will break with DUNE
    for (auto ch : ap->channels()) {
        auto wpid = ap->resolve(ch);
        planechans[wpid.index()].push_back(ch);
    }

    gStyle->SetOptStat(0);
    TCanvas c("c","c",500,500);
    c.Divide(3,1);

    // Desired gain units for showing in the plot
    const double GU = units::mV/units::fC;

    for (int iplane=0; iplane<3; ++iplane) {
        auto& channels = planechans[iplane];
        std::sort(channels.begin(), channels.end());

        /// assume all responses in a plane are the same size.
        const int nsamps = cr->channel_response(channels[0]).size();
        const int nchans = channels.size();

        Assert(nsamps>0);
        Assert(nchans>0);


        TH2F* hist = new TH2F(Form("hist%d", iplane),
                              Form("Per Channel Response Plane %d [mV/fC]", iplane),
                              nsamps, 0, nsamps,
                              nchans, 0, nchans);
        hist->GetXaxis()->SetTitle("ticks");
        hist->GetYaxis()->SetTitle("channel indices");

        for (int ich=0; ich<nchans; ++ich) {
            const auto& resp = cr->channel_response(channels[ich]);
            for (int isamp=0; isamp<nsamps; ++isamp) {
                hist->Fill(isamp+0.5, ich+0.5, resp[isamp]/GU);
            }
        }

        c.cd(iplane+1);
        hist->Draw("colz");
    }

    cerr << "Now ROOT makes the PDF" << endl;
    c.Print(Form("%s.pdf", argv[0]), "pdf");

    return 0;
}
