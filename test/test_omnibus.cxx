/**
   This the OmnibusNoiseFilter class outside of the context of
   the DFP, Factory and Configurable patterns.  As such it is an
   example of tightly coupling exact input through the SST FrameSource
   to the exact OmnibusNoiseFilter class

 */

#include "WireCellSigProc/OmnibusNoiseFilter.h"
#include "WireCellSigProc/Microboone.h"

#include "WireCellSigProc/SimpleChannelNoiseDB.h"
#include "WireCellSst/FrameSource.h"

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/ExecMon.h"

#include <iostream>
#include <string>
#include <numeric>		// iota
#include <string>

#include "TCanvas.h"
#include "TProfile.h"
#include "TH2F.h"

using namespace WireCell;
using namespace std;

const string url_alpha = "http://www.phy.bnl.gov/wire-cell/examples/celltree/3.0/celltree_alpha-new.root";
const string url_cosmic = "http://www.phy.bnl.gov/wire-cell/examples/celltree/3.0/celltree-cosmic-new.root";

void rms_plot(TCanvas& canvas, IFrame::pointer frame, const string& title)
{
    //TProfile h("h", title.c_str(), 9600, 0, 9600);
    TH2F h("h", title.c_str(), 9600, 0, 9600,100,0,1000);

    cerr << title << endl;

    auto traces = frame->traces();
    for (auto trace : *traces.get()) {
	//int tbin = trace->tbin();
	int ch = trace->channel();
	auto charges = trace->charge();
	//cerr << "ch:" << ch <<", tbin:" << tbin <<", " << charges.size() << " charges\n";
	for (auto q : charges) {
	    h.Fill(ch, q);
	}
    }

    //h.Draw();
    h.Draw("colz");
    canvas.Print("test_omnibus.pdf","pdf");
}


int main(int argc, char* argv[])
{
    string url = url_alpha;	// small file
    if (argc > 1) {
	url = argv[1];
    }

    WireCellSst::FrameSource fs;
    auto cfg = fs.default_configuration();
    put(cfg, "filename", url);
    fs.configure(cfg);

    // S&C microboone sampling parameter database
    const double tick = 0.5*units::microsecond;
    const int nsamples = 9600;

    // Q&D microboone channel map
    vector<int> uchans(2400), vchans(2400), wchans(3456);
    const int nchans = uchans.size() + vchans.size() + wchans.size();
    std::iota(uchans.begin(), uchans.end(), 0);
    std::iota(vchans.begin(), vchans.end(), vchans.size());
    std::iota(wchans.begin(), wchans.end(), vchans.size() + uchans.size());


    // Q&D nominal baseline
    const double unombl=2048.0, vnombl=2048.0, wnombl=400.0;

    // Q&D miss-configured channel database
    vector<int> miscfgchan;
    const double from_gain_mVfC=7.8, to_gain_mVfC=14.0,
	from_shaping=1.0*units::microsecond, to_shaping=2.0*units::microsecond;
    for (int ind=2016; ind<= 2096; ++ind) { miscfgchan.push_back(ind); }
    for (int ind=2192; ind<= 2303; ++ind) { miscfgchan.push_back(ind); }
    for (int ind=2352; ind<= 2400; ++ind) { miscfgchan.push_back(ind); }
    
    // Q&D RC+RC time constant - all have same.
    const double rcrc = 1.0*units::millisecond;
    vector<int> rcrcchans(nchans);
    std::iota(rcrcchans.begin(), rcrcchans.end(), 0);

    // Load up components.  Note, in a real app this is done as part
    // of factory + configurable and driven by user configuration.

    auto noise = new SigProc::SimpleChannelNoiseDB;
    noise->set_nominal_baseline(uchans, unombl);
    noise->set_nominal_baseline(vchans, vnombl);
    noise->set_nominal_baseline(wchans, wnombl);
    noise->set_gains_shapings(miscfgchan, from_gain_mVfC, to_gain_mVfC, from_shaping, to_shaping);
    noise->set_sampling(tick, nsamples);
    noise->set_rcrc_constant(rcrcchans, rcrc);
    shared_ptr<WireCell::IChannelNoiseDatabase> noise_sp(noise);

    auto one = new SigProc::Microboone::OneChannelNoise;
    one->set_channel_noisedb(noise_sp);
    shared_ptr<WireCell::IChannelFilter> one_sp(one);

    auto many = new SigProc::Microboone::CoherentNoiseSub;
    shared_ptr<WireCell::IChannelFilter> many_sp(many);

    SigProc::OmnibusNoiseFilter bus;
    bus.set_channel_filters({one_sp});
    bus.set_grouped_filters({many_sp});
    bus.set_channel_noisedb(noise_sp);

    TCanvas canvas("c","canvas",500,500);

    canvas.Print("test_omnibus.pdf[","pdf");

    ExecMon em("starting");

    // This might be done in a DFP graph in a real app 
    IFrame::pointer frame;
    int count = 0;
    while (fs(frame)) {
	++count;
	if (!frame) {
	    cerr << "Hist end of stream, bye." << endl;
	    break;
	}
	rms_plot(canvas, frame, Form("Raw frame #%d", count));
	
	IFrame::pointer quiet;

	cerr << em(Form("Removing noise from %d....", count)) << endl;
	bus(frame, quiet);
	cerr << em("...done") << endl;

	rms_plot(canvas, quiet, Form("Quiet frame #%d", count));
	Assert(quiet);
    }
    canvas.Print("test_omnibus.pdf]","pdf");

    cerr << em.summary() << endl;   

    return 0;
}
