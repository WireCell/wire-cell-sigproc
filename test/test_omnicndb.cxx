#include "WireCellSigProc/OmniChannelNoiseDB.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TH1F.h"

const std::string config_text = R"JSONNET(
// example OmniChannelNoiseDB configuration.

local wc = import "wirecell.jsonnet";
local anodes = import "multi/anodes.jsonnet"; 
{
    anode: wc.tn(anodes.nominal),  
    tick: 0.5*wc.us,
    nsamples: 9600,

    // channel groups is a 2D list.  Each element is one group of
    // channels which should be considered together for coherent noise
    // filtering.  
    groups: [std.range(g*48, (g+1)*48-1) for g in std.range(0,171)],

    // Default channel info which is used if not overriden by one of
    // the channel_info objects (below) for a given channel.
    default_info : {
	nominal_baseline: 0.0,  // adc count
        gain_correction: 1.0,     // unitless
        response_offset: 79,      // ticks?
        pad_window_front: 20,     // ticks?
	pad_window_back: 10,      // ticks?
	min_rms_cut: 1.0,         // units???
	max_rms_cut: 5.0,         // units???

        // parameter used to make "rcrc" spectrum
        rcrc: 1.0*wc.millisecond,

        // parameters used to make "config" spectrum
        reconfig : {},

        // list to make "noise" spectrum mask
        freqmasks: [],

        // field response waveform to make "response" spectrum
        response: {},

    },

    // overide defaults for specific channels.  If an info is
    // mentioned for a particular channel in multiple objects, last
    // one wins.
    channel_info: [             

        {
            channels: 4,        // single channel
            nominal_baseline: 400,
        },
        {
            channels: [1,42,69],      // explicit list
	    nominal_baseline: 2048.0, // adc count
            reconfig : {
                from: {gain: 7.8*wc.mV/wc.fC, shaping: 1.0*wc.us},
                to: {gain: 14.0*wc.mV/wc.fC, shaping: 2.0*wc.us}, 
            }
        },
        {
            channels: { first: 1, last: 4 }, // inclusive range.
            // could also use Jsonnet's std.range(first,last):
            // channels: std.range(1,4),
	    nominal_baseline: 400.0,
        },
        {
            channels: {first: 0, last: 2047,},
            freqmasks: [            // masks in frequency domain.
                { value: 1.0, lobin: 0, hibin: $.nsamples-1 },
                { value: 0.0, lobin: 169, hibin: 173 },
                { value: 0.0, lobin: 513, hibin: 516 },
            ]
        },
        {
            // All channels in a given anode wire plane.  Note, if
            // used, the anode plane must be added to top level
            // configuration sequence.  The "channels" can be
            // specified by the other means.
            channels: { wpid: wc.WirePlaneId(wc.Ulayer) },
            response: {
                // the average response to use that of the wpid.
                wpid: wc.WirePlaneId(wc.Ulayer)
                // or as a raw waveform:
                // waveform: [time-domain samples assumed to be at tick sampling period]
                // waveformid: <uniquenumber>
            }
        }
    ],

}

)JSONNET";




using namespace WireCell;
using namespace std;


typedef std::unordered_set<const WireCell::IChannelNoiseDatabase::filter_t*> filter_bag_t;


void plot_spec(const filter_bag_t& specs, const std::string& name)
{

    std::vector<TGraph*> graphs;
    std::vector<float> tmp;
    for (auto spec: specs) {
        TGraph* graph = new TGraph();
        graphs.push_back(graph);
        for (size_t ind=0; ind<spec->size(); ++ind) {
            double amp = std::abs(spec->at(ind));
            graph->SetPoint(ind, ind, amp);
            tmp.push_back(amp);
        }
    }
    auto mme = std::minmax_element(tmp.begin(), tmp.end());
    float ymin = *mme.first;
    float ymax = *mme.second;

    const int ncolors=5;
    int colors[ncolors] = {1,2,4,6,8};
    for (size_t igraph = 0; igraph<graphs.size(); ++igraph) {
        TGraph* graph = graphs[igraph];
        graph->SetLineColor(colors[igraph%ncolors]);
        graph->SetLineWidth(2);

        if (!igraph) {
            auto frame = graph->GetHistogram();
            frame->SetTitle(name.c_str());
            frame->GetXaxis()->SetTitle("frequency bins");
            frame->GetYaxis()->SetTitle("amplitude");
            graph->Draw("AL");
            frame->SetMinimum(ymin);
            frame->SetMaximum(ymax);
            cerr << name << " ymin=" << ymin << ", ymax=" << ymax << endl;
        }
        else {
            graph->Draw("L");
        }
    }
}


int main(int argc, char* argv[])
{
    /// User code should never do this.
    const std::string pcr_filename = "calib_resp_v1.json.bz2";
    const std::string fr_filename = "ub-10-wnormed.json.bz2";
    const std::string wires_filename = "microboone-celltree-wires-v2.1.json.bz2";
    try {
        auto& pm = PluginManager::instance();
        pm.add("WireCellSigProc");
        pm.add("WireCellGen");

        auto iapcfg = Factory::lookup<IConfigurable>("AnodePlane");
        auto cfg = iapcfg->default_configuration();
        cfg["fields"] = fr_filename;
        cfg["wires"] = wires_filename;
        iapcfg->configure(cfg);
    }
    catch (Exception& e) {
        cerr << "caught Exception: " << errstr(e) << endl;
        return 1;
    }


    Configuration cfg;
    if (argc > 1) {
        cerr << "testing with " << argv[1] << endl;
        cfg = Persist::load(argv[1]);
    }
    else {
        cerr << "testing with build in config text\n";
        cfg = Persist::loads(config_text);
    }

    
    SigProc::OmniChannelNoiseDB db;
    auto def = db.default_configuration();
    cfg = update(def, cfg);
    db.configure(cfg);

    
    auto anode = Factory::find<IAnodePlane>("AnodePlane");
    int nchannels = anode->channels().size();

    gStyle->SetOptStat(0);
    TCanvas canvas("canvas","canvas",500,500);
    canvas.Print(Form("%s.pdf[",argv[0]),"pdf");
    canvas.SetGridx(1);
    canvas.SetGridy(1);
    int nticks = db.number_samples();
    double tick = db.sample_time();
    cerr << nticks << " at " << tick/units::us << " us.\n";

    std::vector<std::string> scalar_names{
        "nominal baseline", "gain correction", "response offset", "pad window front", "pad window back",
            "min rms cut", "max rms cut"};
            

    std::vector<TGraph> scalars(7);
    for (int ch=0; ch<nchannels; ++ch) {
        scalars[0].SetPoint(ch, ch, db.nominal_baseline(ch));
        scalars[1].SetPoint(ch, ch, db.gain_correction(ch));
        scalars[2].SetPoint(ch, ch, db.response_offset(ch));
        scalars[3].SetPoint(ch, ch, db.pad_window_front(ch));
        scalars[4].SetPoint(ch, ch, db.pad_window_back(ch));
        scalars[5].SetPoint(ch, ch, db.min_rms_cut(ch));
        scalars[6].SetPoint(ch, ch, db.max_rms_cut(ch));
    }

    for (size_t ind=0; ind<scalars.size(); ++ind) {
        auto& graph = scalars[ind];
        graph.SetName(scalar_names[ind].c_str());
        graph.SetLineColor(2);
        graph.SetLineWidth(3);
        graph.Draw("AL");

        auto frame = graph.GetHistogram();
        frame->SetTitle(scalar_names[ind].c_str());
        frame->GetXaxis()->SetTitle("channels");
        canvas.Print(Form("%s.pdf",argv[0]),"pdf");
    }


    // find all unique
    std::vector<filter_bag_t> specs(4);

    for (int ch=0; ch<nchannels; ++ch) {
        auto const& f0 = db.rcrc(ch);
        if (specs[0].find(&f0) == specs[0].end()) {
            specs[0].insert(&f0);
        }
        auto const& f1 = db.config(ch);
        if (specs[1].find(&f1) == specs[1].end()) {
            specs[1].insert(&f1);
        }
        auto const& f2 = db.noise(ch);
        if (specs[2].find(&f2) == specs[2].end()) {
            specs[2].insert(&f2);
        }
        auto const& f3 = db.response(ch);
        if (specs[3].find(&f3) == specs[3].end()) {
            specs[3].insert(&f3);
        }
    }

    const std::string spec_names[4] = { "RCRC", "Misconfig", "Noise Mask", "Response" };

    for (int ind=0; ind<4; ++ind) {
        cerr << spec_names[ind] << ": " << specs[ind].size() << " spectra\n";
        plot_spec(specs[ind], spec_names[ind]);
        canvas.Print(Form("%s.pdf",argv[0]),"pdf");
    }



    canvas.Print(Form("%s.pdf]",argv[0]),"pdf");

    
    return 0;
}
