#include "WireCellSigProc/OmnibusSigProc.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Waveform.h"



#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellUtil/ExecMon.h"

#include <iostream>
#include <string>
#include <numeric>		// iota

#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"


using namespace WireCell;
using namespace std;


int ch2plane(int ch, int nuvw[])
{
    int nchans=0;
    for (int iplane = 0; iplane<3; ++iplane) {
        nchans += nuvw[iplane];
        if (ch < nchans) {
            return iplane;
        }
    }
    return -1;
}


void save_into_file(const char* input_filename, const char* output_filename,
                    IFrame::pointer frame_decon, int nrebin)
{
    TFile *tfile_in = new TFile(input_filename);

    TFile *tfile_out = new TFile(output_filename,"RECREATE");
    TTree *Trun = (TTree*)tfile_in->Get("Trun");
    if (Trun) {
        Trun = Trun->CloneTree();
        Trun->SetDirectory(tfile_out);
    }

    int nuvwt[4] = {0};

    TH1I* h_threshold[3] = {nullptr};
    for (int iplane=0; iplane<3; ++iplane) {
        std::vector<std::string> names{"orig", "raw", "baseline", "threshold"};
        for (auto name: names) {
            TH1 *hist = (TH1*)tfile_in->Get(Form("h%c_%s", 'u'+iplane, name.c_str()));
            if (hist) {
                hist->SetDirectory(tfile_out);

            }
            if (name == "raw") {
                if (!hist) {
                    std::cerr << "No raw data to steal from input file " << input_filename << std::endl;
                    exit(1);
                }
                nuvwt[iplane] = hist->GetNbinsX();
                nuvwt[3] = hist->GetNbinsY();
            }
            if (name == "threshold") {
                h_threshold[iplane] = (TH1I*)hist;
            }
        }
    }
  
    const int nticks = nuvwt[3];

    TH2F* h_decon[3] = {nullptr};

    int nbefore=0;
    for (int iplane=0; iplane<3; ++iplane) {
        if (!h_threshold[iplane]) { // in case input does not provide thresholds for us to cruelly overwrite
            const std::string name = Form("h%c_threshold", 'u'+iplane);
            h_threshold[iplane] = new TH1I(name.c_str(), name.c_str(),
                                           nuvwt[iplane], nbefore-0.5, nbefore+nuvwt[iplane]-0.5);
        }
        {
            const std::string name = Form("h%c_decon", 'u'+iplane);
            h_decon[iplane] = new TH2F(name.c_str(), name.c_str(),
                                       nuvwt[iplane], nbefore-0.5, nbefore+nuvwt[iplane]-0.5,
                                       nticks/nrebin, 0, nticks);
        }
        nbefore += nuvwt[iplane];
    }
  
    auto traces = frame_decon->traces();
    //for (auto trace : *traces.get()) {
    for (size_t index : frame_decon->tagged_traces("wiener")){
	auto trace = traces->at(index);
        int tbin = trace->tbin();
        int ch = trace->channel();
        auto charges = trace->charge();

        int iplane = ch2plane(ch, nuvwt);
        for (size_t ind=0; ind<charges.size(); ++ind) {
            h_decon[iplane]->Fill(ch, tbin + ind, charges[ind]);
        }
    }


  
    // save bad channels 
    TTree *T_bad = new TTree("T_bad","T_bad");
    int chid, plane, start_time,end_time;
    T_bad->Branch("chid",&chid,"chid/I");
    T_bad->Branch("plane",&plane,"plane/I");
    T_bad->Branch("start_time",&start_time,"start_time/I");
    T_bad->Branch("end_time",&end_time,"end_time/I");
    T_bad->SetDirectory(tfile_out);

    TTree *T_lf = new TTree("T_lf","T_lf");
    int channel;
    T_lf->Branch("channel",&channel,"channel/I");
  

    Waveform::ChannelMaskMap input_cmm = frame_decon->masks();
    for (auto const& it: input_cmm) {

        if (it.first == "bad"){ // save bad ... 
            //std::cout << "Xin1: " << it.first << " " << it.second.size() << std::endl;
            for (auto const &it1 : it.second){
                chid = it1.first;
                plane = ch2plane(chid, nuvwt);
                //std::cout << "Xin1: " << chid << " " << plane << " " << it1.second.size() << std::endl;
                for (size_t ind = 0; ind < it1.second.size(); ++ind){
                    start_time = it1.second[ind].first;
                    end_time = it1.second[ind].second;
                    T_bad->Fill();
                }
            }
        }
        else if (it.first =="lf_noisy"){
            for (auto const &it1 : it.second){
                channel = it1.first;
                T_lf->Fill();
            }
      
        }
        else if (it.first=="threshold"){
            for (auto const &it1 : it.second){
                chid = it1.first;
                float threshold = it1.second[0].first/it1.second[0].second;

                const int iplane = ch2plane(chid, nuvwt);
                int nbefore = 0;
                for (int ind=0; ind<iplane; ++ind) { nbefore += nuvwt[iplane]; }
                h_threshold[iplane]->SetBinContent(chid+1-nbefore, threshold*nrebin*3);
            }
        }
    }

    tfile_out->Write();
    tfile_out->Close();
  
}


class XinFileIterator {
    TH2* hist[3];		// per plane
    WireCell::Waveform::ChannelMaskMap ret;
    TFile *file;
public:
    XinFileIterator(const char* filename, const char* histtype="raw") {
        file = TFile::Open(filename);
        string uvw = "uvw";
        for (int ind=0; ind<3; ++ind) {
            auto c = uvw[ind];
            std::string name = Form("h%c_%s", c, histtype);
            cerr << "Loading " << name << endl;
            hist[ind] = (TH2*)file->Get(name.c_str());
        }
      
        TTree *T_bad = (TTree*)file->Get("T_bad");
        int chid, plane, start_time,end_time;
        T_bad->SetBranchAddress("chid",&chid);
        T_bad->SetBranchAddress("plane",&plane);
        T_bad->SetBranchAddress("start_time",&start_time);
        T_bad->SetBranchAddress("end_time",&end_time);
      
        for (int i=0;i!=T_bad->GetEntries();i++){
            T_bad->GetEntry(i);
            WireCell::Waveform::BinRange chirped_bins;
            chirped_bins.first = start_time;
            chirped_bins.second = end_time;
            ret["bad"][chid].push_back(chirped_bins);
        }
      
      
        TTree *T_lf = (TTree*)file->Get("T_lf");
        int channel;
        T_lf->SetBranchAddress("channel",&channel);
        for (int i=0;i!=T_lf->GetEntries();i++){
            T_lf->GetEntry(i);
            WireCell::Waveform::BinRange chirped_bins;
            chirped_bins.first = 0;
            chirped_bins.second = hist[0]->GetNbinsY();
            ret["lf_noisy"][channel].push_back(chirped_bins);
        }
        delete T_lf;
        delete T_bad;
      
	//file->Close();
	//delete file;
    }

    int plane(int ch) {
	if (ch < 2400) return 0;
	if (ch < 2400+2400) return 1;
	return 2;
    }
    int index(int ch) {
	if (ch < 2400) return ch;
	if (ch < 2400+2400) return ch-2400;
	return ch-2400-2400;
    }

    vector<float> at(int ch) {
	TH2* h = hist[plane(ch)];
	int ind = index(ch);
	vector<float> ret(9600);
	for (int itick=0; itick<9600; ++itick) {
	    ret[itick] = h->GetBinContent(ind+1, itick+1);
	}
	return ret;
    }

    void clear(){
        delete hist[0];
        delete hist[1];
        delete hist[2];
    
        file->Close();
        delete file;
    }
  
    /// Return a frame, the one and only in the file.
    IFrame::pointer frame() {
	ITrace::vector traces;

	int chindex=0;
	for (int iplane=0; iplane<3; ++iplane) {
	    TH2* h = hist[iplane];

	    int nchannels = h->GetNbinsX();
	    int nticks = h->GetNbinsY();

	    cerr << "plane " << iplane << ": " << nchannels << " X " << nticks << endl;

	    double qtot = 0.0;
	    for (int ich = 1; ich <= nchannels; ++ich) {
		ITrace::ChargeSequence charges;
		for (int itick = 1; itick <= nticks; ++itick) {
		    auto q = h->GetBinContent(ich, itick);
		    charges.push_back(q);
		    qtot += q;
		}
		SimpleTrace* st = new SimpleTrace(chindex, 0.0, charges);
		traces.push_back(ITrace::pointer(st));
		++chindex;
		//cerr << "qtot in plane/ch/index "
		//     << iplane << "/" << ich << "/" << chindex << " = " << qtot << endl;
	    }
	}
	SimpleFrame* sf = new SimpleFrame(0, 0, traces, 0.5*units::microsecond, ret);
	return IFrame::pointer(sf);
    }
    
};




int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "This test needs an input data file in \"Magnify\" ROOT format with \"raw\" histograms." << endl;
	return 1;
    }
    std::string url_in = argv[1];
    std::string url_out = Form("%s.root", argv[0]);
    if (argc > 2) {
      url_out = argv[2];
    }
    cerr << "Output too: " << url_out << std::endl;

    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");
    pm.add("WireCellSigProc");

    string filenames[4] = {
        "microboone-noise-spectra-v2.json.bz2",
        "garfield-1d-3planes-21wires-6impacts-v6.json.bz2",
        "microboone-celltree-wires-v2.json.bz2",
        "ub-10-wnormed.json.bz2",
    };
    
    // do the geometry ... 
    {
        auto anodecfg = Factory::lookup<IConfigurable>("AnodePlane");
        auto cfg = anodecfg->default_configuration();
        cfg["fields"] = filenames[3];
        cfg["wires"] = filenames[2];
        anodecfg->configure(cfg);
    }

    // add the response function ...
    {
        auto ifrcfg = Factory::lookup<IConfigurable>("FieldResponse");
        auto cfg = ifrcfg->default_configuration();
        cfg["filename"] = filenames[3];
        ifrcfg->configure(cfg);
    }

    Int_t nticks = 9592;
    
    // add the filters
    {
        // Tight Gaussian filters
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Gaus_tight");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 1.11408e-01 * units::megahertz;
            cfg["power"] = 2;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }

        // Tight Wiener filters for U for ROI finding
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wiener_tight_U");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 5.75416e+01/800.*2 * units::megahertz;
            cfg["power"] = 4.10358e+00;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }
      
        // Tight Wiener filters for V for ROI finding
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wiener_tight_V");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 5.99306e+01/800.*2* units::megahertz;
            cfg["power"] = 4.20820e+00;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }
      
        // Tight Wiener filters for W 
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wiener_tight_W");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 5.88802e+01/800.*2  * units::megahertz;
            cfg["power"] = 4.17455e+00;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }
      
      
        // Wide Wiener filters for U for hit
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wiener_wide_U");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 1.78695e+01/200.*2.  * units::megahertz;
            cfg["power"] = 5.33129e+00;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }

      
        // Wide Wiener filters for V for hit
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wiener_wide_V");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 1.84666e+01/200.*2.  * units::megahertz;
            cfg["power"] = 5.60489e+00;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }
      
        // Wide Wiener filters for W for hit
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wiener_wide_W");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 1.83044e+01/200.*2. * units::megahertz;
            cfg["power"] = 5.44945e+00;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }
      
        // Wide Gaussian filters for charge
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Gaus_wide");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["sigma"] = 0.14 * units::megahertz;
            cfg["power"] = 2;
            cfg["flag"] = true;
            incrcfg->configure(cfg);
        }
      
        // Tight low frequency filter for ROI
        {
            auto incrcfg = Factory::lookup<IConfigurable>("LfFilter","ROI_tight_lf");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["tau"] = 0.02 * units::megahertz;
            incrcfg->configure(cfg);
        }

        {
            auto incrcfg = Factory::lookup<IConfigurable>("LfFilter","ROI_tighter_lf");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["tau"] = 0.1 * units::megahertz;
            incrcfg->configure(cfg);
        }
      
        // Loose low frequency filter for ROI
        {
            auto incrcfg = Factory::lookup<IConfigurable>("LfFilter","ROI_loose_lf");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = nticks;
            cfg["max_freq"] = 1 * units::megahertz;
            cfg["tau"] = 0.0025 * units::megahertz;
            incrcfg->configure(cfg);
        }
      
        // Wire Filter for induction planes
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wire_ind");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = 2400;
            cfg["max_freq"] = 1 ;
            cfg["sigma"] = 1./sqrt(3.1415926)*1.4 ;
            cfg["power"] = 2;
            cfg["flag"] = false;
            incrcfg->configure(cfg);
        }
      
        // Wire Filter for collection planes 
        {
            auto incrcfg = Factory::lookup<IConfigurable>("HfFilter","Wire_col");
            auto cfg = incrcfg->default_configuration();
            cfg["nbins"] = 3456;
            cfg["max_freq"] = 1 ;
            cfg["sigma"] = 1.0/sqrt(3.1415926)*3.0 ;
            cfg["power"] = 2;
            cfg["flag"] = false;
            incrcfg->configure(cfg);
        }
    }

    // per channel response 
    {
        const std::string cr_tn = "PerChannelResponse";
        const std::string pcr_filename = "calib_resp_v2.json.bz2";
        auto icrcfg = Factory::lookup<IConfigurable>(cr_tn);
        auto cfg = icrcfg->default_configuration();
        cfg["filename"] = pcr_filename;
        icrcfg->configure(cfg);
    }
    

    // std::cout << "asd " << std::endl;
    


    // various software filters ...
    // one HF for ROI,  two LF for ROI finding
    // one HF for charge
    // one HF for hit
    // two HF filters for wire dimension 
    
    
    int nrebin = 4;


    ExecMon em("starting");

    XinFileIterator fs(url_in.c_str());

    cerr <<  em("loading rootfiles") << endl;

    IFrame::pointer frame = fs.frame();

    //    cerr << em("fill the frame") << endl;

    fs.clear();
    
    cerr <<  em("close the file") << endl;
    
    
    SigProc::OmnibusSigProc bus;
    bus.configure(bus.default_configuration());


    IFrame::pointer frame_decon;
    
    cerr << em("Do deconvolution") << endl;
    bus(frame, frame_decon);
    cerr << em(" ... done") << endl;
    
    save_into_file(url_in.c_str(), url_out.c_str(), frame_decon, nrebin);
    
    
}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
