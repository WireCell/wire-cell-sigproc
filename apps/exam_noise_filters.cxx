
#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellSigProc/OmnibusNoiseFilter.h"
#include "WireCellSigProc/OneChannelNoise.h"
#include "WireCellSigProc/CoherentNoiseSub.h"
#include "WireCellSigProc/SimpleChannelNoiseDB.h"

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/ExecMon.h"

#include <iostream>
#include <string>
#include <numeric>		// iota
#include <string>

#include "TCanvas.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"

using namespace WireCell;
using namespace std;

const string url_test = "/data0/bviren/data/uboone/test_3455_0.root"; // big!

void save_into_file(const char* filename,IFrame::pointer frame_orig,IFrame::pointer frame_raw, int nwire_u, int nwire_v, int nwire_w, int nticks){
  TFile *file1 = new TFile(filename);

  TFile *file = new TFile("temp.root","RECREATE");
  TTree *Trun = ((TTree*)file1->Get("Trun"))->CloneTree();
  Trun->SetDirectory(file);
  
  TH2I *hu_orig = new TH2I("hu_orig","hu_orig",nwire_u,-0.5,nwire_u-0.5,nticks,0,nticks);
  TH2I *hv_orig = new TH2I("hv_orig","hv_orig",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,nticks,0,nticks);
  TH2I *hw_orig = new TH2I("hw_orig","hw_orig",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,nticks,0,nticks);

  TH2F *hu_raw = new TH2F("hu_raw","hu_raw",nwire_u,-0.5,nwire_u-0.5,nticks,0,nticks);
  TH2F *hv_raw = new TH2F("hv_raw","hv_raw",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,nticks,0,nticks);
  TH2F *hw_raw = new TH2F("hw_raw","hw_raw",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,nticks,0,nticks);

  TH2F *hu_decon = new TH2F("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,int(nticks/6.),0,nticks);
  TH2F *hv_decon = new TH2F("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,int(nticks/6.),0,nticks);
  TH2F *hw_decon = new TH2F("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,int(nticks/6.),0,nticks);
  
  TH1F *hu_baseline = (TH1F*)file1->Get("hu_baseline");
  TH1F *hv_baseline = (TH1F*)file1->Get("hv_baseline");
  TH1F *hw_baseline = (TH1F*)file1->Get("hw_baseline");
  
  hu_baseline->SetDirectory(file);
  hv_baseline->SetDirectory(file);
  hw_baseline->SetDirectory(file);

  TH1F *hu_threshold = (TH1F*)file1->Get("hu_threshold");
  TH1F *hv_threshold = (TH1F*)file1->Get("hv_threshold");
  TH1F *hw_threshold = (TH1F*)file1->Get("hw_threshold");
  
  hu_threshold->SetDirectory(file);
  hv_threshold->SetDirectory(file);
  hw_threshold->SetDirectory(file);


  auto traces = frame_orig->traces();
  for (auto trace : *traces.get()) {
    int tbin = trace->tbin();
    int ch = trace->channel();
    auto charges = trace->charge();
    if (ch < nwire_u){
      int counter = 0;
      for (auto q : charges) {
	counter ++;
	hu_orig->SetBinContent(ch+1,tbin+counter,q); 
      }
    }else if (ch < nwire_v + nwire_u){
      int counter = 0;
      for (auto q : charges) {
	counter ++;
	hv_orig->SetBinContent(ch+1-nwire_u,tbin+counter,q); 
      }
    }else{
      int counter = 0;
      for (auto q : charges) {
	counter ++;
	hw_orig->SetBinContent(ch+1-nwire_u-nwire_v,tbin+counter,q); 
      }
    }
  }

  traces = frame_raw->traces();
  for (auto trace : *traces.get()) {
    int tbin = trace->tbin();
    int ch = trace->channel();
    auto charges = trace->charge();
    if (ch < nwire_u){
      int counter = 0;
      for (auto q : charges) {
	counter ++;
	hu_raw->SetBinContent(ch+1,tbin+counter,q); 
      }
    }else if (ch < nwire_v + nwire_u){
      int counter = 0;
      for (auto q : charges) {
	counter ++;
	hv_raw->SetBinContent(ch+1-nwire_u,tbin+counter,q); 
      }
    }else{
      int counter = 0;
      for (auto q : charges) {
	counter ++;
	hw_raw->SetBinContent(ch+1-nwire_u-nwire_v,tbin+counter,q); 
      }
    }
  }

  // save bad channels 
  TTree *T_bad = new TTree("T_bad","T_bad");
  int chid, plane, start_time,end_time;
  T_bad->Branch("chid",&chid,"chid/I");
  T_bad->Branch("plane",&plane,"plane/I");
  T_bad->Branch("start_time",&start_time,"start_time/I");
  T_bad->Branch("end_time",&end_time,"end_time/I");
  T_bad->SetDirectory(file);

  Waveform::ChannelMaskMap input_cmm = frame_raw->masks();
  for (auto const& it: input_cmm) {
    //std::cout << "Xin1: " << it.first << " " << it.second.size() << std::endl;
    for (auto const &it1 : it.second){
      chid = it1.first;
      if (chid < nwire_u){
	plane = 0;
      }else if (chid < nwire_v){
	plane = 1;
      }else{
	plane = 2;
      }
      //std::cout << "Xin1: " << chid << " " << plane << " " << it1.second.size() << std::endl;
      for (int ind = 0; ind < it1.second.size(); ++ind){
	start_time = it1.second[ind].first;
	end_time = it1.second[ind].second;
	T_bad->Fill();
      }
    }
  }


  file->Write();
  file->Close();
}

void rms_plot(TCanvas& canvas, IFrame::pointer frame, const string& title)
{
    //TProfile h("h", title.c_str(), 9600, 0, 9600);
    TH2F h("h", title.c_str(), 9600, 0, 9600,100,0,1000);

    cerr << title << endl;

    auto traces = frame->traces();
    for (auto trace : *traces.get()) {
	int tbin = trace->tbin();
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

class XinFileIterator {
    TH2* hist[3];		// per plane
public:
    XinFileIterator(const char* filename, const char* histtype="orig") {
	TFile* file = TFile::Open(filename);
	string uvw = "uvw";
	for (int ind=0; ind<3; ++ind) {
	    auto c = uvw[ind];
	    std::string name = Form("h%c_%s", c, histtype);
	    cerr << "Loading " << name << endl;
	    hist[ind] = (TH2*)file->Get(name.c_str());
	}
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
	SimpleFrame* sf = new SimpleFrame(0, 0, traces);
	return IFrame::pointer(sf);
    }
    
};


int main(int argc, char* argv[])
{
    string url = url_test;
    if (argc > 1) {
	url = argv[1];
    }
    XinFileIterator fs(url.c_str());

    // S&C microboone sampling parameter database
    const double tick = 0.5*units::microsecond;
    const int nsamples = 9594;

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
    for (int ind=2016; ind<= 2095; ++ind) { miscfgchan.push_back(ind); }
    for (int ind=2192; ind<= 2303; ++ind) { miscfgchan.push_back(ind); }
    for (int ind=2352; ind< 2400; ++ind) { miscfgchan.push_back(ind); }
    
    // hard-coded bad channels
    vector<int> bad_channels;
    for (int i=0;i!=wchans.size();i++){
      if (i>=7136 - 4800 && i <=7263 - 4800){
	if (i != 7200- 4800 && i!=7215 - 4800)
	  bad_channels.push_back(i+4800);
      }
    }

    // Q&D RC+RC time constant - all have same.
    const double rcrc = 1.0*units::millisecond;
    vector<int> rcrcchans(nchans);
    std::iota(rcrcchans.begin(), rcrcchans.end(), 0);
    
    //harmonic noises
    vector<int> harmonicchans(uchans.size() + vchans.size());
    std::iota(harmonicchans.begin(), harmonicchans.end(), 0);
    
    vector<int> special_chans;
    special_chans.push_back(2240);

    WireCellSigProc::SimpleChannelNoiseDB::mask_t h36kHz(0,169,173);
    WireCellSigProc::SimpleChannelNoiseDB::mask_t h108kHz(0,513,516);
    WireCellSigProc::SimpleChannelNoiseDB::mask_t hspkHz(0,17,19);
    WireCellSigProc::SimpleChannelNoiseDB::multimask_t hharmonic;
    hharmonic.push_back(h36kHz);
    hharmonic.push_back(h108kHz);
    WireCellSigProc::SimpleChannelNoiseDB::multimask_t hspecial;
    hspecial.push_back(h36kHz);
    hspecial.push_back(h108kHz);
    hspecial.push_back(hspkHz);

    // do the coherent subtraction
    
    std::vector< std::vector<int> > channel_groups;
    for (int i=0;i!=172;i++){
      std::vector<int> channel_group;
      for (int j=0;j!=48;j++){
	channel_group.push_back(i*48+j);
      }
      channel_groups.push_back(channel_group);
    }

    // Load up components.  Note, in a real app this is done as part
    // of factory + configurable and driven by user configuration.

    auto noise = new WireCellSigProc::SimpleChannelNoiseDB;
    // initialize
    noise->set_sampling(tick, nsamples);
    // set nominal baseline
    noise->set_nominal_baseline(uchans, unombl);
    noise->set_nominal_baseline(vchans, vnombl);
    noise->set_nominal_baseline(wchans, wnombl);
    // set misconfigured channels
    noise->set_gains_shapings(miscfgchan, from_gain_mVfC, to_gain_mVfC, from_shaping, to_shaping);
    // do the RCRC
    noise->set_rcrc_constant(rcrcchans, rcrc);
    // set initial bad channels
    noise->set_bad_channels(bad_channels);
    // set the harmonic filter
    noise->set_filter(harmonicchans,hharmonic);
    noise->set_filter(special_chans,hspecial);
    noise->set_channel_groups(channel_groups);

    shared_ptr<WireCell::IChannelNoiseDatabase> noise_sp(noise);

    auto one = new WireCellSigProc::OneChannelNoise;
    one->set_channel_noisedb(noise_sp);
    shared_ptr<WireCell::IChannelFilter> one_sp(one);

    auto many = new WireCellSigProc::CoherentNoiseSub;
    shared_ptr<WireCell::IChannelFilter> many_sp(many);


    WireCellSigProc::OmnibusNoiseFilter bus;
    bus.set_channel_filters({one_sp});
    bus.set_grouped_filters({many_sp});
    bus.set_channel_noisedb(noise_sp);

    //TCanvas canvas("c","canvas",500,500);

    //canvas.Print("test_omnibus.pdf[","pdf");

    ExecMon em("starting");

    // This might be done in a DFP graph in a real app 
    IFrame::pointer frame = fs.frame();
    // rms_plot(canvas, frame, "Raw frame");
	
    IFrame::pointer quiet;

    cerr << em("Removing noise") << endl;
    bus(frame, quiet);
    cerr << em("...done") << endl;

    // rms_plot(canvas, quiet, "Quiet frame");
    Assert(quiet);

    save_into_file(url.c_str(),frame,quiet,uchans.size(),vchans.size(),wchans.size(),nsamples);
    //    canvas.Print("test_omnibus.pdf]","pdf");

    cerr << em.summary() << endl;   

    return 0;
}
