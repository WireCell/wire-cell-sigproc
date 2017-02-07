
#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellSigProc/OmnibusNoiseFilter.h"
#include "WireCellSigProc/Microboone.h"

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
    if (argc < 2) {
	cerr << "This test needs an input data file.  Legend has it that one is found at " << url_test << endl;
	return 1;
    }

    std::string url = argv[1];

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
    const double from_gain_mVfC=4.7, to_gain_mVfC=14.0,
	from_shaping=1.0*units::microsecond, to_shaping=2.0*units::microsecond;
    
    
    // these are before hardware fix 
    for (int ind=2016; ind<= 2095; ++ind) { miscfgchan.push_back(ind); }
    for (int ind=2192; ind<= 2303; ++ind) { miscfgchan.push_back(ind); }
    for (int ind=2352; ind< 2400; ++ind) { miscfgchan.push_back(ind); }
    
    // these are after the hardware fix 
    // for (int ind=768; ind<864; ++ind) { miscfgchan.push_back(ind); }
    // for (int ind=2016; ind<2160; ++ind) { miscfgchan.push_back(ind); }
    // for (int ind=2176; ind<2400; ++ind) { miscfgchan.push_back(ind); }
    // for (int ind=2400+1440; ind<2400+1536; ++ind) { miscfgchan.push_back(ind); }
    // for (int ind=4800+1536; ind<4800+1728; ++ind) { miscfgchan.push_back(ind); }
    

    // hard-coded bad channels
    vector<int> bad_channels;
    for (int i=0;i!=wchans.size();i++){
      if (i>=7136 - 4800 && i <=7263 - 4800){
	if (i != 7200- 4800 && i!=7215 - 4800)
	  bad_channels.push_back(i+4800);
      }
    }
    
    // the following are additional channels after the hardware fix
    // bad_channels.push_back(300);
    // bad_channels.push_back(308);
    // bad_channels.push_back(310);
    // bad_channels.push_back(438);
    // bad_channels.push_back(336);
    // bad_channels.push_back(337);



    // Q&D RC+RC time constant - all have same.
    const double rcrc = 1.0*units::millisecond;
    vector<int> rcrcchans(nchans);
    std::iota(rcrcchans.begin(), rcrcchans.end(), 0);
    
    //harmonic noises
    vector<int> harmonicchans(uchans.size() + vchans.size());
    std::iota(harmonicchans.begin(), harmonicchans.end(), 0);
    
    vector<int> special_chans;
    special_chans.push_back(2240);

    SigProc::SimpleChannelNoiseDB::mask_t h36kHz(0,169,173);
    SigProc::SimpleChannelNoiseDB::mask_t h108kHz(0,513,516);
    SigProc::SimpleChannelNoiseDB::mask_t hspkHz(0,17,19);
    SigProc::SimpleChannelNoiseDB::multimask_t hharmonic;
    hharmonic.push_back(h36kHz);
    hharmonic.push_back(h108kHz);
    SigProc::SimpleChannelNoiseDB::multimask_t hspecial;
    hspecial.push_back(h36kHz);
    hspecial.push_back(h108kHz);
    hspecial.push_back(hspkHz);

    float u_resp_array[120]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.364382, 0.387949, 0.411053, 0.433979, 0.456863, 0.479746, 0.502641, 0.52554, 0.548441, 0.57134, 0.591765, 0.609448, 0.626848, 0.644094, 0.661364, 0.678859, 0.695231, 0.710462, 0.726147, 0.742373, 0.761332, 0.783313, 0.806325, 0.830412, 0.857676, 0.888412, 0.920705, 0.954624, 0.990242, 1.02766, 1.06121, 1.09027, 1.12037, 1.15157, 1.18392, 1.21748, 1.25229, 1.28824, 1.32509, 1.36256, 1.40051, 1.43907, 1.47857, 1.51933, 1.56134, 1.60404, 1.72665, 1.94005, 2.16994, 2.42041, 2.69475, 3.07222, 3.67375, 4.60766, 5.91864, 7.30178, 8.3715, 8.94736, 8.93705, 8.40339, 7.2212, 5.76382, 3.8931, 1.07893, -3.52481, -11.4593, -20.4011, -29.1259, -34.9544, -36.9358, -35.3303, -31.2068, -25.8614, -20.3613, -15.3794, -11.2266, -7.96091, -5.50138, -3.71143, -2.44637, -1.57662, -0.99733, -0.62554, -0.393562, -0.249715, -0.15914, -0.100771, -0.062443, -0.037283, -0.0211508, -0.0112448, -0.00552085, -0.00245133, -0.000957821, -0.000316912, -8.51679e-05, -2.21299e-05, -1.37496e-05, -1.49806e-05, -1.36935e-05, -9.66758e-06, -5.20773e-06, -7.4787e-07, 3.71199e-06, 8.17184e-06, 1.26317e-05, 1.70916e-05, 2.15514e-05, 2.60113e-05, 3.04711e-05};
    float v_resp_array[120]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0865303, 0.0925559, 0.0983619, 0.104068, 0.109739, 0.115403, 0.121068, 0.126735, 0.132403, 0.138072, 0.143739, 0.149408, 0.155085, 0.160791, 0.166565, 0.172454, 0.178514, 0.184795, 0.191341, 0.198192, 0.205382, 0.212944, 0.220905, 0.229292, 0.238129, 0.247441, 0.257256, 0.267601, 0.278502, 0.28999, 0.298745, 0.304378, 0.310105, 0.315921, 0.321818, 0.327796, 0.333852, 0.339967, 0.346098, 0.352169, 0.358103, 0.363859, 0.36945, 0.374915, 0.380261, 0.385401, 0.39016, 0.394378, 0.39804, 0.401394, 0.405145, 0.410714, 0.4205, 0.437951, 0.467841, 0.516042, 0.587738, 0.694157, 0.840763, 1.01966, 1.22894, 1.5612, 2.12348, 3.31455, 5.59355, 9.10709, 14.1756, 18.4603, 19.9517, 17.4166, 10.6683, 1.40656, -10.0638, -19.034, -23.654, -24.0558, -21.4418, -17.3229, -12.9485, -9.08912, -6.05941, -3.86946, -2.38669, -1.43678, -0.853335, -0.503951, -0.296551, -0.173029, -0.0990099, -0.0547172, -0.0287882, -0.0142758, -0.00661815, -0.00284757, -0.00115702, -0.000456456, -0.000183439, -8.04214e-05, -4.20533e-05, -2.62903e-05, -1.64098e-05, -6.68039e-06, 3.04903e-06, 1.27784e-05, 2.25079e-05, 3.22373e-05, 4.19667e-05, 5.16961e-05, 6.14255e-05, 7.11549e-05};
    WireCell::Waveform::realseq_t u_resp(nsamples);
    WireCell::Waveform::realseq_t v_resp(nsamples);
    for (int i=0;i!=120;i++){
      u_resp.at(i) = u_resp_array[i];
      v_resp.at(i) = v_resp_array[i];
    }
    WireCell::Waveform::compseq_t u_resp_freq = WireCell::Waveform::dft(u_resp);
    WireCell::Waveform::compseq_t v_resp_freq = WireCell::Waveform::dft(v_resp);
    

    int uplane_time_shift = 79;
    int vplane_time_shift = 82;

    
    // do the coherent subtraction
    
    std::vector< std::vector<int> > channel_groups;
    for (int i=0;i!=172;i++){
    //for (int i=150;i!=151;i++){
      std::vector<int> channel_group;
      for (int j=0;j!=48;j++){
	channel_group.push_back(i*48+j);
      }
      channel_groups.push_back(channel_group);
    }

    // Load up components.  Note, in a real app this is done as part
    // of factory + configurable and driven by user configuration.

    auto noise = new SigProc::SimpleChannelNoiseDB;
    // initialize
    noise->set_sampling(tick, nsamples);
    // set nominal baseline
    noise->set_nominal_baseline(uchans, unombl);
    noise->set_nominal_baseline(vchans, vnombl);
    noise->set_nominal_baseline(wchans, wnombl);

    noise->set_response(uchans,u_resp_freq);
    noise->set_response(vchans,v_resp_freq);

    noise->set_response_offset(uchans,uplane_time_shift);
    noise->set_response_offset(vchans,vplane_time_shift);

    noise->set_pad_window_front(uchans,20);
    noise->set_pad_window_back(uchans,10);
    noise->set_pad_window_front(vchans,10);
    noise->set_pad_window_back(vchans,10);
    noise->set_pad_window_front(wchans,10);
    noise->set_pad_window_back(wchans,10);
	
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

    // these are before Hardware Fix
    for (int i=0;i!=uchans.size();i++){
      if (uchans.at(i)<100){
    	noise->set_min_rms_cut_one(uchans.at(i),1);
    	noise->set_max_rms_cut_one(uchans.at(i),5);
      }else if (uchans.at(i)<2000){
    	noise->set_min_rms_cut_one(uchans.at(i),1.9);
    	noise->set_max_rms_cut_one(uchans.at(i),11);
      }else{
    	noise->set_min_rms_cut_one(uchans.at(i),0.9);
    	noise->set_max_rms_cut_one(uchans.at(i),5);
      }
    }
    for (int i=0;i!=vchans.size();i++){
      if (vchans.at(i)<290+uchans.size()){
    	noise->set_min_rms_cut_one(vchans.at(i),1);
    	noise->set_max_rms_cut_one(vchans.at(i),5);
      }else if (vchans.at(i)<2200+uchans.size()){
    	noise->set_min_rms_cut_one(vchans.at(i),1.9);
    	noise->set_max_rms_cut_one(vchans.at(i),11);
      }else{
    	noise->set_min_rms_cut_one(vchans.at(i),1);
    	noise->set_max_rms_cut_one(vchans.at(i),5);
      }
    }
    
    // these are the one after the Hardware Fix
    // for (int i=0;i!=uchans.size();i++){
    //   if (uchans.at(i)<600){
    // 	noise->set_min_rms_cut_one(uchans.at(i),1+(1.7-1)/600.*i);
    // 	noise->set_max_rms_cut_one(uchans.at(i),5);
    //   }else if (uchans.at(i)<1800){
    // 	noise->set_min_rms_cut_one(uchans.at(i),1.7);
    // 	noise->set_max_rms_cut_one(uchans.at(i),11);
    //   }else{
    // 	noise->set_min_rms_cut_one(uchans.at(i),1+ (1.7-1)/600.*(2399-i));
    // 	noise->set_max_rms_cut_one(uchans.at(i),5);
    //   }
    // }
    // for (int i=0;i!=vchans.size();i++){
    //   if (vchans.at(i)<600+uchans.size()){
    // 	noise->set_min_rms_cut_one(vchans.at(i),0.8+(1.7-0.8)/600.*i);
    // 	noise->set_max_rms_cut_one(vchans.at(i),5);
    //   }else if (vchans.at(i)<1800+uchans.size()){
    // 	noise->set_min_rms_cut_one(vchans.at(i),1.7);
    // 	noise->set_max_rms_cut_one(vchans.at(i),11);
    //   }else{
    // 	noise->set_min_rms_cut_one(vchans.at(i),0.8+ (1.7-0.8)/600.*(2399-i));
    // 	noise->set_max_rms_cut_one(vchans.at(i),5);
    //   }
    // }

    // Cut in W plane is the same before and after the Hardware Fix 
    noise->set_min_rms_cut(wchans,1.3);
    noise->set_max_rms_cut(wchans,8.0);
    
    
    shared_ptr<WireCell::IChannelNoiseDatabase> noise_sp(noise);

    auto one = new SigProc::Microboone::OneChannelNoise;
    one->set_channel_noisedb(noise_sp);
    shared_ptr<WireCell::IChannelFilter> one_sp(one);
    

    auto many = new SigProc::Microboone::CoherentNoiseSub;
    many->set_channel_noisedb(noise_sp);
    shared_ptr<WireCell::IChannelFilter> many_sp(many);


    SigProc::OmnibusNoiseFilter bus;
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
