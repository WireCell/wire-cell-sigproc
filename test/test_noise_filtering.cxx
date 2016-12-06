


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

int main(int argc, char* argv[])
{
  std::vector<std::vector<float>> horigs;
#include "example-chirping-48.h"
  assert(horigs.size()==48);
  assert(horigs.at(0).size()==9594);
#include "example-noisy-48.h"
  assert(horigs.size()==96);
#include "example-misconfig-48.h"
  assert(horigs.size()==144);
#include "example-rcrc-48.h"
  assert(horigs.size()==192);
  
  // std::cout << horigs.size() << " " << horigs.at(0).size() << std::endl;
  
  ITrace::vector traces;
  int chindex = 624;
  for (int ich = 0; ich!=48; ich++){
    ITrace::ChargeSequence charges;
    for (int itick =0; itick!=9594;itick++){
      auto q = horigs.at(ich+48).at(itick);
      charges.push_back(q);
    }
    SimpleTrace *st = new SimpleTrace(chindex+ich,0.0,charges);
    traces.push_back(ITrace::pointer(st));
  }
  
  chindex = 720;
  for (int ich = 0; ich!=48; ich++){
    ITrace::ChargeSequence charges;
    for (int itick =0; itick!=9594;itick++){
      auto q = horigs.at(ich).at(itick);
      charges.push_back(q);
    }
    SimpleTrace *st = new SimpleTrace(chindex+ich,0.0,charges);
    traces.push_back(ITrace::pointer(st));
  }
  chindex = 2016;
  for (int ich = 0; ich!=48; ich++){
    ITrace::ChargeSequence charges;
    for (int itick =0; itick!=9594;itick++){
      auto q = horigs.at(ich+96).at(itick);
      charges.push_back(q);
    }
    SimpleTrace *st = new SimpleTrace(chindex+ich,0.0,charges);
    traces.push_back(ITrace::pointer(st));
  }
  chindex = 7728;
  for (int ich = 0; ich!=48; ich++){
    ITrace::ChargeSequence charges;
    for (int itick =0; itick!=9594;itick++){
      auto q = horigs.at(ich+144).at(itick);
      charges.push_back(q);
    }
    SimpleTrace *st = new SimpleTrace(chindex+ich,0.0,charges);
    traces.push_back(ITrace::pointer(st));
  }


  SimpleFrame* sf = new SimpleFrame(0, 0, traces);


  

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
  
  auto one = new SigProc::Microboone::OneChannelNoise;
  one->set_channel_noisedb(noise_sp);
  shared_ptr<WireCell::IChannelFilter> one_sp(one);
  
  auto many = new SigProc::Microboone::CoherentNoiseSub;
  shared_ptr<WireCell::IChannelFilter> many_sp(many);
  
  
  SigProc::OmnibusNoiseFilter bus;
  bus.set_channel_filters({one_sp});
  bus.set_grouped_filters({many_sp});
  bus.set_channel_noisedb(noise_sp);
  
  IFrame::pointer frame = IFrame::pointer(sf);
  
  IFrame::pointer quiet;
  bus(frame, quiet);
  Assert(quiet);

  std::vector<std::vector<float>> hfilts;
#include "example-noisy-48-filtered.h"
  assert(hfilts.size()==48);
  assert(hfilts.at(0).size()==9594);
#include "example-chirping-48-filtered.h"
  assert(hfilts.size()==96);
#include "example-misconfig-48-filtered.h"
  assert(hfilts.size()==144);
#include "example-rcrc-48-filtered.h"
  assert(hfilts.size()==192);
  // test ...
  auto traces1 = quiet->traces();
  int counter = 0;
 
  for (auto trace : *traces1.get()) {
    int tbin = trace->tbin();
    int ch = trace->channel();
    auto charges = trace->charge();
    // std::cout << ch << " " << counter << " " << charges.size() << " " << hfilts.at(counter).size() << " " << charges.at(0) << " " << hfilts.at(counter).at(0) << std::endl;
    for (int i=0;i!=9594;i++){
      assert( fabs(charges.at(i) - hfilts.at(counter).at(i)) < 1.0 );
    }


    counter ++;
  }

  return 0;
}
