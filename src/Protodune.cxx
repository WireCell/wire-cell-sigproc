/** FIXME: this file is full of magic numbers and likely totally not
 * usable for detectors other than MicroBooNE. 
 * Modified from Microboone.cxx
 */
 
#include "WireCellSigProc/Protodune.h"
#include "WireCellSigProc/Derivations.h"

#include "WireCellUtil/NamedFactory.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <set>

WIRECELL_FACTORY(pdStickyCodeMitig,
                 WireCell::SigProc::Protodune::StickyCodeMitig,
                 WireCell::IChannelFilter, WireCell::IConfigurable)

using namespace WireCell::SigProc;

bool Protodune::LinearInterpSticky(WireCell::Waveform::realseq_t& signal,
								   std::vector<std::pair<int,int> >& st_ranges){
	const int nsiglen = signal.size();
    // identify ranges of sticky codes (st_ranges)
    // std::vector<std::pair<int,int> > st_ranges;
    for(int i=0; i<nsiglen; i++){
      int val = signal.at(i);
      int mod = val % 64;
      if(mod==0 || mod==1 || mod==63){
        if (st_ranges.size()==0){
          st_ranges.push_back(std::make_pair(i,i));
        }else if ( (st_ranges.back().second + 1) == i){
          st_ranges.back().second = i;
        }else{
          st_ranges.push_back(std::make_pair(i,i));
        }
      }
    }

    // linear-interpolation correction first
    for (size_t j = 0; j<st_ranges.size();j++){
      int start = st_ranges.at(j).first -1 ;
      int end = st_ranges.at(j).second + 1 ;
      if (start >=0 && end <nsiglen){
        double start_content = signal.at(start);
        double end_content = signal.at(end);
        for (int i = start+1; i<end; i++){
    	    double content = start_content + (end_content - start_content) /(end-start) *(i-start);
    	    signal.at(i) = content;
        }
      }
    }
    return true;
}

bool Protodune::FftInterpSticky(WireCell::Waveform::realseq_t& signal,
	                 std::vector<std::pair<int,int> >& st_ranges){
	const int nsiglen = signal.size();
    // group into two subsamples ("even" & "odd")
    int nsublen = nsiglen/2;
    int nsublen2 = nsiglen - nsublen;
    WireCell::Waveform::realseq_t signal_even(nsublen); // 0-th bin, 2, 4, etc.
    WireCell::Waveform::realseq_t signal_odd(nsublen2);
    for(int j=0; j<nsiglen; j++){
    	if(j%2==0) signal_even.at(j/2) = signal.at(j);
    	else signal_odd.at((j-1)/2) = signal.at(j);
    }

    // dft resampling for "even", see example in test_zero_padding.cxx
    auto tran_even = WireCell::Waveform::dft(signal_even);
    tran_even.resize(nsublen*2);
    if(nsublen%2==0){
    	std::rotate(tran_even.begin()+nsublen/2, tran_even.begin()+nsublen, tran_even.end());
    }
    else{
    	std::rotate(tran_even.begin()+(nsublen+1)/2, tran_even.begin()+nsublen, tran_even.end());
    }
    // inverse FFT
    auto signal_even_fc = WireCell::Waveform::idft(tran_even);
    float scale = tran_even.size() / nsublen;
    WireCell::Waveform::scale(signal_even_fc, scale);

    // similar for "odd"
    auto tran_odd = WireCell::Waveform::dft(signal_odd);
    tran_odd.resize(nsublen2*2);
    if(nsublen2%2==0){
    	std::rotate(tran_odd.begin()+nsublen2/2, tran_odd.begin()+nsublen2, tran_odd.end());
    }
    else{
    	std::rotate(tran_odd.begin()+(nsublen2+1)/2, tran_odd.begin()+nsublen2, tran_odd.end());
    }
    auto signal_odd_fc = WireCell::Waveform::idft(tran_odd);
    float scale2 = tran_odd.size() / nsublen2;
	WireCell::Waveform::scale(signal_odd_fc, scale2);

	// replace the linear interpolation with dft interpolation
    for (size_t j = 0; j<st_ranges.size();j++){
      int start = st_ranges.at(j).first -1 ;
      int end = st_ranges.at(j).second + 1 ;
      if (start >=0 && end <=nsiglen){
        for (int i = start+1; i<end; i++){
        	if(i%2==0){// predict "even" with "odd"
        		signal.at(i) = signal_odd_fc.at(i-1);
        	}
        	else{
        		signal.at(i) = signal_even_fc.at(i);
        	}
        }
      }
    }

	return true;
}


bool Protodune::FftShiftSticky(WireCell::Waveform::realseq_t& signal,
	                double toffset,
	                std::vector<std::pair<int,int> >& st_ranges){
	const int nsiglen = signal.size();
    // group into two subsamples ("even" & "odd")
    int nsublen = nsiglen/2;
    int nsublen2 = nsiglen - nsublen;
    WireCell::Waveform::realseq_t signal_even(nsublen); // 0-th bin, 2, 4, etc.
    WireCell::Waveform::realseq_t signal_odd(nsublen2);
    for(int j=0; j<nsiglen; j++){
    	if(j%2==0) signal_even.at(j/2) = signal.at(j);
    	else signal_odd.at((j-1)/2) = signal.at(j);
    }

    // dft shift for "even"
    auto tran_even = WireCell::Waveform::dft(signal_even);
    double f0 = 1./nsublen;
    const double PI = std::atan(1.0)*4;
    for(size_t i=0; i<tran_even.size(); i++){
    	double fi = i * f0;
    	double omega = 2*PI*fi;
    	auto z = tran_even.at(i);
    	std::complex<float> z1(0, omega*toffset);
    	// std::complex<double> z2 = z * std::exp(z1);
    	tran_even.at(i) = z * std::exp(z1);
    }
    // inverse FFT
    auto signal_even_fc = WireCell::Waveform::idft(tran_even);
    // float scale = 1./tran_even.size();
    // WireCell::Waveform::scale(signal_even_fc, 1./nsublen);    

    // similar to "odd"
    auto tran_odd = WireCell::Waveform::dft(signal_odd);
    f0 = 1./nsublen2;
    for(size_t i=0; i<tran_odd.size(); i++){
    	double fi = i * f0;
    	double omega = 2*PI*fi;
    	auto z = tran_odd.at(i);
    	std::complex<float> z1(0, omega*toffset);
    	// std::complex z2 = z * std::exp(z1);
    	tran_odd.at(i) = z * std::exp(z1);
    }
    //
    auto signal_odd_fc = WireCell::Waveform::idft(tran_odd);
    // float scale = 1./tran_odd.size();
    // WireCell::Waveform::scale(signal_odd_fc, 1./nsublen2);

	// replace the linear interpolation with dft interpolation
    for (size_t j = 0; j<st_ranges.size();j++){
      int start = st_ranges.at(j).first -1 ;
      int end = st_ranges.at(j).second + 1 ;
      if (start >=0 && end <=nsiglen){
        for (int i = start+1; i<end; i++){
        	if(i%2==0){ // predict "even" with "odd"
        		int ind = (i-1)/2. - toffset;
        		if(ind>=0 && ind<nsublen2) signal.at(i) = signal_odd_fc.at(ind);
        		// std::cerr << "lc:fc= " << signal_lc.at(i) << " " << signel_even_fc.at() << std::endl;
        	}
        	else{
        		int ind = i/2. - toffset;
        		if(ind>=0 && ind<nsublen) signal.at(i) = signal_even_fc.at(ind);
        	}
        }
      }
    }

	return true;
}



/* 
 * Classes
 */


/*
 * Configuration base class used for a couple filters
 */
Protodune::ConfigFilterBase::ConfigFilterBase(const std::string& anode, const std::string& noisedb)
    : m_anode_tn(anode)
    , m_noisedb_tn(noisedb) {}
Protodune::ConfigFilterBase::~ConfigFilterBase() {}
void Protodune::ConfigFilterBase::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);
    //std::cerr << "ConfigFilterBase: \n" << cfg << "\n";
}
WireCell::Configuration Protodune::ConfigFilterBase::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn; 
    cfg["noisedb"] = m_noisedb_tn; 
    return cfg;
}



Protodune::StickyCodeMitig::StickyCodeMitig(const std::string& anode, const std::string& noisedb)
    : ConfigFilterBase(anode, noisedb)
//    , m_check_chirp() // fixme, there are magic numbers hidden here
//    , m_check_partial() // fixme, here too.
{
}
Protodune::StickyCodeMitig::~StickyCodeMitig()
{
}

WireCell::Waveform::ChannelMaskMap Protodune::StickyCodeMitig::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;


    auto signal_lc = signal; // copy, need to keep original signal
    std::vector<std::pair<int,int> > st_ranges; // ranges of sticky codes
    LinearInterpSticky(signal_lc, st_ranges);
    FftInterpSticky(signal_lc, st_ranges);
    // FftShiftSticky(signal_lc, 0.5, st_ranges); // shift by 0.5 tick

    signal = signal_lc;

    return ret;
}



WireCell::Waveform::ChannelMaskMap Protodune::StickyCodeMitig::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}




// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
