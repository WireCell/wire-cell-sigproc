/** FIXME: this file is full of magic numbers and likely totally not
 * usable for detectors other than MicroBooNE. 
 * Modified from Microboone.cxx
 */

#include "WireCellSigProc/Microboone.h" 
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
WIRECELL_FACTORY(pdFembClockReSmp,
                 WireCell::SigProc::Protodune::FembClockReSmp,
                 WireCell::IChannelFilter, WireCell::IConfigurable)

using namespace WireCell::SigProc;

bool atLocalMinimum(WireCell::Waveform::realseq_t& signal,
                    int ind, int sideband){
    int nsiglen = signal.size();
    int left_ind = ind - sideband;
    int right_ind = ind + sideband;
    if(left_ind<0) left_ind = 0;
    if(right_ind>=nsiglen) right_ind = nsiglen-1;
    for(int i=left_ind; i<=right_ind; i++){
        if(signal.at(i) < signal.at(ind)) return false;
    }
    return true;
}

bool atLocalMaximum(WireCell::Waveform::realseq_t& signal,
                    int ind, int sideband){
    int nsiglen = signal.size();
    int left_ind = ind - sideband;
    int right_ind = ind + sideband;
    if(left_ind<0) left_ind =0;
    if(right_ind>=nsiglen) right_ind = nsiglen-1;
    for(int i=left_ind; i<=right_ind; i++){
        if(signal.at(i) > signal.at(ind)) return false;
    }
    return true;
}

int longestSequence(std::vector<float> arr, float fepsilon=1e-3)
{
  int n = arr.size();
  if(n==0) return 0;
  int longest = 0;
  int length = 1;
  for(int i = 1; i < n; i++){
      if(std::fabs(arr[i] - arr[i-1]) < fepsilon)
          length++;
      else
      {
          if(length > longest)
              longest = length;
          if(longest > n-1-i) // longer than the remaining
              return longest;
          length = 1;
      }
  }
  return (length > longest) ? length : longest;
}

bool Protodune::LinearInterpSticky(WireCell::Waveform::realseq_t& signal,
								   std::vector<std::pair<int,int> >& st_ranges, int ch){
	const int nsiglen = signal.size();
    // find ranges of sticky codes
    for(int i=0; i<nsiglen; i++){
      int val = signal.at(i);
      int mod = val % 64;
      if(mod==0 || mod==1 || mod==63
        || (ch==4 && mod==6) || (ch==159 && mod==6) || (ch==168 && mod==7) || (ch==323 && mod==24) 
        || (ch==164 && mod==36) || (ch==451 && mod==25) ){
        if (st_ranges.size()==0){
          st_ranges.push_back(std::make_pair(i,i));
        }else if ( (st_ranges.back().second + 1) == i){
          st_ranges.back().second = i;
        }else{
          st_ranges.push_back(std::make_pair(i,i));
        }
      }
    }

    // protect the ticks close to local minima/maxima
    // do not apply linear interpolation on them
    // std::vector<int> min_protect;
    // std::vector<int> max_protect;
    // int protect_length = 5;
    // for(auto const& rng: st_ranges){
    //     int pstart = rng.first  - protect_length; // backward protection
    //     int pend   = rng.second + protect_length;
    //     if(pstart<0) pstart=0;
    //     if(pend >= nsiglen) pend = nsiglen-1;
    //     for(int pind=pstart; pind<=pend; pind++){
    //         if(atLocalMinimum(signal, pind, protect_length)){
    //             for(int n=pind-protect_length; n<=pind+protect_length; n++){
    //                 if(min_protect.empty())
    //                     min_protect.push_back(n);
    //                 else if(n>min_protect.back()){
    //                     min_protect.push_back(n);
    //                 }
    //             }
    //         }
    //         else if(atLocalMaximum(signal, pind, protect_length)){
    //             for(int n=pind-protect_length; n<=pind+protect_length; n++){
    //                 if(max_protect.empty())
    //                     max_protect.push_back(n);
    //                 else if(n>max_protect.back()){
    //                     max_protect.push_back(n);
    //                 }
    //             }
    //         }
    //     }
    // }

    // if(ch==19){
    //     std::cerr << " ch " << ch << " , local min protection: " << min_protect.size() << std::endl;
    //     for(auto& p: min_protect){
    //         if(p>980 && p<1000) std::cerr << p << " ";
    //     }
    //     std::cerr << std::endl;
    // }

    std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);
    // linear-interpolation correction
    // int maxlen_equal_sticky = 5;
    for(auto const& rng: st_ranges){
      int start = rng.first -1;
      int end = rng.second + 1;
      // FIXME: do we need linear interp if rng.first==0?
      if (start >=0 && end <nsiglen){
        std::vector<float> digits;
        for(int i=start+1; i<end; i++){
            digits.push_back(signal.at(i));
        }
        // int length_equal_sticky = longestSequence(digits);
        auto min = std::max_element(digits.begin(), digits.end());
        auto max = std::min_element(digits.begin(), digits.end());
        double max_value = *max;
        double min_value = *min;

        double start_content = signal.at(start);
        double end_content = signal.at(end);
        for (int i = start+1; i<end; i++){
            bool fProtect = false; // if false, do linear interp
            // auto it1= std::find(max_protect.begin(), max_protect.end(), i);
            // auto it2= std::find(min_protect.begin(), min_protect.end(), i);
            // if(it1 != max_protect.end() || it2 != min_protect.end()){ // found protection in either min/max
            //     fProtect = true;
            // }

            // int length_sticky = rng.second - rng.first +1;
            // if(length_sticky>5) fProtect = false;
            // else if(std::fabs(signal.at(i) - temp.first)>40)
            //     fProtect = true;
            // else if(length_equal_sticky>=2)
            //     fProtect = false;
            // else
            //     fProtect = true;
            

            // bool fSkip = true;
            // if(signal.at(i) > temp.first){ // above baseline
            //     auto it= std::find(max_protect.begin(), max_protect.end(), i);
            //     if(it == max_protect.end()) fSkip = false; // no protection
            // }
            // else{
            //     auto it = std::find(min_protect.begin(), min_protect.end(), i);
            //     if(it == min_protect.end()) fSkip = false;
            // }

            if(rng.second == rng.first) fProtect = true; //single sticky, do FFT
            else{

                if(max_value - temp.first > 15){ // peak > 15, nearby > 2 rms
                    if( (start_content - temp.first > 2*temp.second)
                        && (end_content - temp.first > 2*temp.second) ){
                        fProtect = true;
                    }
                }
                else if(temp.first - min_value > 15){
                    if( (temp.first - end_content > 2*temp.second)
                        && (temp.first - end_content > 2*temp.second) ){
                        fProtect = true;
                    }

                }

            }


            if(!fProtect){
                double content = start_content + (end_content - start_content) /(end-start) *(i-start);
                signal.at(i) = content;
            }

        }
      }
      else if(start<0 && end <nsiglen){// sticky at the first tick
        for(int i=0; i<end; i++){
            signal.at(i) = signal.at(end);
        }

      }
      else if(start>=0 && end==nsiglen){// sticky at the last tick
        for(int i=start+1; i<end; i++){
            signal.at(i) = signal.at(start);
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

bool Protodune::FftScaling(WireCell::Waveform::realseq_t& signal, int nsamples){
	const int nsiglen = signal.size();
	auto tran = WireCell::Waveform::dft(signal);
	tran.resize(nsamples);
    if(nsiglen%2==0){ // ref test_zero_padding.cxx
    	std::rotate(tran.begin()+nsiglen/2, tran.begin()+nsiglen, tran.end());
    }
    else{
    	std::rotate(tran.begin()+(nsiglen+1)/2, tran.begin()+nsiglen, tran.end());
    }
    // inverse FFT
    auto signal_fc = WireCell::Waveform::idft(tran);
    WireCell::Waveform::scale(signal_fc, nsamples / nsiglen);
    signal = signal_fc;

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
    LinearInterpSticky(signal_lc, st_ranges, ch);
    FftInterpSticky(signal_lc, st_ranges);
    // FftShiftSticky(signal_lc, 0.5, st_ranges); // shift by 0.5 tick

    signal = signal_lc;

    return ret;
}



WireCell::Waveform::ChannelMaskMap Protodune::StickyCodeMitig::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}


Protodune::FembClockReSmp::FembClockReSmp(const std::string& anode, const std::string& noisedb)
    : ConfigFilterBase(anode, noisedb)
    , m_check_partial() // fixme, here too.
{
}
Protodune::FembClockReSmp::~FembClockReSmp()
{
}

WireCell::Waveform::ChannelMaskMap Protodune::FembClockReSmp::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    if( (ch>=2128 && ch<=2175) // W plane
    ||  (ch>=1520 && ch<=1559) // V plane
    ||  (ch>=440  && ch<=479)  // U plane
    ){
    	signal.resize(5996);
    	FftScaling(signal, 6000);
    }

    // correct rc undershoot
    auto spectrum = WireCell::Waveform::dft(signal);
    bool is_partial = m_check_partial(spectrum); // Xin's "IS_RC()"

    if(!is_partial){
        auto const& spec = m_noisedb->rcrc(ch);
        WireCell::Waveform::shrink(spectrum, spec);
    }

    // remove the DC component 
    spectrum.front() = 0;
    signal = WireCell::Waveform::idft(spectrum);    

    //Now calculate the baseline ...
    std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);
    auto temp_signal = signal;
    for (size_t i=0;i!=temp_signal.size();i++){
    if (fabs(temp_signal.at(i)-temp.first)>6*temp.second){
        temp_signal.at(i) = temp.first;
    }
    }
    float baseline = WireCell::Waveform::median_binned(temp_signal);
    //correct baseline
    WireCell::Waveform::increase(signal, baseline *(-1));


    // Now do the adaptive baseline for the bad RC channels
    if (is_partial) {
    // add something
    WireCell::Waveform::BinRange temp_chirped_bins;
    temp_chirped_bins.first = 0;
    temp_chirped_bins.second = signal.size();

    // auto wpid = m_anode->resolve(ch);
    // const int iplane = wpid.index();
    // if (iplane!=2) {        // not collection
    //     ret["lf_noisy"][ch].push_back(temp_chirped_bins);
    //     //std::cout << "Partial " << ch << std::endl;
    // }
    Microboone::SignalFilter(signal);
    Microboone::RawAdapativeBaselineAlg(signal);
    Microboone::RemoveFilterFlags(signal);
    }

    return ret;
}



WireCell::Waveform::ChannelMaskMap Protodune::FembClockReSmp::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
