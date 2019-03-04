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
WIRECELL_FACTORY(pdOneChannelNoise,
                 WireCell::SigProc::Protodune::OneChannelNoise,
                 WireCell::IChannelFilter, WireCell::IConfigurable)

using namespace WireCell::SigProc;

// adapted from WCP
// FIXME: some hardcoded 6000 ticks
bool LedgeIdentify(WireCell::Waveform::realseq_t& signal/*TH1F* h2*/, double baseline, int & LedgeStart, int & LedgeEnd){
    int ledge = false;
        int UNIT = 5;    // rebin unit
        int CONTIN = 20; // length of the continuous region
        int JUMPS = 4;   // how many bins can accidental jump
        std::vector<int> averaged; // store the rebinned waveform
        int up = signal.size()/UNIT;// h2->GetNbinsX()/UNIT;
    int nticks =  signal.size();// h2->GetNbinsX();
    // rebin 
    for(int i=0;i<up;i++){ 
                int temp = 0;
                for(int j=0;j<UNIT;j++){
                   temp += signal.at(i*UNIT+j); //h2->GetBinContent(i*UNIT+1+j);
                }
                averaged.push_back(temp);
        }
    // refine the selection cuts if there is a large signal
    auto imax = std::min_element(signal.begin(), signal.end());
    double max_value = *imax;

    // if(h2->GetMaximum()-baseline>1000) { CONTIN = 16; JUMPS = 5; }
    if(max_value-baseline>1000) { CONTIN = 16; JUMPS = 5; }
    // start judging
    int decreaseD = 0, tolerence=0;
        int start = 0;
        // int end = 0; // never used?
        for(int i=1;i<up-1;i++){
                if(averaged.at(i)<averaged.at(i-1)) {
                        if(decreaseD==0) start = i;
                        decreaseD +=1;
                }
         else {
                        if(averaged.at(i+1)<averaged.at(i-1)&&tolerence<JUMPS&&decreaseD>0){ // we can ignore several jumps in the decreasing region
                                decreaseD+=2;
                                tolerence++;
                                i = i+1;
                        }
                        else{
                                if(decreaseD>CONTIN){
                                        ledge = true;
                                        LedgeStart = start*UNIT;
                                        break;
                                }
                                else{
                                        decreaseD = 0;
                                        tolerence=0;
                                        start = 0;
                                        // end = 0;// end never used?
                                }
                        }
                }
        }
    // find the sharp start edge
     if(ledge &&LedgeStart>30){ 
                int edge = 0;
                int i = LedgeStart/UNIT-1;
                if(averaged.at(i)>averaged.at(i-1)&&averaged.at(i-1)>averaged.at(i-2)){ // find a edge
                        edge = 1;
                }
                if(edge == 0) ledge = false; // if no edge, this is not ledge
                if((averaged.at(i)-averaged.at(i-2)<10*UNIT)&&(averaged.at(i)-averaged.at(i-3)<10*UNIT)) // slope cut
                        ledge = false;
                if(averaged.at(LedgeStart/UNIT)-baseline*UNIT>200*UNIT) ledge = false; // ledge is close to the baseline
        }
    // test the decay time
    if(ledge &&LedgeStart>20){
                double height = 0;
                if(LedgeStart<5750) { // calculate the height of edge
                        // double tempHeight = h2 ->GetBinContent(LedgeStart+1+200) +  h2 ->GetBinContent(LedgeStart+1+220) +  h2 ->GetBinContent(LedgeStart+1+180) +  h2 ->GetBinContent(LedgeStart+1+240);
                        // height = h2 ->GetBinContent(LedgeStart+1) - tempHeight/4;
                        double tempHeight = signal.at(LedgeStart+200) +  signal.at(LedgeStart+220) +  signal.at(LedgeStart+180) +  signal.at(LedgeStart+240);
                        height = signal.at(LedgeStart) - tempHeight/4;                        
            height /= 0.7;
                }
                // else height =  h2 ->GetBinContent(LedgeStart+1) -  h2 ->GetBinContent(6000);
                else height =  signal.at(LedgeStart) -  signal.back();
                if(height<0) height = 80; // norminal value
                if(height>30&&LedgeStart<5900){ // test the decay with a relatively large height
                        double height50 = 0, height100 = 0;
                        // height50 =  h2 ->GetBinContent(LedgeStart+51);
                        // height100 =  h2 ->GetBinContent(LedgeStart+101);
                        // double height50Pre =   h2 ->GetBinContent(LedgeStart+1)- height*(1-exp(-50/100.)); // minimum 100 ticks decay time
                        // double height100Pre =   h2 ->GetBinContent(LedgeStart+1) - height*(1-exp(-100./100)); // minimum 100 ticks decay time

                        height50 =  signal.at(LedgeStart+50);
                        height100 =  signal.at(LedgeStart+100);
                        double height50Pre =   signal.at(LedgeStart)- height*(1-std::exp(-50/100.)); // minimum 100 ticks decay time
                        double height100Pre =   signal.at(LedgeStart) - height*(1-std::exp(-100./100)); // minimum 100 ticks decay time                        
            // if the measured is much smaller than the predicted, this is not ledge
                        if(height50-height50Pre<-8) ledge = false; 
                        if(height100-height100Pre<-8)  ledge = false;
                }
        }

    // determine the end of ledge
    // case 1: find a jump of 10 ADC in the rebinned waveform
    // case 2: a continuous 20 ticks has an average close to baseline, and a RMS larger than 3 ADC
    // case 3: reaching the tick 6000 but neither 1 nor 2 occurs
    if(ledge){
        LedgeEnd = 0;
        for(int i = LedgeStart/UNIT; i<up-1; i++){ // case 1
            if(averaged.at(i+1)-averaged.at(i)>50) { 
                LedgeEnd = i*UNIT+5;
                break;
            }
        }
        if(LedgeEnd == 0) { // not find a jump, case 2
            // double tempA[20];
            WireCell::Waveform::realseq_t tempA(20);
            for(int i = LedgeStart+80;i<nticks-20;i+=20){
                for(int j=i;j<20+i;j++){
                    // tempA[j-i] = h2->GetBinContent(j+1);
                    tempA.at(j-i) = signal.at(j);
                }
                auto wfinfo = WireCell::Waveform::mean_rms(tempA);
                // if(TMath::Mean(20,tempA)-baseline<2&&TMath::RMS(20,tempA)>3){
                if(wfinfo.first - baseline < 2 && wfinfo.second > 3){
                    LedgeEnd = i;
                    break;
                }
            }
        }
        if(LedgeEnd == 0) LedgeEnd = 6000;
    }
    // done, release the memory
    // vector<int>(averaged).swap(averaged); // is it necessary?
    return ledge;

}



// adapted from WCP
// int judgePlateau(int channel, TH1F* h2,double baseline, double & PlateauStart, double & PlateauStartEnd){
//         int continueN = 0;
//         int threshold = 200;
//         int maximumF  = 50;
//         int maxBin = h2->GetMaximumBin();
//         for(int i=maxBin+10;i<5880&&i<maxBin+500;i++){
//                 int plateau = 1;
//                 int max = 0, min = 10000;
//                 for(int j=i;j<i+20;j++){
//                         int binC = h2->GetBinContent(j+1);
//                         if(binC<baseline+threshold||binC>h2->GetMaximum()-500) {
//                                 plateau = 0;
//                                 break;
//                         }
//                         if(binC>max) max = binC;
//                         if(binC<min) min = binC;
//                 }
//                 if(plateau==1&&max-min<maximumF){ // plateau found
//                         PlateauStart = i;
//                         PlateauStartEnd = i+20;
//                         for(int k = i+20; k<6000;k++){
//                                 if( h2->GetBinContent(k+1)<baseline+threshold){
//                                         PlateauStartEnd = k-1;
//                                         break;
//                                 }
//                         }
//                         return 1;
//                 }
//         }
//         return 0;
// }

// bool atLocalMinimum(WireCell::Waveform::realseq_t& signal,
//                     int ind, int sideband){
//     int nsiglen = signal.size();
//     int left_ind = ind - sideband;
//     int right_ind = ind + sideband;
//     if(left_ind<0) left_ind = 0;
//     if(right_ind>=nsiglen) right_ind = nsiglen-1;
//     for(int i=left_ind; i<=right_ind; i++){
//         if(signal.at(i) < signal.at(ind)) return false;
//     }
//     return true;
// }

// bool atLocalMaximum(WireCell::Waveform::realseq_t& signal,
//                     int ind, int sideband){
//     int nsiglen = signal.size();
//     int left_ind = ind - sideband;
//     int right_ind = ind + sideband;
//     if(left_ind<0) left_ind =0;
//     if(right_ind>=nsiglen) right_ind = nsiglen-1;
//     for(int i=left_ind; i<=right_ind; i++){
//         if(signal.at(i) > signal.at(ind)) return false;
//     }
//     return true;
// }

// int longestSequence(std::vector<float> arr, float fepsilon=1e-3)
// {
//   int n = arr.size();
//   if(n==0) return 0;
//   int longest = 0;
//   int length = 1;
//   for(int i = 1; i < n; i++){
//       if(std::fabs(arr[i] - arr[i-1]) < fepsilon)
//           length++;
//       else
//       {
//           if(length > longest)
//               longest = length;
//           if(longest > n-1-i) // longer than the remaining
//               return longest;
//           length = 1;
//       }
//   }
//   return (length > longest) ? length : longest;
// }

bool Protodune::LinearInterpSticky(WireCell::Waveform::realseq_t& signal,
								   WireCell::Waveform::BinRangeList& rng_list, int ch){
	const int nsiglen = signal.size();
    // // find ranges of sticky codes
    // for(int i=0; i<nsiglen; i++){
    //   int val = signal.at(i);
    //   int mod = val % 64;
    //   if(mod==0 || mod==1 || mod==63
    //     || (ch==4 && mod==6) || (ch==159 && mod==6) || (ch==168 && mod==7) || (ch==323 && mod==24) 
    //     || (ch==164 && mod==36) || (ch==451 && mod==25) ){
    //     if (st_ranges.size()==0){
    //       st_ranges.push_back(std::make_pair(i,i));
    //     }else if ( (st_ranges.back().second + 1) == i){
    //       st_ranges.back().second = i;
    //     }else{
    //       st_ranges.push_back(std::make_pair(i,i));
    //     }
    //   }
    // }

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
    for(auto const& rng: rng_list){
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
	                 WireCell::Waveform::BinRangeList& rng_list){
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
    for (size_t j = 0; j<rng_list.size();j++){
      int start = rng_list.at(j).first -1 ;
      int end = rng_list.at(j).second + 1 ;
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
    : m_anode_tn(anode)
    , m_noisedb_tn(noisedb)
//    , m_check_chirp() // fixme, there are magic numbers hidden here
//    , m_check_partial() // fixme, here too.
{
}
Protodune::StickyCodeMitig::~StickyCodeMitig()
{
}

void Protodune::StickyCodeMitig::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
    }

    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    m_extra_stky.clear();
    auto jext = cfg["extra_stky"];
    if(!jext.isNull()){
        for(auto jone: jext) {
            int channel = jone["channel"].asInt();
            // std::cerr << "[wgu] ch# " << channel << " has " << jone["bits"].size() << " extra stky bits:" << std::endl;
            for(auto bit: jone["bits"]){
                // std::cerr << "[wgu] " << bit.asInt() << std::endl;
                m_extra_stky[channel].push_back((short int)bit.asInt());
            }
        }
    }

}
WireCell::Configuration Protodune::StickyCodeMitig::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;    
    return cfg;
}

WireCell::Waveform::ChannelMaskMap Protodune::StickyCodeMitig::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    const int nsiglen = signal.size();
    // tag sticky bins
    WireCell::Waveform::BinRange sticky_rng;
    WireCell::Waveform::BinRangeList sticky_rng_list;
    std::vector<short int> extra;
    if (m_extra_stky.find(ch) != m_extra_stky.end()){
        extra = m_extra_stky.at(ch);
    }

    for(int i=0; i<nsiglen; i++){
      int val = signal.at(i);
      int mod = val % 64;
      auto it = std::find(extra.begin(), extra.end(), mod);
      bool is_stky = (mod==0 || mod==1 || mod==63 || it!=extra.end());
      if(is_stky){
        if (sticky_rng_list.empty()){
          sticky_rng_list.push_back(std::make_pair(i,i));
        }else if ( (sticky_rng_list.back().second + 1) == i){
          sticky_rng_list.back().second = i;
        }else{
          sticky_rng_list.push_back(std::make_pair(i,i));
        }
      }
    }


    // auto signal_lc = signal; // copy, need to keep original signal
    LinearInterpSticky(signal, sticky_rng_list, ch);
    FftInterpSticky(signal, sticky_rng_list);
    // FftShiftSticky(signal_lc, 0.5, st_ranges); // alternative approach, shift by 0.5 tick
    // signal = signal_lc;

    for(auto rng: sticky_rng_list){
        if(rng.second-rng.first>5){
            ret["sticky"][ch].push_back(rng);
        }
    }

    return ret;
}



WireCell::Waveform::ChannelMaskMap Protodune::StickyCodeMitig::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}


Protodune::OneChannelNoise::OneChannelNoise(const std::string& anode, const std::string& noisedb)
    : ConfigFilterBase(anode, noisedb)
    , m_check_partial() // fixme, here too.
{
}
Protodune::OneChannelNoise::~OneChannelNoise()
{
}

WireCell::Waveform::ChannelMaskMap Protodune::OneChannelNoise::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    // do we need a nominal baseline correction?
    // float baseline = m_noisedb->nominal_baseline(ch);

    // correct the FEMB 302 clock issue
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
        auto const& spec = m_noisedb->rcrc(ch); // rc_layers set to 1 in nf.jsonnet
        WireCell::Waveform::shrink(spectrum, spec);
    }

    // remove the "50kHz" noise in some collection channels
    // FIXME: do we need a channel list input?
    auto wpid = m_anode->resolve(ch);      
    const int iplane = wpid.index();
    if(iplane==2){
    auto mag = WireCell::Waveform::magnitude(spectrum);
    Microboone::RawAdapativeBaselineAlg(mag); // subtract "linear" background in spectrum

    for(int i=0; i<57; i++){ // 150 - 3000th freq bin
        int nslice = 50;
        int istart = 150 + nslice*i;
        int iend = istart + nslice;
        // std::cerr << istart << " " << iend << std::endl;
        WireCell::Waveform::realseq_t mag_slice(nslice); // slice of magnitude spectrum
        std::copy(mag.begin() + istart, mag.begin() + iend, mag_slice.begin());
        std::pair<double,double> stat = WireCell::Waveform::mean_rms(mag_slice);
        // std::cerr << "[wgu] mean/rms: " << stat.first << " " << stat.second << std::endl;
        for(int j=0; j<nslice; j++){
            float content = mag_slice.at(j) - stat.first;
             
            if(iend<1000){
                if(content>2000 && content>5.*stat.second){
                int tbin = istart + j;
                spectrum.at(tbin).real(0);
                spectrum.at(tbin).imag(0);
                spectrum.at(6000+1-tbin).real(0); // FIXME: assuming 6000 ticks
                spectrum.at(6000+1-tbin).imag(0);
                // std::cerr << "[wgu] chan: " << ch << " , freq tick: " << tbin << " , amp: " << content << std::endl;
                }
            }
            else if(content>250 && content>10.*stat.second){
                spectrum.at(j).real(0);
                spectrum.at(j).imag(0);
                spectrum.at(6000+1-j).real(0); // FIXME: assuming 6000 ticks
                spectrum.at(6000+1-j).imag(0); 
            }
        }
    }

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

    if (iplane!=2) {        // not collection
        ret["lf_noisy"][ch].push_back(temp_chirped_bins);
        //std::cout << "Partial " << ch << std::endl;
    }
    Microboone::SignalFilter(signal);
    Microboone::RawAdapativeBaselineAlg(signal);
    Microboone::RemoveFilterFlags(signal);
    }


    // ledge identification
    WireCell::Waveform::BinRange ledge_bins;
    bool is_ledge = LedgeIdentify(signal, 0., ledge_bins.first, ledge_bins.second);
    if(is_ledge){
        // FIXME: do we need collection plane only?
        ret["ledge"][ch].push_back(ledge_bins);
        // std::cerr << "[wgu] ledge found in ch "<< ch << " , bins= [" << ledge_bins.first << " , " << ledge_bins.second << " ]"<< std::endl;
    }

    return ret;
}



WireCell::Waveform::ChannelMaskMap Protodune::OneChannelNoise::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
