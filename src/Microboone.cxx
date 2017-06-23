/** FIXME: this file is full of magic numbers and likely totally not
 * usable for detectors other than MicroBooNE. 
 */

#include "WireCellSigProc/Microboone.h"
#include "WireCellSigProc/Derivations.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <set>

using namespace WireCell::SigProc;

double filter_time(double freq){
    double a = 0.143555;
    double b = 4.95096;
    return (freq>0)*exp(-0.5*pow(freq/a,b));
}

double filter_low(double freq){
    return  1-exp(-pow(freq/0.06,2));
}


bool Microboone::Subtract_WScaling(WireCell::IChannelFilter::channel_signals_t& chansig,
				   const WireCell::Waveform::realseq_t& medians)
{
    double ave_coef = 0;
    double_t ave_coef1 = 0;
    std::map<int,double> coef_all;

    const int nbin = medians.size();

    for (auto it: chansig) {
	int ch = it.first;
	WireCell::IChannelFilter::signal_t& signal = it.second;
	
	
	double sum2 = 0;
	double sum3 = 0;
	double coef = 0;

	std::pair<double,double> temp = Derivations::CalcRMS(signal);
	//std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);
	
	// if ( abs(ch-6117)<5)
	//     std::cout << ch << " " << temp.first << " " << temp.second << " "  << std::endl;
	
	for (int j=0;j!=nbin;j++){
	    if (fabs(signal.at(j)) < 4 * temp.second){
		sum2 += signal.at(j) * medians.at(j);
		sum3 += medians.at(j) * medians.at(j);
	    }
	}
	if (sum3 >0) {
	    coef = sum2/sum3;
	}
	// protect against the extreme cases
	// if (coef < 0.6) coef = 0.6;
	// if (coef > 1.5) coef = 1.5;

	coef_all[ch] = coef;
	if (coef != 0){
	    ave_coef += coef;
	    ave_coef1 ++;
	}
    }
    if (ave_coef1 > 0) {
	ave_coef = ave_coef / ave_coef1;
    }
  
    for (auto it: chansig) {
	int ch = it.first;
	WireCell::IChannelFilter::signal_t& signal = it.second;
	float scaling;
	if (ave_coef != 0 ){
	    scaling = coef_all[ch]/ave_coef;
	    // add some protections ... 
	    if (scaling < 0.5 ) scaling = 0.5;
	    if (scaling > 1.5) scaling = 1.5;
	}else{
	    scaling = 0;
	}
	// if ( abs(ch-6117)<5)
	//     std::cout << ch << " " << scaling << " "  << std::endl;
	//scaling = 1.0;
	for (int i=0;i!=nbin;i++){
	    if (fabs(signal.at(i)) > 0.001) {
		signal.at(i) = signal.at(i) - medians.at(i) * scaling;
	    }
	}
	// std::cout << signal.at(0) << std::endl;
	//std::cout << it.second.at(0) << std::endl;
	chansig[ch] = signal;
    }

    // for (auto it: chansig){
    //   std::cout << "Xin2 " << it.second.at(0) << std::endl;
    // }
  
    return true;
}

bool Microboone::SignalProtection(WireCell::Waveform::realseq_t& medians, const WireCell::Waveform::compseq_t& respec, int res_offset, int pad_f, int pad_b)
{
   
  
    // WireCell::Waveform::realseq_t temp1;
    // for (int i=0;i!=medians.size();i++){
    //   if (fabs(medians.at(i) - mean) < 4.5*rms)
    //     temp1.push_back(medians.at(i));
    // }
    // temp = WireCell::Waveform::mean_rms(temp1);
    // mean = temp.first;
    // rms = temp.second;
  
    //std::cout << temp.first << " " << temp.second << std::endl;
    int nbin = medians.size();

    int protection_factor = 5.0;
    float upper_decon_limit = 0.05;
    float upper_adc_limit = 15;
    float min_adc_limit = 50;


    std::vector<int> signals;
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin,0);

    
    

     // calculate the RMS 
    std::pair<double,double> temp = Derivations::CalcRMS(medians);
    double mean = temp.first;
    double rms = temp.second;
    
    float limit;
    if (protection_factor*rms > upper_adc_limit){
	limit = protection_factor*rms;
    }else{
	limit = upper_adc_limit;
    }
    if (min_adc_limit < limit){
	limit = min_adc_limit;
    }

    for (int j=0;j!=nbin;j++) {
	float content = medians.at(j);
	if (fabs(content-mean)>limit){
	    //protection_factor*rms) {
	    medians.at(j) = 0; 
	    signalsBool.at(j) = 1;
	    // add the front and back padding
	    for (int k=0;k!=pad_b;k++){
		int bin = j+k+1;
		if (bin > nbin-1) bin = nbin-1;
		signalsBool.at(bin) = 1;
	    }
	    for (int k=0;k!=pad_f;k++){
		int bin = j-k-1;
		if (bin <0) { bin = 0; }
		signalsBool.at(bin) = 1;
	    }
	}
    }
  
    // the deconvolution protection code ... 
    if (respec.size()>0){
	//std::cout << nbin << std::endl;

     	WireCell::Waveform::compseq_t medians_freq = WireCell::Waveform::dft(medians);
     	WireCell::Waveform::shrink(medians_freq,respec);
	
    	for (size_t i=0;i!=medians_freq.size();i++){
    	    double freq;
    	    // assuming 2 MHz digitization
    	    if (i <medians_freq.size()/2.){
    		freq = i/(1.*medians_freq.size())*2.;
    	    }else{
    		freq = (medians_freq.size() - i)/(1.*medians_freq.size())*2.;
    	    }
    	    std::complex<float> factor = filter_time(freq)*filter_low(freq);
    	    medians_freq.at(i) = medians_freq.at(i) * factor;
    	}
    	WireCell::Waveform::realseq_t medians_decon = WireCell::Waveform::idft(medians_freq);
	
    	temp = Derivations::CalcRMS(medians_decon);
    	mean = temp.first;
    	rms = temp.second;
	
    	if (protection_factor*rms > upper_decon_limit){
    	    limit = protection_factor*rms;
    	}else{
    	    limit = upper_decon_limit;
    	}
	
    	for (int j=0;j!=nbin;j++) {
    	    float content = medians_decon.at(j);
    	    if ((content-mean)>limit){
    		int time_bin = j + res_offset;
		if (time_bin >= nbin) time_bin -= nbin;
    		medians.at(time_bin) = 0; 
    		signalsBool.at(time_bin) = 1;
    		// add the front and back padding
    		for (int k=0;k!=pad_b;k++){
    		    int bin = time_bin+k+1;
    		    if (bin > nbin-1) bin = nbin-1;
    		    signalsBool.at(bin) = 1;
    		}
    		for (int k=0;k!=pad_f;k++){
    		    int bin = time_bin-k-1;
    		    if (bin <0) { bin = 0; }
    		    signalsBool.at(bin) = 1;
    		}
    	    }
    	}
    }

    // std::cout << "haha" << std::endl;

    for (int j=0;j!=nbin;j++) {
	if( signalsBool.at(j) == 1 ) {
	    signals.push_back(j);
	}
    }
  
    // std::cout << "haha" << " " << signals.size() << std::endl;

  
    for (size_t j=0;j!=signals.size();j++) {
	int bin = signals.at(j);
	int prev_bin=bin;
	int next_bin=bin;
    
	int flag = 1;
	while(flag) {
	    prev_bin--;
	    if (find(signals.begin(),signals.end(),prev_bin)==signals.end() || prev_bin <=0) {
		flag = 0;
	    }
	}
    
	flag =1;
	while(flag) {
	    next_bin++;
	    if (find(signals.begin(),signals.end(),next_bin)==signals.end() || next_bin >=nbin-1) {
		flag = 0;
	    }
	}

	//std::cout << bin << " " << std::endl;
	if (prev_bin <0 ) { prev_bin = 0; }
	if (next_bin > nbin - 1) { next_bin = nbin - 1; }

	float prev_content = medians.at(prev_bin);//h44->GetBinContent(prev_bin+1);
	float next_content = medians.at(next_bin);

	float content = prev_content + (bin - prev_bin)/ (next_bin - prev_bin*1.0) 
	    * (next_content - prev_content);
	medians.at(bin) = content;
	//h44->SetBinContent(bin+1,content);
    }
  
    //std::cout << "haha" << std::endl;

    return true;
}

bool Microboone::NoisyFilterAlg(WireCell::Waveform::realseq_t& sig, float min_rms, float max_rms)
{
    const double rmsVal = Microboone::CalcRMSWithFlags(sig);
    
    // int planeNum,channel_no;
    // if (ch < 2400) {
    // 	planeNum = 0;
    // 	channel_no = ch;
    // }
    // else if (ch < 4800) {
    // 	planeNum = 1;
    // 	channel_no = ch - 2400;
    // }
    // else {
    // 	planeNum = 2;
    // 	channel_no = ch - 4800;
    // }

    // //std::cout << rmsVal << " " << ch << std::endl;
  
    // double maxRMSCut[3] = {10.0,10.0,5.0};
    // double minRMSCut[3] = {2,2,2};

    // if (planeNum == 0) {
    // 	if (channel_no < 100) {
    // 	    maxRMSCut[0] = 5;
    // 	    minRMSCut[0] = 1;
    // 	}
    // 	else if (channel_no >= 100 && channel_no<2000) {
    // 	    maxRMSCut[0] = 11; // increase the threshold slightly ... 
    // 	    minRMSCut[0] = 1.9;
    // 	}
    // 	else if (channel_no >= 2000 && channel_no < 2400) {
    // 	    maxRMSCut[0] = 5;
    // 	    minRMSCut[0] = 0.9; // take into account FT-1 channel ... 
    // 	}
    // }
    // else if (planeNum == 1) {
    // 	if (channel_no <290){
    // 	    maxRMSCut[1] = 5;
    // 	    minRMSCut[1] = 1;
    // 	}
    // 	else if (channel_no>=290 && channel_no < 2200) {
    // 	    maxRMSCut[1] = 11;
    // 	    minRMSCut[1] = 1.9;
    // 	}
    // 	else if (channel_no >=2200) {
    // 	    maxRMSCut[1] = 5;
    // 	    minRMSCut[1] = 1;
    // 	}
    // }
    // else if (planeNum == 2) {
    // 	maxRMSCut[2] = 8;
    // 	minRMSCut[2] = 1.3; // reduce threshold to take into account the adaptive baseline ... 
    // }
  
    //if(rmsVal > maxRMSCut[planeNum] || rmsVal < minRMSCut[planeNum]) {
    if(rmsVal > max_rms || rmsVal < min_rms) {
	int numBins = sig.size();
	for(int i = 0; i < numBins; i++) {
	    sig.at(i) = 10000.0;
	}
      
	return true;
    }
  
    return false;
}



bool Microboone::Chirp_raise_baseline(WireCell::Waveform::realseq_t& sig, int bin1, int bin2)
{
    if (bin1 < 0 ) bin1 = 0;
    if (bin2 > (int)sig.size()) bin2 = sig.size();
    for (int i=bin1; i<bin2;i++) {
	sig.at(i) = 10000.0;
    }
    return true;
}

float Microboone::CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig)
{
    float theRMS = 0.0;
    //int waveformSize = sig.size();
  
    WireCell::Waveform::realseq_t temp;
    for (size_t i=0;i!=sig.size();i++){
	if (sig.at(i) < 4096) temp.push_back(sig.at(i));
    }
    float par[3];
    if (temp.size()>0) {
	par[0] = WireCell::Waveform::percentile_binned(temp,0.5 - 0.34);
	par[1] = WireCell::Waveform::percentile_binned(temp,0.5);
	par[2] = WireCell::Waveform::percentile_binned(temp,0.5 + 0.34);
    
	//    std::cout << par[0] << " " << par[1] << " " << par[2] << std::endl;

	theRMS = sqrt((pow(par[2]-par[1],2)+pow(par[1]-par[0],2))/2.);
    }
  
    return theRMS;
}

bool Microboone::SignalFilter(WireCell::Waveform::realseq_t& sig)
{
    const double sigFactor = 4.0;
    const int padBins = 8;
  
    float rmsVal = Microboone::CalcRMSWithFlags(sig);
    float sigThreshold = sigFactor*rmsVal;
  
    float ADCval;
    std::vector<bool> signalRegions;
    int numBins = sig.size();

    for (int i = 0; i < numBins; i++) {
	ADCval = sig.at(i);
	if (((ADCval > sigThreshold) || (ADCval < -1.0*sigThreshold)) && (ADCval < 4096.0)) {
	    signalRegions.push_back(true);
	}
	else {
	    signalRegions.push_back(false);
	}
    }
  
    for(int i = 0; i < numBins; i++) {
	if(signalRegions[i] == true) {
	    int bin1 = i - padBins;
	    if (bin1 < 0 ) {
		bin1 = 0;
	    }
	    int bin2 = i + padBins;
	    if (bin2 > numBins) {
		bin2 = numBins;
	    }
	  
	    for(int j = bin1; j < bin2; j++) {
		ADCval = sig.at(j);
		if(ADCval < 4096.0) {
		    sig.at(j) = sig.at(j)+20000.0;
		}
	    }
	}
    }  
    return true;
}


bool Microboone::RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig)
{
    const int windowSize = 20;
    const int numBins = sig.size();
    int minWindowBins = windowSize/2;
  
    std::vector<double> baselineVec(numBins, 0.0);
    std::vector<bool> isFilledVec(numBins, false);

    int numFlaggedBins = 0;
  
    for(int j = 0; j < numBins; j++) {
	if(sig.at(j) == 10000.0) {
	    numFlaggedBins++;
	}
    }
    if(numFlaggedBins == numBins) {
	return true; // Eventually replace this with flag check
    }
  
    double baselineVal = 0.0;
    int windowBins = 0;
    //int index;
    double ADCval = 0.0;
    for(int j = 0; j <= windowSize/2; j++) {
	ADCval = sig.at(j);
	if(ADCval < 4096.0) {
	    baselineVal += ADCval;
	    windowBins++;
	}
    }

    if(windowBins == 0) {
	baselineVec[0] = 0.0;
    }
    else {
	baselineVec[0] = baselineVal/((double) windowBins);
    }
  
    if(windowBins < minWindowBins) {
	isFilledVec[0] = false;
    }
    else {
	isFilledVec[0] = true;
    }
  
    for(int j = 1; j < numBins; j++) {
	int oldIndex = j-windowSize/2-1;
	int newIndex = j+windowSize/2;
	
	if(oldIndex >= 0) {
	    ADCval = sig.at(oldIndex);
	    if(ADCval < 4096.0) {
		baselineVal -= sig.at(oldIndex);
		windowBins--;
	    }
	}
	if(newIndex < numBins) {  
	    ADCval = sig.at(newIndex);
	    if(ADCval < 4096) {
		baselineVal += sig.at(newIndex);
		windowBins++;
	    }
	}
      
	if(windowBins == 0) {
	    baselineVec[j] = 0.0;
	}
	else {
	    baselineVec[j] = baselineVal/windowBins;
	}
      
	if(windowBins < minWindowBins) {
	    isFilledVec[j] = false;
	}
	else {
	    isFilledVec[j] = true;
	}
    }
  
    for(int j = 0; j < numBins; j++) {
	bool downFlag = false;
	bool upFlag = false;
      
	ADCval = sig.at(j);
	if(ADCval != 10000.0) {
	    if(isFilledVec[j] == false) {
		int downIndex = j;
		while((isFilledVec[downIndex] == false) && (downIndex > 0) && (sig.at(downIndex) != 10000.0)) {
		    downIndex--;
		}
	      
		if(isFilledVec[downIndex] == false) { 
		    downFlag = true;
		}
	      
		int upIndex = j;
		while((isFilledVec[upIndex] == false) && (upIndex < numBins-1) && (sig.at(upIndex) != 10000.0)) {
		    upIndex++;
		}
	      
		if(isFilledVec[upIndex] == false) {
		    upFlag = true;
		}
	      
		if((downFlag == false) && (upFlag == false)) {
		    baselineVec[j] = ((j-downIndex)*baselineVec[downIndex]+(upIndex-j)*baselineVec[upIndex])/((double) upIndex-downIndex);
		}
		else if((downFlag == true) && (upFlag == false)) {
		    baselineVec[j] = baselineVec[upIndex];
		}
		else if((downFlag == false) && (upFlag == true)) {
		    baselineVec[j] = baselineVec[downIndex];
		}
		else {
		    baselineVec[j] = 0.0;
		}
	    }
	  
	    sig.at(j) = ADCval -baselineVec[j];
	}
    }

    return true;
  
}


bool Microboone::RemoveFilterFlags(WireCell::Waveform::realseq_t& sig)
{
    int numBins = sig.size();
    for(int i = 0; i < numBins; i++) {
	double ADCval = sig.at(i);
	if (ADCval > 4096.0) {
	    if(ADCval > 10000.0) {
		sig.at(i) = ADCval-20000.0;
	    }
	    else {
		sig.at(i) = 0.0;
	    }

	}
    }
  
    return true;
}



Microboone::CoherentNoiseSub::CoherentNoiseSub()
{
}
Microboone::CoherentNoiseSub::~CoherentNoiseSub()
{
}

WireCell::Waveform::ChannelMaskMap
Microboone::CoherentNoiseSub::apply(channel_signals_t& chansig) const
{
    //std::cout << "Xin2: " << std::endl;
    // find the median among all 
    WireCell::Waveform::realseq_t medians = Derivations::CalcMedian(chansig);
    //std::cout << medians.size() << " " << medians.at(0) << " " << medians.at(1) << std::endl;
  

    // For Xin: here is how you can get the response spectrum for this group.
    const int achannel = chansig.begin()->first;

    //std::cout << m_noisedb->response_offset(achannel) << std::endl;

    const Waveform::compseq_t& respec = m_noisedb->response(achannel);
    const int res_offset = m_noisedb->response_offset(achannel);
    const int pad_f = m_noisedb->pad_window_front(achannel);
    const int pad_b = m_noisedb->pad_window_back(achannel);


    //    if (respec.size()) {
    // now, apply the response spectrum to deconvolve the median
    // and apply the special protection or pass respec into
    // SignalProtection().
    //}

    // do the signal protection and adaptive baseline
    Microboone::SignalProtection(medians,respec,res_offset,pad_f,pad_b);
    
    //std::cout <<"abc " << " " << chansig.size() << " " << medians.size() << std::endl;

    // calculate the scaling coefficient and subtract
    Microboone::Subtract_WScaling(chansig, medians);

    //std::cout <<"abc1 " << std::endl;

    // for (auto it: chansig){
    //   std::cout << "Xin3 " << it.second.at(0) << std::endl;
    // }
  
    return WireCell::Waveform::ChannelMaskMap();		// not implemented
}
WireCell::Waveform::ChannelMaskMap
Microboone::CoherentNoiseSub::apply(int channel, signal_t& sig) const
{
    return WireCell::Waveform::ChannelMaskMap();		// not implemented
}




Microboone::OneChannelNoise::OneChannelNoise()
    : m_check_chirp() // fixme, there are magic numbers hidden here
    , m_check_partial() // fixme, here too.
{
}
Microboone::OneChannelNoise::~OneChannelNoise()
{
}

void Microboone::OneChannelNoise::configure(const WireCell::Configuration& config)
{
    // fixme!
}
WireCell::Configuration Microboone::OneChannelNoise::default_configuration() const
{
    // fixme!
    Configuration cfg;
    return cfg;
}


WireCell::Waveform::ChannelMaskMap Microboone::OneChannelNoise::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    // fixme: some channels are just bad can should be skipped.

    // get signal with nominal baseline correction
    float baseline = m_noisedb->nominal_baseline(ch);
    WireCell::Waveform::increase(signal, baseline *(-1));

    // get signal with nominal gain correction 
    float gc = m_noisedb->gain_correction(ch);
    // if (ch < 2400)
    //   std::cout << gc << " " << ch << std::endl;
    auto signal_gc = signal; // copy, need to keep original signal
    
    WireCell::Waveform::scale(signal_gc, gc);
    
    // determine if chirping
    WireCell::Waveform::BinRange chirped_bins;
    bool is_chirp = m_check_chirp(signal_gc, chirped_bins.first, chirped_bins.second);
    if (is_chirp) {
      ret["chirp"][ch].push_back(chirped_bins);
      // for (int i=chirped_bins.first;i!=chirped_bins.second;i++){
      // 	signal.at(i) = 0;
      // }
    }

    auto spectrum = WireCell::Waveform::dft(signal);
    bool is_partial = m_check_partial(spectrum); // Xin's "IS_RC()"
    // if (is_partial){
    //   std::cout << ch << std::endl;
    // }
    
    if (!is_partial) {
	//std::cout << "Xin: " << spectrum.front().real() << " " ;
	WireCell::Waveform::shrink(spectrum, m_noisedb->rcrc(ch));
	//std::cout << spectrum.front().real() << std::endl;
    }

    // if (ch==2000) std::cout << "2000" << " " << m_noisedb->config(ch).at(1) << " " << m_noisedb->gain_correction(ch) << std::endl;
    // if (ch==2016) std::cout << "2016" << " " << m_noisedb->config(ch).at(1) << " " << m_noisedb->gain_correction(ch) << std::endl;

   
    WireCell::Waveform::scale(spectrum, m_noisedb->config(ch));
    
    WireCell::Waveform::scale(spectrum, m_noisedb->noise(ch));

    // remove the DC component 
    spectrum.front() = 0;
    signal = WireCell::Waveform::idft(spectrum);

    //Now calculate the baseline ...
    baseline = WireCell::Waveform::median_binned(signal);
    //correct baseline
    WireCell::Waveform::increase(signal, baseline *(-1));


    //*** from here down is where things become way to microboone
    //*** specific and are basically copy-pasted from the prototype


    // Now do adaptive baseline for the chirping channels
    if (is_chirp) {
	Microboone::Chirp_raise_baseline(signal,chirped_bins.first, chirped_bins.second);
	Microboone::SignalFilter(signal);
	Microboone::RawAdapativeBaselineAlg(signal);
    }
    // Now do the adaptive baseline for the bad RC channels
    if (is_partial) {
	Microboone::SignalFilter(signal);
	Microboone::RawAdapativeBaselineAlg(signal);
    }

    // Identify the Noisy channels ... 
    Microboone::SignalFilter(signal);

    //
    const float min_rms = m_noisedb->min_rms_cut(ch);
    const float max_rms = m_noisedb->max_rms_cut(ch);

    bool is_noisy = Microboone::NoisyFilterAlg(signal,min_rms,max_rms);
    Microboone::RemoveFilterFlags(signal);

    // if (is_noisy){
    //   std::cout << "Xin: " << signal.at(1) << std::endl;
    // }

    // std::cout << ch << " " << is_chirp << " " << is_partial << " " << is_noisy << std::endl;

    if (is_noisy) {
	chirped_bins.first = 0;
	chirped_bins.second = signal.size();
	ret["noisy"][ch].push_back(chirped_bins);
    }

    return ret;
}



WireCell::Waveform::ChannelMaskMap Microboone::OneChannelNoise::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}






Microboone::ADCBitShift::ADCBitShift(int nbits, int exam_nticks, double threshold_sigma,  double threshold_fix)
    : m_nbits(nbits)
    , m_exam_nticks(exam_nticks)
    , m_threshold_sigma(threshold_sigma)
    , m_threshold_fix(threshold_fix)
{
}
Microboone::ADCBitShift::~ADCBitShift()
{
}

void Microboone::ADCBitShift::configure(const WireCell::Configuration& cfg)
{
    m_nbits = get<int>(cfg, "Number_of_ADC_bits", m_nbits);
    m_exam_nticks = get<int>(cfg, "Exam_number_of_ticks_test", m_exam_nticks);
    m_threshold_sigma = get<double>(cfg, "Threshold_sigma_test", m_threshold_sigma);
    m_threshold_fix = get<double>(cfg, "Threshold_fix", m_threshold_fix);
}
WireCell::Configuration Microboone::ADCBitShift::default_configuration() const
{
    Configuration cfg;
    cfg["Number_of_ADC_bits"] = m_nbits;
    cfg["Exam_number_of_ticks_test"] = m_exam_nticks;
    cfg["Threshold_sigma_test"] = m_threshold_sigma;
    cfg["Threshold_fix"] = m_threshold_fix;
	
    return cfg;
}


// Not useful ... 

WireCell::Waveform::ChannelMaskMap Microboone::ADCBitShift::apply(channel_signals_t& chansig) const
{

    
    return WireCell::Waveform::ChannelMaskMap();
}

// ADC Bit Shift problem ... 
WireCell::Waveform::ChannelMaskMap Microboone::ADCBitShift::apply(int ch, signal_t& signal) const
{
    std::vector<int> counter(m_nbits,0);
    std::vector<int> s,s1(m_nbits,0);
    int nbin = m_exam_nticks;

    for (int i=0;i!=nbin;i++){
	int x = signal.at(i);
	s.clear();
	do
	    {
		s.push_back( (x & 1));
	    } while (x >>= 1);
	s.resize(m_nbits);
	
	for (int j=0;j!=m_nbits;j++){
	    counter.at(j) += abs(s.at(j) - s1.at(j));
	}
	s1=s;
    }
    
    int nshift = 0;
    for (int i=0;i!=m_nbits;i++){
	if (counter.at(i) < nbin/2. - nbin/2. *sqrt(1./nbin) * m_threshold_sigma){
	    nshift ++;
	}else{
	    break;
	}
    }
    
    WireCell::Waveform::ChannelMaskMap ret;
    if (nshift!=0 && nshift < 11){
	WireCell::Waveform::BinRange ADC_bit_shifts;
	ADC_bit_shifts.first = nshift;
	ADC_bit_shifts.second = nshift;
	
	ret["ADCBitShift"][ch].push_back(ADC_bit_shifts);

	// do the correction ...
	const int nl = signal.size();
	int x[nl], x_orig[nl];
	for (int i=0;i!=nl;i++){
	    x_orig[i] = signal.at(i);
	    x[i] = signal.at(i);
	}
	
	std::set<int> fillings;
	int filling;
	int mean = 0;
	for (int i=0;i!=nl;i++){
	    filling = WireCell::Bits::lowest_bits(x_orig[i], nshift);
	    int y = WireCell::Bits::shift_right(x_orig[i], nshift, filling, 12);
	    fillings.insert(filling);
	    // cout << y << " ";
	    // filling = lowest_bits(x[i], nshift);
	    // filling = x[i] & int(pow(2, nshift)-1);
	    x[i] = y;
	    mean += x[i];
	}
	mean = mean/nl;

	int exp_diff = pow(2,m_nbits-nshift)*m_threshold_fix;

	// examine the results ... 
	int prev1_bin_content = mean;
	int prev_bin_content = mean;
	int next_bin_content = mean;
	int next1_bin_content = mean;
	int curr_bin_content;
	
	for (int i=0;i<nl;i=i+1){
	    curr_bin_content = x[i];
	    // when to judge if one is likely bad ... 
	    if (abs(curr_bin_content-mean)>exp_diff && 
		(abs(curr_bin_content - prev_bin_content) > exp_diff ||
		 abs(curr_bin_content - next_bin_content) > exp_diff)
		){
		int exp_value = ( (2*prev_bin_content - prev1_bin_content) +
				    (prev_bin_content + next_bin_content)/2. + 
				    (prev_bin_content * 2./3. + next1_bin_content/3.))/3.;
		for (auto it = fillings.begin(); it!=fillings.end(); it++){
		    int y = WireCell::Bits::shift_right(x_orig[i], nshift, (*it), 12);
		    // when to switch ... 
		    if (fabs(y-exp_value) < fabs(x[i] - exp_value)){
			x[i] = y;//hs->SetBinContent(i+1,y);
		    }
		}
	    }
	    // if (chid == 1281 && i < 10){
	    // 	std::cout <<x_orig[i] << " " << x[i] << std::endl;
	    // }
	    //hs->SetBinContent(i+1,mean);
	    
	    prev1_bin_content = prev_bin_content;
	    prev_bin_content = x[i];
	    if (i+2 < nl){
		next_bin_content = x[i+2];
	    }else{
		next_bin_content = mean;
	    }
	    if (i+3 < nl){
		next1_bin_content = x[i+3];
	    }else{
		next1_bin_content = mean;
	    }
	}
	
	
	// change the histogram ...
	for (int i=0;i!=nl;i++){
	    signal.at(i) = x[i];
	}
	
	
    }
   
    
    
    return ret;
}



// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
