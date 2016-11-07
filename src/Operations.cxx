/** FIXME: this file is full of magic numbers and likely totally not
 * usable for detectors other than MicroBooNE. 
 */

#include "WireCellSigProc/Operations.h"
#include "WireCellSigProc/Derivations.h"

#include <cmath>
#include <complex>
#include <iostream>

using namespace WireCellSigProc;

bool Operations::Subtract_WScaling(WireCell::IChannelFilter::channel_signals_t& chansig,
				   const WireCell::Waveform::realseq_t& medians)
{
  double ave_coef = 0;
  double_t ave_coef1 = 0;
  std::map<int,double> coef_all;

  const int nbin = medians.size();

  for (auto it: chansig){
    int ch = it.first;
    WireCell::IChannelFilter::signal_t& signal = it.second;
    double sum2 = 0;
    double sum3 = 0;
    double coef = 0;
    std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);
    
    for (int j=0;j!=nbin;j++){
      if (fabs(signal.at(j)) < 4 * temp.second){
	sum2 += signal.at(j) * medians.at(j);
	sum3 += medians.at(j) * medians.at(j);
      }
    }
    if (sum3 >0) {
      coef = sum2/sum3;
    }
    coef_all[ch] = coef;
    if (coef != 0){
      ave_coef += coef;
      ave_coef1 ++;
    }
  }
  if (ave_coef1>0) {
    ave_coef = ave_coef / ave_coef1;
  }
  
  for (auto it: chansig){
    int ch = it.first;
    WireCell::IChannelFilter::signal_t& signal = it.second;
    float scaling;
    if (ave_coef != 0 ){
      scaling = coef_all[ch]/ave_coef;
    }else{
      scaling = 0;
    }
    //std::cout << scaling << " " << signal.at(0) << " " << medians.at(0) << std::endl;
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

bool Operations::SignalProtection(WireCell::Waveform::realseq_t& medians){
  // calculate the RMS 
  std::pair<double,double> temp = Derivations::CalcRMS(medians);
  double mean = temp.first;
  double rms = temp.second;
  
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

  int pad_window = 5;
  int protection_factor = 5.0;

  std::vector<int> signals;
  std::vector<bool> signalsBool;
  for (int j=0;j!=nbin;j++)
    signalsBool.push_back(0);
  
  for (int j=0;j!=nbin;j++){
    float content = medians.at(j);
    if (fabs(content-mean)>protection_factor*rms){
      medians.at(j) = 0; 
      signalsBool.at(j) = 1;
      // add the front and back padding
      for (int k=0;k!=pad_window;k++){
	int bin = j+k+1;
	if (bin > nbin-1) bin = nbin-1;
	signalsBool.at(bin) = 1;
	bin = j-k-1;
	if (bin <0) { bin = 0; }
	signalsBool.at(bin) = 1;
      }
    }
  }
  
  // std::cout << "haha" << std::endl;

  for (int j=0;j!=nbin;j++)
    if( signalsBool.at(j) == 1 )
      signals.push_back(j);
  
  // std::cout << "haha" << " " << signals.size() << std::endl;

  
  for (int j=0;j!=signals.size();j++){
    int bin = signals.at(j);
    int prev_bin=bin;
    int next_bin=bin;
    
    int flag = 1;
    while(flag){
      prev_bin--;
      if (find(signals.begin(),signals.end(),prev_bin)==signals.end() || prev_bin <=0){
	flag = 0;
      }
    }
    
    flag =1;
    while(flag){
      next_bin++;
      if (find(signals.begin(),signals.end(),next_bin)==signals.end() || next_bin >=nbin-1){
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

bool Operations::NoisyFilterAlg(WireCell::Waveform::realseq_t& sig, int ch){
  const double rmsVal = Operations::CalcRMSWithFlags(sig);
  int planeNum,channel_no;
  if (ch < 2400){
    planeNum = 0;
    channel_no = ch;
  }else if (ch < 4800){
    planeNum = 1;
    channel_no = ch - 2400;
  }else{
    planeNum = 2;
    channel_no = ch - 4800;
  }

  //std::cout << rmsVal << " " << ch << std::endl;
  
  double maxRMSCut[3] = {10.0,10.0,5.0};
  double minRMSCut[3] = {2,2,2};

  if (planeNum == 0){
    if (channel_no < 100){
      maxRMSCut[0] = 5;
      minRMSCut[0] = 1;
    }else if (channel_no >= 100 && channel_no<2000){
      maxRMSCut[0] = 11; // increase the threshold slightly ... 
      minRMSCut[0] = 1.9;
    }else if (channel_no >= 2000 && channel_no < 2400){
      maxRMSCut[0] = 5;
      minRMSCut[0] = 0.9; // take into account FT-1 channel ... 
    }
  }else if (planeNum == 1){
    if (channel_no <290){
      maxRMSCut[1] = 5;
      minRMSCut[1] = 1;
    }else if (channel_no>=290 && channel_no < 2200){
      maxRMSCut[1] = 11;
      minRMSCut[1] = 1.9;
    }else if (channel_no >=2200){
      maxRMSCut[1] = 5;
      minRMSCut[1] = 1;
    }
  }else if (planeNum == 2){
    maxRMSCut[2] = 8;
    minRMSCut[2] = 1.3; // reduce threshold to take into account the adaptive baseline ... 
  }
  
  if(rmsVal > maxRMSCut[planeNum] || rmsVal < minRMSCut[planeNum])
    {
      int numBins = sig.size();
      for(int i = 0; i < numBins; i++)
	{
	  sig.at(i) = 10000.0;
	}         
      
      return true;
    }
  
  return false;
}



bool Operations::Chirp_raise_baseline(WireCell::Waveform::realseq_t& sig, int bin1, int bin2){
  if (bin1 < 0 ) bin1 = 0;
  if (bin2 > sig.size()) bin2 = sig.size();
  for (int i=bin1; i<bin2;i++){
    sig.at(i) = 10000.0;
  }
  return true;
}

float Operations::CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig){
  float theRMS = 0.0;
  int waveformSize = sig.size();
  
  WireCell::Waveform::realseq_t temp;
  for (int i=0;i!=sig.size();i++){
    if (sig.at(i) < 4096) temp.push_back(sig.at(i));
  }
  float par[3];
  if (temp.size()>0){
    par[0] = WireCell::Waveform::percentile(temp,0.5 - 0.34);
    par[1] = WireCell::Waveform::percentile(temp,0.5);
    par[2] = WireCell::Waveform::percentile(temp,0.5 + 0.34);
    
    //    std::cout << par[0] << " " << par[1] << " " << par[2] << std::endl;

    theRMS = sqrt((pow(par[2]-par[1],2)+pow(par[1]-par[0],2))/2.);
  }
  
  return theRMS;
}

bool Operations::SignalFilter(WireCell::Waveform::realseq_t& sig){
  const double sigFactor = 4.0;
  const int padBins = 8;
  
  float rmsVal = Operations::CalcRMSWithFlags(sig);
  float sigThreshold = sigFactor*rmsVal;
  
  float ADCval;
  std::vector<bool> signalRegions;
  int numBins = sig.size();

  for (int i = 0; i < numBins; i++){
    ADCval = sig.at(i);
    if (((ADCval > sigThreshold) || (ADCval < -1.0*sigThreshold)) && (ADCval < 4096.0)){
      signalRegions.push_back(true);
    }else{
      signalRegions.push_back(false);
    }
  }
  
  for(int i = 0; i < numBins; i++)
    {
      if(signalRegions[i] == true)
	{
	  int bin1 = i - padBins;
	  if (bin1 < 0 ) bin1 = 0;
	  int bin2 = i + padBins;
	  if (bin2 > numBins) bin2 = numBins;
	  
	  for(int j = bin1; j < bin2; j++)
	    {
	      ADCval = sig.at(j);
	      if(ADCval < 4096.0)
		{
		  sig.at(j) = sig.at(j)+20000.0;
		}
	    }
	}
    }  
  return true;
}


bool Operations::RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig){
  const int windowSize = 20;
  const int numBins = sig.size();
  int minWindowBins = windowSize/2;
  
  double baselineVec[numBins];
  bool isFilledVec[numBins];

  int numFlaggedBins = 0;
  
  for(int j = 0; j < numBins; j++)
    {
      if(sig.at(j) == 10000.0)
	{
	  numFlaggedBins++;
	}
    }
  if(numFlaggedBins == numBins) return true; // Eventually replace this with flag check
  
  double baselineVal = 0.0;
  int windowBins = 0;
  int index;
  double ADCval;
  for(int j = 0; j <= windowSize/2; j++)
    {
      ADCval = sig.at(j);
      if(ADCval < 4096.0)
	{
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
  
  int oldIndex;
  int newIndex;
  for(int j = 1; j < numBins; j++)
    {
      oldIndex = j-windowSize/2-1;
      newIndex = j+windowSize/2;
      
      if(oldIndex >= 0)
	{
	  ADCval = sig.at(oldIndex);
	  if(ADCval < 4096.0)
	    {
	      baselineVal -= sig.at(oldIndex);
	      windowBins--;
	    }
	}
      if(newIndex < numBins)
	{  
	  ADCval = sig.at(newIndex);
	  if(ADCval < 4096)
	    {
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
  
  int downIndex;
  int upIndex;
  bool downFlag;
  bool upFlag;
  for(int j = 0; j < numBins; j++)
    {
      downFlag = false;
      upFlag = false;
      
      ADCval = sig.at(j);
      if(ADCval != 10000.0)
	{
	  if(isFilledVec[j] == false)
	    {
	      downIndex = j;
	      while((isFilledVec[downIndex] == false) && (downIndex > 0) && (sig.at(downIndex) != 10000.0))
		{
		  downIndex--;
		}
	      
	      if(isFilledVec[downIndex] == false)
		downFlag = true;
	      
	      upIndex = j;
	      while((isFilledVec[upIndex] == false) && (upIndex < numBins-1) && (sig.at(upIndex) != 10000.0))
		{
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


bool Operations::RemoveFilterFlags(WireCell::Waveform::realseq_t& sig){
  double ADCval;
  int numBins = sig.size();
  for(int i = 0; i < numBins; i++)
    {
      ADCval = sig.at(i);
      if (ADCval > 4096.0)
	{
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


// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
