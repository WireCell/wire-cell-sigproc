#include "WireCellSigProc/Operations.h"

#include <cmath>
#include <complex>
#include <iostream>

using namespace WireCellSigProc;

bool Operations::NoisyFilterAlg(WireCell::Waveform::realseq_t& sig, int ch){
  double rmsVal = Operations::CalcRMSWithFlags(sig);
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

  //  std::cout << rmsVal << " " << ch << std::endl;
  
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
  if (bin2 > sig.size()-1) bin2 = sig.size()-1;
  for (int i=bin1; i<=bin2;i++){
    sig.at(i) = sig.at(i) + 10000.0;
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

  if(windowBins == 0)
    baselineVec[0] = 0.0;
  else
    baselineVec[0] = baselineVal/((double) windowBins);
  
  if(windowBins < minWindowBins)
    isFilledVec[0] = false;
  else
    isFilledVec[0] = true;
  
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
      
      if(windowBins == 0)
	baselineVec[j] = 0.0;
      else
	baselineVec[j] = baselineVal/windowBins;
      
      if(windowBins < minWindowBins)
	isFilledVec[j] = false;
      else
	isFilledVec[j] = true;
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
	      
	      if(isFilledVec[upIndex] == false)
		upFlag = true;
	      
	      if((downFlag == false) && (upFlag == false))
		baselineVec[j] = ((j-downIndex)*baselineVec[downIndex]+(upIndex-j)*baselineVec[upIndex])/((double) upIndex-downIndex);
	      else if((downFlag == true) && (upFlag == false))
		baselineVec[j] = baselineVec[upIndex];
	      else if((downFlag == false) && (upFlag == true))
		baselineVec[j] = baselineVec[downIndex];
	      else
		baselineVec[j] = 0.0;
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
	  if(ADCval > 10000.0)
	    sig.at(i) = ADCval-20000.0;
	  else
	    sig.at(i) = 0.0;

	}
    }
  
  
  return true;
}
