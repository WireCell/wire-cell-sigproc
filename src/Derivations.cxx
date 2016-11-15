#include "WireCellSigProc/Derivations.h"

#include <iostream>
using namespace WireCellSigProc;

std::pair<double,double> Derivations::CalcRMS(const WireCell::Waveform::realseq_t& signal)
{
    std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);
    double mean = temp.first;
    double rms = temp.second;
    WireCell::Waveform::realseq_t temp1;
    for (int i=0;i!=signal.size();i++){
	if (fabs(signal.at(i)-mean) < 4.5 * rms){
	    temp1.push_back(signal.at(i));
	}
    }
    temp = WireCell::Waveform::mean_rms(temp1);
    return temp;
}

WireCell::Waveform::realseq_t Derivations::CalcMedian(const WireCell::IChannelFilter::channel_signals_t& chansig)
{
    float max_rms = 0;

    //std::cout << "Xin3: " << chansig.size() << std::endl;
    const int nchannel = chansig.size();
    const int nbins = (chansig.begin()->second).size();
    float content[nchannel][nbins];
    int start_ch = 0;
    for (auto it: chansig){
     	int ch = it.first;
     	WireCell::IChannelFilter::signal_t& signal = it.second;
    	std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);
    	if (temp.second > max_rms) {
    	    max_rms = temp.second;
    	}
	for (int i=0;i!=nbins;i++){
	    content[start_ch][i] = signal.at(i);
	}
	start_ch ++;
    	//std::cout << ch << " " << signal.size() << std::endl;
	//nbins = signal.size();
	//std::cout << ch << " " << signal.size() << std::endl;
    }
    //std::cout << max_rms << std::endl;
	
    
    
  

    // std::cout << nbins << " " << chansig.size() << std::endl;
    WireCell::Waveform::realseq_t medians(nbins);//(chansig.size());
    for (int i=0;i!=nbins;i++){
	WireCell::Waveform::realseq_t temp;
	for (int j=0;j!=nchannel;j++){
	    if (fabs(content[j][i]) < 5 * max_rms && 
	 	fabs(content[j][i]) > 0.001){
	 	temp.push_back(content[j][i]);
	    } 
	}
// for (auto it: chansig){
	//     WireCell::IChannelFilter::signal_t& signal = it.second;
	//     // std::cout << signal.at(i) << " " << max_rms << std::endl;
	//     if (fabs(signal.at(i)) < 5 * max_rms && 
	// 	fabs(signal.at(i)) > 0.001){
	// 	temp.push_back(signal.at(i));
	//     } 
	//     // if (i==10){
	//     // 	std::cout << it.first << " " << signal.at(i) << std::endl;
	//     // }

	// }
	// if (i==10)
	//   std::cout << temp.size() << std::endl;
	if (temp.size()>0){
	    medians.at(i)=WireCell::Waveform::median_binned(temp);
	}else{
	    medians.at(i)=0;
	}
    }

    return medians;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
