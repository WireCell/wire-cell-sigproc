
#include "ROI_formation.h"

#include <iostream>

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_formation::ROI_formation(Waveform::ChannelMaskMap& cmm,int nwire_u, int nwire_v, int nwire_w, int nbins, float th_factor_ind, float th_factor_col, int pad, float asy, int rebin , double l_factor, double l_max_th, double l_factor1, int l_short_length , double l_fixed_threshold )
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
  , nbins(nbins)
  , th_factor_ind(th_factor_ind)
  , th_factor_col(th_factor_col)
  , pad(pad)
  , asy(asy)
  , rebin(rebin)
  , l_factor(l_factor)
  , l_max_th(l_max_th)
  , l_factor1(l_factor1)
  , l_short_length(l_short_length)
  , l_fixed_threshold(l_fixed_threshold)
{
  self_rois_u.resize(nwire_u);
  self_rois_v.resize(nwire_v);
  self_rois_w.resize(nwire_w);
  
  loose_rois_u.resize(nwire_u);
  loose_rois_v.resize(nwire_v);
  loose_rois_w.resize(nwire_w);

  uplane_rms.resize(nwire_u);
  vplane_rms.resize(nwire_v);
  wplane_rms.resize(nwire_w);

  for (auto it = cmm["bad"].begin(); it!=cmm["bad"].end(); it++){
    int ch = it->first;
    std::vector<std::pair<int,int>> temps;
    bad_ch_map[ch] = temps;
    for (size_t ind = 0; ind < it->second.size(); ind++){
      bad_ch_map[ch].push_back(std::make_pair(it->second.at(ind).first, it->second.at(ind).second));
      //std::cout << ch << " " <<  << std::endl;
    }
  }
  //std::cout << bad_ch_map.size() << std::endl;
  
}

ROI_formation::~ROI_formation(){
  
}

void ROI_formation::Clear(){
  self_rois_u.clear();
  self_rois_v.clear();
  self_rois_w.clear();
  
  loose_rois_u.clear();
  loose_rois_v.clear();
  loose_rois_w.clear();

  uplane_rms.clear();
  vplane_rms.clear();
  wplane_rms.clear();
}




void ROI_formation::extend_ROI_self(int plane){
  if (plane==0){
    for (size_t i=0;i!=self_rois_u.size();i++){
      std::vector<std::pair<int,int>> temp_rois;
      int temp_begin=0, temp_end=0;
      for (size_t j=0;j!=self_rois_u.at(i).size();j++){
	temp_begin = self_rois_u.at(i).at(j).first - pad;
	if (temp_begin < 0 ) temp_begin = 0;
	temp_end = self_rois_u.at(i).at(j).second + pad;
	if (temp_end >= nbins) temp_end = nbins - 1;
	// merge 
	if (temp_rois.size() == 0){
	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	}else{
	  if (temp_begin < temp_rois.back().second){
	    if (temp_end > temp_rois.back().second){
	      temp_rois.back().second = temp_end;
	    }
	  }else{
	    temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	  }
	}
      }
      self_rois_u.at(i) = temp_rois;
    }
  }else if (plane==1){
    for (size_t i=0;i!=self_rois_v.size();i++){
      std::vector<std::pair<int,int>> temp_rois;
      int temp_begin=0, temp_end=0;
      for (size_t j=0;j!=self_rois_v.at(i).size();j++){
	temp_begin = self_rois_v.at(i).at(j).first - pad;
	if (temp_begin < 0 ) temp_begin = 0;
	temp_end = self_rois_v.at(i).at(j).second + pad;
	if (temp_end >= nbins) temp_end = nbins - 1;
	// merge 
	if (temp_rois.size() == 0){
	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	}else{
	  if (temp_begin < temp_rois.back().second){
	    if (temp_end > temp_rois.back().second){
	      temp_rois.back().second = temp_end;
	    }
	  }else{
	    temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	  }
	}
      }
      self_rois_v.at(i) = temp_rois;
    }
  }else if (plane==2){
    for (size_t i=0;i!=self_rois_w.size();i++){
      std::vector<std::pair<int,int>> temp_rois;
      int temp_begin=0, temp_end=0;
      for (size_t j=0;j!=self_rois_w.at(i).size();j++){
	temp_begin = self_rois_w.at(i).at(j).first - pad;
	if (temp_begin < 0 ) temp_begin = 0;
	temp_end = self_rois_w.at(i).at(j).second + pad;
	if (temp_end >= nbins) temp_end = nbins - 1;
	// merge 
	if (temp_rois.size() == 0){
	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	}else{
	  if (temp_begin < temp_rois.back().second){
	    if (temp_end > temp_rois.back().second){
	      temp_rois.back().second = temp_end;
	    }
	  }else{
	    temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	  }
	}
      }
      self_rois_w.at(i) = temp_rois;
    }
  }

}

void ROI_formation::create_ROI_connect_info(int plane){

  if (plane==0){
    // u 
    for (int i=0;i!=nwire_u-2;i++){
      for (size_t j=0; j!=self_rois_u.at(i).size();j++){
	int start1 = self_rois_u.at(i).at(j).first;
	int end1 = self_rois_u.at(i).at(j).second;
	int length1 = end1-start1+1;
	for (size_t k=0; k!=self_rois_u.at(i+2).size();k++){
	  int start2 = self_rois_u.at(i+2).at(k).first;
	  int end2 = self_rois_u.at(i+2).at(k).second;
	  int length2 = end2 - start2 + 1;
	  if ( fabs(length2 - length1) < (length2 + length1) * asy){
	    int start3 = (start1+start2)/2.;
	    int end3 = (end1+end2)/2.;
	    if (start3 < end3 && start3 <= end1 && start3 <=end2 && end3 >= start1 && end3 >=start2){
	      // go through existing ones to make sure there is no overlap
	      int flag = 0; 
	      for (size_t i1 = 0; i1!=self_rois_u.at(i+1).size();i1++){
		int max_start = start3;
		if (self_rois_u.at(i+1).at(i1).first > max_start)
		  max_start = self_rois_u.at(i+1).at(i1).first;
		int min_end = end3;
		if (self_rois_u.at(i+1).at(i1).second < min_end)
		  min_end = self_rois_u.at(i+1).at(i1).second ;
		if (max_start < min_end){
		  flag = 1;
		  break;
		}
	      }
	      if (flag == 0)
		self_rois_u.at(i+1).push_back(std::make_pair(start3,end3));
	    }
	  }
	} 
      }
    }
  }else if (plane==1){
    // v
    for (int i=0;i!=nwire_v-2;i++){
      for (size_t j=0; j!=self_rois_v.at(i).size();j++){
	int start1 = self_rois_v.at(i).at(j).first;
	int end1 = self_rois_v.at(i).at(j).second;
	int length1 = end1-start1+1;
	for (size_t k=0; k!=self_rois_v.at(i+2).size();k++){
	  int start2 = self_rois_v.at(i+2).at(k).first;
	  int end2 = self_rois_v.at(i+2).at(k).second;
	  int length2 = end2 - start2 + 1;
	  if ( fabs(length2 - length1) < (length2 + length1) * asy){
	    int start3 = (start1+start2)/2.;
	    int end3 = (end1+end2)/2.;
	    if (start3 < end3 && start3 <= end1 && start3 <=end2 && end3 >= start1 && end3 >=start2){
	      // go through existing ones to make sure there is no overlap
	      int flag = 0; 
	      for (size_t i1 = 0; i1!=self_rois_v.at(i+1).size();i1++){
		int max_start = start3;
		if (self_rois_v.at(i+1).at(i1).first > max_start)
		  max_start = self_rois_v.at(i+1).at(i1).first;
		int min_end = end3;
		if (self_rois_v.at(i+1).at(i1).second < min_end)
		  min_end = self_rois_v.at(i+1).at(i1).second ;
		if (max_start < min_end){
		  flag = 1;
		  break;
		}
	      }
	      if (flag == 0)
		self_rois_v.at(i+1).push_back(std::make_pair(start3,end3));
	    }
	  }
	} 
      }
    }
  }else if (plane==2){
    // w?
    for (int i=0;i!=nwire_w-2;i++){
      for (size_t j=0; j!=self_rois_w.at(i).size();j++){
	int start1 = self_rois_w.at(i).at(j).first;
	int end1 = self_rois_w.at(i).at(j).second;
	int length1 = end1-start1+1;
	for (size_t k=0; k!=self_rois_w.at(i+2).size();k++){
	  int start2 = self_rois_w.at(i+2).at(k).first;
	  int end2 = self_rois_w.at(i+2).at(k).second;
	  int length2 = end2 - start2 + 1;
	  if ( fabs(length2 - length1) < (length2 + length1) * asy){
	    int start3 = (start1+start2)/2.;
	    int end3 = (end1+end2)/2.;
	    if (start3 < end3 && start3 <= end1 && start3 <=end2 && end3 >= start1 && end3 >=start2){
	      // go through existing ones to make sure there is no overlap
	      int flag = 0; 
	      for (size_t i1 = 0; i1!=self_rois_w.at(i+1).size();i1++){
		int max_start = start3;
		if (self_rois_w.at(i+1).at(i1).first > max_start)
		  max_start = self_rois_w.at(i+1).at(i1).first;
		int min_end = end3;
		if (self_rois_w.at(i+1).at(i1).second < min_end)
		  min_end = self_rois_w.at(i+1).at(i1).second ;
		if (max_start < min_end){
		  flag = 1;
		  break;
		}
	      }
	      if (flag == 0)
		self_rois_w.at(i+1).push_back(std::make_pair(start3,end3));
	    }
	  }
	} 
      }
    }
  }
}

double ROI_formation::cal_RMS(Waveform::realseq_t signal){
  double result = 0;
  if (signal.size()>0){
    // do quantile ... 
    float par[3];
    par[0] = WireCell::Waveform::percentile_binned(signal,0.5 - 0.34);
    par[1] = WireCell::Waveform::percentile_binned(signal,0.5);
    par[2] = WireCell::Waveform::percentile_binned(signal,0.5 + 0.34);
    float rms = sqrt((pow(par[2]-par[1],2)+pow(par[1]-par[0],2))/2.);

    float rms2 = 0;
    float rms1 = 0;
    for(size_t i =0; i!=signal.size();i++){
      if (fabs(signal.at(i)) < 5.0 * rms){
	rms1 += pow(signal.at(i),2);
	rms2 ++;
      }
    }
    if (rms2 >0){
      result = sqrt(rms1/rms2);
    }
  }
  
  return result;
}

void ROI_formation::find_ROI_by_decon_itself(int plane, const Array::array_xxf& r_data, const Array::array_xxf& r_data_tight){

  int offset=0;
  if (plane==0){
    offset = 0;
  }else if (plane==1){
    offset = nwire_u;
  }else if (plane==2){
    offset = nwire_u + nwire_v;
  }
  
  for (int irow = 0; irow!=r_data.rows(); irow++){
    // calclulate rms for a row of r_data
    Waveform::realseq_t signal(nbins);
    Waveform::realseq_t signal1(nbins);
    Waveform::realseq_t signal2(nbins);

    for (int icol = 0; icol!=r_data.cols();icol++){
     
      
    }
    
    if (bad_ch_map.find(irow+offset)!=bad_ch_map.end()){
      int ncount = 0;
      for (int icol=0;icol!=r_data.cols();icol++){
	bool flag = true;
	for (size_t i=0; i!=bad_ch_map[irow+offset].size(); i++){
	  if (icol >= bad_ch_map[irow+offset].at(i).first &&
	      icol <= bad_ch_map[irow+offset].at(i).second){
	    flag = false;
	    break;
	  }
	}
	if (flag){
	  signal.at(ncount) = r_data(irow,icol);
	  signal1.at(icol) = r_data(irow,icol);
	  signal2.at(icol) = r_data_tight(irow,icol);
	  ncount ++;
	}else{
	  signal1.at(icol) = 0;
	  signal2.at(icol) = 0;
	}
      }
      signal.resize(ncount);
    }else{
      for (int icol = 0; icol!= r_data.cols(); icol++){
	signal.at(icol) = r_data(irow,icol);
	signal1.at(icol) = r_data(irow,icol);
	signal2.at(icol) = r_data_tight(irow,icol);
      }
    }
    // do threshold and fill rms 
    double rms = cal_RMS(signal);
    double threshold = 0;
    if (plane==0){
      threshold = th_factor_ind * rms + 1;
      uplane_rms.at(irow) = rms;
    }else if (plane==1){
      threshold = th_factor_ind * rms + 1;
      vplane_rms.at(irow) = rms;
    }else if (plane==2){
      threshold = th_factor_col * rms + 1;
      wplane_rms.at(irow) = rms;
    }
    // std::cout << plane << " " << signal.size() << " " << irow << " " << rms << std::endl;
    
    // create rois
    int roi_begin=-1;
    int roi_end=-1;
    
    std::vector<std::pair<int,int>> temp_rois;
    // now find ROI, above five sigma, and pad with +- six time ticks
    for (int j=0;j<int(signal1.size())-1;j++){
      double content = signal1.at(j);
      double content_tight = signal2.at(j);

      if (content > threshold || 
    	  (content_tight > threshold )){
    	roi_begin = j;
    	roi_end = j;
    	for (int k=j+1;k< int(signal1.size());k++){
    	  if (signal1.at(k) > threshold ||
    	      (signal2.at(k) > threshold)){
    	    roi_end = k;
    	  }else{
    	    break;
    	  }
    	}
    	int temp_roi_begin = roi_begin ; // filter_pad;
    	if (temp_roi_begin <0 ) temp_roi_begin = 0;
    	int temp_roi_end = roi_end ; // filter_pad;
    	if (temp_roi_end >int(signal1.size())-1) temp_roi_end = int(signal1.size())-1;

    	//if (chid == 1151) std::cout << temp_roi_begin << " " << temp_roi_end << std::endl;

	
    	if (temp_rois.size() == 0){
    	  temp_rois.push_back(std::make_pair(temp_roi_begin,temp_roi_end));
    	}else{
    	  if (temp_roi_begin <= temp_rois.back().second){
    	    temp_rois.back().second = temp_roi_end;
    	  }else{
    	    temp_rois.push_back(std::make_pair(temp_roi_begin,temp_roi_end));
    	  }
    	}
    	j = roi_end + 1;
      }
    }

    
    // fill rois ...
    if (plane==0){
      self_rois_u.at(irow) = temp_rois;
    }else if (plane==1){
      self_rois_v.at(irow) = temp_rois;
    }else{
      self_rois_w.at(irow) = temp_rois;
    }
    //    std::cout << plane << " " << irow << " " << temp_rois.size() << std::endl;
  }
  
  extend_ROI_self(plane);
  create_ROI_connect_info(plane);
}

void ROI_formation::find_ROI_by_decon_itself(int plane, const Array::array_xxf& r_data){
  find_ROI_by_decon_itself(plane, r_data,r_data);
}


void ROI_formation::extend_ROI_loose(int plane){

  if (plane==0){
    // compare the loose one with tight one 
    for(int i=0;i!=nwire_u;i++){
      std::vector<std::pair<int,int>> temp_rois;
      for (size_t j=0;j!=loose_rois_u.at(i).size();j++){
	int start = loose_rois_u.at(i).at(j).first;
	int end = loose_rois_u.at(i).at(j).second;
	for (size_t k=0;k!=self_rois_u.at(i).size();k++){
	  int temp_start = self_rois_u.at(i).at(k).first;
	  int temp_end = self_rois_u.at(i).at(k).second;
	  if (start > temp_start && start < temp_end)
	    start = temp_start;
	  // loop through all the tight one to examine start
	  if (end > temp_start && end < temp_end)
	    end = temp_end; 
	  // loop through all the tight one to examine the end
	}
	if (temp_rois.size()==0){
	  temp_rois.push_back(std::make_pair(start,end));
	}else{
	  if (start < temp_rois.back().second){
	    temp_rois.back().second = end;
	  }else{
	    temp_rois.push_back(std::make_pair(start,end));
	  }
	}
      }
      loose_rois_u.at(i) = temp_rois;
    }
  }else if (plane==1){
    for(int i=0;i!=nwire_v;i++){
      std::vector<std::pair<int,int>> temp_rois;
      for (size_t j=0;j!=loose_rois_v.at(i).size();j++){
	int start = loose_rois_v.at(i).at(j).first;
	int end = loose_rois_v.at(i).at(j).second;
	for (size_t k=0;k!=self_rois_v.at(i).size();k++){
	  int temp_start = self_rois_v.at(i).at(k).first;
	  int temp_end = self_rois_v.at(i).at(k).second;
	  if (start > temp_start && start < temp_end)
	    start = temp_start;
	  // loop through all the tight one to examine start
	  if (end > temp_start && end < temp_end)
	    end = temp_end; 
	  // loop through all the tight one to examine the end
	}
	if (temp_rois.size()==0){
	  temp_rois.push_back(std::make_pair(start,end));
	}else{
	  if (start < temp_rois.back().second){
	    temp_rois.back().second = end;
	  }else{
	    temp_rois.push_back(std::make_pair(start,end));
	  }
	}
      }
      loose_rois_v.at(i) = temp_rois;
    }
  }
  
}


double ROI_formation::local_ave(Waveform::realseq_t& signal, int bin, int width){
  double sum1 = 0;
  double sum2 = 0;
  
  for (int i=-width;i<width+1;i++){
    int current_bin = bin + i;

    while (current_bin <0)
      current_bin += signal.size();
    while (current_bin >= int(signal.size()))
      current_bin -= signal.size();
    
    sum1 += signal.at(current_bin);
    sum2 ++;
  }

  if (sum2>0){
    return sum1/sum2;
  }else{
    return 0;
  }
}


int ROI_formation::find_ROI_end(Waveform::realseq_t& signal, int bin, double th ){
  int end = bin;
  double content = signal.at(end);
  while(content>th){
    end ++;
    if (end >=int(signal.size())){
      content = signal.at(end-signal.size());
    }else{
      content = signal.at(end);
    }
    if (end == int(signal.size())) break;
  }

  while(local_ave(signal,end+1,1) < local_ave(signal,end,1)){
    end++;
    if (end == int(signal.size())) break;
  } 
  return end;

}
int ROI_formation::find_ROI_begin(Waveform::realseq_t& signal, int bin, double th ){
  // find the first one before bin and is below threshold ... 
  int begin = bin;
  double content = signal.at(begin);
  while(content > th){
    begin --;
    if (begin <0){
      content = signal.at(begin + int(signal.size()));
    }else{
      content = signal.at(begin);
    }
    if (begin == 0) break;
  }
  
  // calculate the local average
  // keep going and find the minimum
  while( local_ave(signal,begin-1,1) < local_ave(signal,begin,1)){
    begin --;
    if (begin == 0) break;
  }
  
  return begin;
}


void ROI_formation::find_ROI_loose(int plane, const Array::array_xxf& r_data, int rebin){

  
  extend_ROI_loose(plane);
}
