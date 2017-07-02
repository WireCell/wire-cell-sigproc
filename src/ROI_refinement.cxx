#include "ROI_refinement.h"
#include <iostream>
#include <set>

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_refinement::ROI_refinement(Waveform::ChannelMaskMap& cmm,int nwire_u, int nwire_v, int nwire_w, float th_factor, float fake_signal_low_th, float fake_signal_high_th)
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
  , th_factor(th_factor)
  , fake_signal_low_th(fake_signal_low_th)
  , fake_signal_high_th(fake_signal_high_th)
{
  rois_u_tight.resize(nwire_u);
  rois_u_loose.resize(nwire_u);

  rois_v_tight.resize(nwire_v);
  rois_v_loose.resize(nwire_v);

  rois_w_tight.resize(nwire_w);
  
  for (int i=0;i!=nwire_u;i++){
    SignalROIList temp_rois;
    rois_u_tight.at(i) = temp_rois;
  }
  for (int i=0;i!=nwire_v;i++){
    SignalROIList temp_rois;
    rois_v_tight.at(i) = temp_rois;
  }
  for (int i=0;i!=nwire_w;i++){
    SignalROIList temp_rois;
    rois_w_tight.at(i) = temp_rois;
  }
  
  for (int i=0;i!=nwire_u;i++){
    SignalROIList temp_rois;
    rois_u_loose.at(i) = temp_rois;
  }
  for (int i=0;i!=nwire_v;i++){
    SignalROIList temp_rois;
    rois_v_loose.at(i) = temp_rois;
  }

  for (auto it = cmm["bad"].begin(); it!=cmm["bad"].end(); it++){
    int ch = it->first;
    std::vector<std::pair<int,int>> temps;
    bad_ch_map[ch] = temps;
    for (size_t ind = 0; ind < it->second.size(); ind++){
      bad_ch_map[ch].push_back(std::make_pair(it->second.at(ind).first, it->second.at(ind).second));
      //std::cout << ch << " " <<  << std::endl;
    }
  }
  
}

ROI_refinement::~ROI_refinement(){
  Clear();
}

void ROI_refinement::Clear(){
  for (int i=0;i!=nwire_u;i++){
    for (auto it = rois_u_tight.at(i).begin(); it!=rois_u_tight.at(i).end();it++){
      delete *it;
    }
    rois_u_tight.at(i).clear();
    for (auto it = rois_u_loose.at(i).begin(); it!=rois_u_loose.at(i).end();it++){
      delete *it;
    }
    rois_u_loose.at(i).clear();
  }

  for (int i=0;i!=nwire_v;i++){
    for (auto it = rois_v_tight.at(i).begin(); it!=rois_v_tight.at(i).end();it++){
      delete *it;
    }
    rois_v_tight.at(i).clear();
    for (auto it = rois_v_loose.at(i).begin(); it!=rois_v_loose.at(i).end();it++){
      delete *it;
    }
    rois_v_loose.at(i).clear();
  }
  
  for (int i=0;i!=nwire_w;i++){
    for (auto it = rois_w_tight.at(i).begin(); it!=rois_w_tight.at(i).end();it++){
      delete *it;
    }
    rois_w_tight.at(i).clear();
  }
  
  front_rois.clear();
  back_rois.clear();
  contained_rois.clear();
}

void ROI_refinement::unlink(SignalROI* prev_roi, SignalROI* next_roi){
  if (front_rois.find(prev_roi)!=front_rois.end()){
    SignalROISelection& temp_rois = front_rois[prev_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),next_roi);
    if (it != temp_rois.end())
      temp_rois.erase(it);
  }
  if (back_rois.find(next_roi)!=back_rois.end()){
    SignalROISelection& temp_rois = back_rois[next_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),prev_roi);
    if (it != temp_rois.end())
      temp_rois.erase(it);
  }
}

void ROI_refinement::link(SignalROI* prev_roi, SignalROI* next_roi){
  if (front_rois.find(prev_roi)!=front_rois.end()){
    SignalROISelection& temp_rois = front_rois[prev_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),next_roi);
    if (it == temp_rois.end())
      temp_rois.push_back(next_roi);
  }else{
    SignalROISelection temp_rois;
    temp_rois.push_back(next_roi);
    front_rois[prev_roi] = temp_rois;
  }

  if (back_rois.find(next_roi)!=back_rois.end()){
    SignalROISelection& temp_rois = back_rois[next_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),prev_roi);
    if (it == temp_rois.end())
      temp_rois.push_back(prev_roi);
  }else{
    SignalROISelection temp_rois;
    temp_rois.push_back(prev_roi);
    back_rois[next_roi] = temp_rois;
  }
}

void ROI_refinement::load_data(int plane, const Array::array_xxf& r_data, ROI_formation& roi_form){
  // fill RMS 
  std::vector<float> plane_rms;
  int offset = 0;
  if (plane==0){
    offset = 0;
    plane_rms = roi_form.get_uplane_rms();
  }else if (plane==1){
    offset = nwire_u;
    plane_rms = roi_form.get_vplane_rms();
  }else if (plane==2){
    offset = nwire_u+nwire_v;
    plane_rms = roi_form.get_wplane_rms();
  }

  // load data ... 
  for (int irow = 0; irow!=r_data.rows(); irow++){
    Waveform::realseq_t signal(r_data.cols());
    if (bad_ch_map.find(irow+offset)!=bad_ch_map.end()){
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
	  signal.at(icol) = r_data(irow,icol);
	}else{
	  signal.at(icol) = 0;
	}
      }
    }else{
      for (int icol = 0; icol!= r_data.cols(); icol++){
	signal.at(icol) = r_data(irow,icol);
      }
    }

    int chid = irow+offset;
    // load tight rois
    std::vector<std::pair<int,int>>& uboone_rois = roi_form.get_self_rois(irow+offset);
    for (size_t i=0;i!=uboone_rois.size();i++){
      SignalROI *tight_roi = new SignalROI(plane,irow+offset, uboone_rois.at(i).first,uboone_rois.at(i).second, signal);
      float threshold = plane_rms.at(irow) * th_factor;
      if (tight_roi->get_above_threshold(threshold).size()==0) {
	delete tight_roi;
	continue;
      }
      
      if (plane==0){
	rois_u_tight[chid].push_back(tight_roi);
       	if (chid>0){
       	  //form connectivity map
       	  for (auto it = rois_u_tight[chid-1].begin();it!=rois_u_tight[chid-1].end();it++){
       	    SignalROI *prev_roi = *it;
      	    if (tight_roi->overlap(prev_roi)){
      	      if (front_rois.find(prev_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		front_rois[prev_roi] = temp_rois;
      	      }else{
      		front_rois[prev_roi].push_back(tight_roi);
      	      }
      	      if (back_rois.find(tight_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(prev_roi);
      		back_rois[tight_roi] = temp_rois;
      	      }else{
      		back_rois[tight_roi].push_back(prev_roi);
      	      }
      	    }
       	  }
     	}
	if (chid < nwire_u-1){
	  // add the code for the next one to be completed
	  for (auto it = rois_u_tight[chid+1].begin();it!=rois_u_tight[chid+1].end();it++){
	    SignalROI *next_roi = *it;
	    if (tight_roi->overlap(next_roi)){
	      if (back_rois.find(next_roi) == back_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(tight_roi);
		back_rois[next_roi] = temp_rois;
	      }else{
		back_rois[next_roi].push_back(tight_roi);
	      }
	      
	      if (front_rois.find(tight_roi) == front_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(next_roi);
		front_rois[tight_roi] = temp_rois;
	      }else{
		front_rois[tight_roi].push_back(next_roi);
	      }
	    }
	  }
	}
      }else if (plane==1){
	rois_v_tight[chid-nwire_u].push_back(tight_roi);
      	if (chid>nwire_u){
      	  //form connectivity map
      	  for (auto it = rois_v_tight[chid - nwire_u-1].begin();it!=rois_v_tight[chid-nwire_u-1].end();it++){
      	    SignalROI *prev_roi = *it;
      	    if (tight_roi->overlap(prev_roi)){
      	      if (front_rois.find(prev_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		front_rois[prev_roi] = temp_rois;
      	      }else{
      		front_rois[prev_roi].push_back(tight_roi);
      	      }
      	      if (back_rois.find(tight_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(prev_roi);
      		back_rois[tight_roi] = temp_rois;
      	      }else{
      		back_rois[tight_roi].push_back(prev_roi);
      	      }
      	    }
      	  }
      	}
	
	if (chid<nwire_u+nwire_v-1){
      	  //form connectivity map
      	  for (auto it = rois_v_tight[chid - nwire_u+1].begin();it!=rois_v_tight[chid-nwire_u+1].end();it++){
      	    SignalROI *next_roi = *it;
      	    if (tight_roi->overlap(next_roi)){
      	      if (back_rois.find(next_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		back_rois[next_roi] = temp_rois;
      	      }else{
      		back_rois[next_roi].push_back(tight_roi);
      	      }
      	      if (front_rois.find(tight_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(next_roi);
      		front_rois[tight_roi] = temp_rois;
      	      }else{
      		front_rois[tight_roi].push_back(next_roi);
      	      }
      	    }
      	  }
      	}
      }else {
	rois_w_tight[chid-nwire_u-nwire_v].push_back(tight_roi);
	
      	if (chid>nwire_u+nwire_v){
      	  //form connectivity map
      	  for (auto it = rois_w_tight[chid-nwire_u-nwire_v-1].begin();it!=rois_w_tight[chid-nwire_u-nwire_v-1].end();it++){
      	    SignalROI *prev_roi = *it;
      	    if (tight_roi->overlap(prev_roi)){
      	      if (front_rois.find(prev_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		front_rois[prev_roi] = temp_rois;
      	      }else{
      		front_rois[prev_roi].push_back(tight_roi);
      	      }
      	      if (back_rois.find(tight_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(prev_roi);
      		back_rois[tight_roi] = temp_rois;
      	      }else{
      		back_rois[tight_roi].push_back(prev_roi);
      	      }
      	    }
      	  }
      	}


	if (chid<nwire_u+nwire_v+nwire_w-1){
      	  //form connectivity map
      	  for (auto it = rois_w_tight[chid-nwire_u-nwire_v+1].begin();it!=rois_w_tight[chid-nwire_u-nwire_v+1].end();it++){
      	    SignalROI *next_roi = *it;
      	    if (tight_roi->overlap(next_roi)){
      	      if (back_rois.find(next_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		back_rois[next_roi] = temp_rois;
      	      }else{
      		back_rois[next_roi].push_back(tight_roi);
      	      }
      	      if (front_rois.find(tight_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(next_roi);
      		front_rois[tight_roi] = temp_rois;
      	      }else{
      		front_rois[tight_roi].push_back(next_roi);
      	      }
      	    }
      	  }
      	}
      }

      
      
    }// loop over tight rois ... 

    if (plane!=2){
      uboone_rois = roi_form.get_loose_rois(chid);
      for (size_t i = 0; i!=uboone_rois.size();i++){
	SignalROI *loose_roi = new SignalROI(plane,chid,uboone_rois.at(i).first,uboone_rois.at(i).second,signal);
	float threshold = plane_rms.at(irow) * th_factor;
	if (loose_roi->get_above_threshold(threshold).size()==0) {
	  delete loose_roi;
	  continue;
	}
	if (plane==0){
	  rois_u_loose[chid].push_back(loose_roi);
	  
	  if (chid>0){
	    //form connectivity map
	    for (auto it=rois_u_loose[chid-1].begin();it!=rois_u_loose[chid-1].end();it++){
	      SignalROI *prev_roi = *it;
	      if (loose_roi->overlap(prev_roi)){
		if (front_rois.find(prev_roi) == front_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(loose_roi);
		  front_rois[prev_roi] = temp_rois;
		}else{
		  front_rois[prev_roi].push_back(loose_roi);
		}
		if (back_rois.find(loose_roi) == back_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(prev_roi);
		  back_rois[loose_roi] = temp_rois;
		}else{
		  back_rois[loose_roi].push_back(prev_roi);
		}
	      }
	    }
	  }
	  
	  if (chid<nwire_u-1){
	    //form connectivity map
	    for (auto it=rois_u_loose[chid+1].begin();it!=rois_u_loose[chid+1].end();it++){
	      SignalROI *next_roi = *it;
	      if (loose_roi->overlap(next_roi)){
		if (back_rois.find(next_roi) == back_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(loose_roi);
		  back_rois[next_roi] = temp_rois;
		}else{
		  back_rois[next_roi].push_back(loose_roi);
		}
		if (front_rois.find(loose_roi) == front_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(next_roi);
		  front_rois[loose_roi] = temp_rois;
		}else{
		  front_rois[loose_roi].push_back(next_roi);
		}
	      }
	    }
	  }
	  
	  //form contained map ... 
	  for (auto it=rois_u_tight[chid].begin();it!=rois_u_tight[chid].end();it++){
	    SignalROI *tight_roi = *it;
	    if (tight_roi->overlap(loose_roi)){
	      if (contained_rois.find(loose_roi)==contained_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(tight_roi);
		contained_rois[loose_roi] = temp_rois;
	      }else{
		contained_rois[loose_roi].push_back(tight_roi);
	      }
	    }
	  }
	}else if (plane==1){
	  rois_v_loose[chid-nwire_u].push_back(loose_roi);
	  
	  if (chid>nwire_u){
	    //form connectivity map
	    for (auto it = rois_v_loose[chid-nwire_u-1].begin();it!=rois_v_loose[chid-nwire_u-1].end();it++){
	      SignalROI *prev_roi = *it;
	      if (loose_roi->overlap(prev_roi)){
		if (front_rois.find(prev_roi) == front_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(loose_roi);
		  front_rois[prev_roi] = temp_rois;
		}else{
		  front_rois[prev_roi].push_back(loose_roi);
		}
		if (back_rois.find(loose_roi) == back_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(prev_roi);
		  back_rois[loose_roi] = temp_rois;
		}else{
		  back_rois[loose_roi].push_back(prev_roi);
		}
	      }
	    }
	  }
	  
	  if (chid<nwire_u+nwire_v-1){
	    //form connectivity map
	    for (auto it = rois_v_loose[chid-nwire_u+1].begin();it!=rois_v_loose[chid-nwire_u+1].end();it++){
	      SignalROI *next_roi = *it;
	      if (loose_roi->overlap(next_roi)){
		if (back_rois.find(next_roi) == back_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(loose_roi);
		  back_rois[next_roi] = temp_rois;
		}else{
		  back_rois[next_roi].push_back(loose_roi);
		}
		if (front_rois.find(loose_roi) == front_rois.end()){
		  SignalROISelection temp_rois;
		  temp_rois.push_back(next_roi);
		  front_rois[loose_roi] = temp_rois;
		}else{
		  front_rois[loose_roi].push_back(next_roi);
		}
	      }
	    }
	  }
	  
	  //form contained map ... 
	  for (auto it = rois_v_tight[chid-nwire_u].begin();it!=rois_v_tight[chid-nwire_u].end();it++){
	    SignalROI *tight_roi = *it;
	    if (tight_roi->overlap(loose_roi)){
	      if (contained_rois.find(loose_roi)==contained_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(tight_roi);
		contained_rois[loose_roi] = temp_rois;
	      }else{
		contained_rois[loose_roi].push_back(tight_roi);
	      }
	    }
	  }
	}
      }
    }
    
  } // loop over signal rows
  
}

void ROI_refinement::CleanUpROIs(int plane){

  // clean up ROIs
  std::map<SignalROI*, int> ROIsaved_map;

  if (plane==0){
    //int counter = 0;
    for (size_t i=0;i!=rois_u_loose.size();i++){
      // counter += rois_u_loose.at(i).size();
      for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
	SignalROI *roi = *it;
	if (ROIsaved_map.find(roi)==ROIsaved_map.end()){
	  if (contained_rois.find(roi) != contained_rois.end()){
	    // contain good stuff
	    SignalROISelection temp_rois;
	    temp_rois.push_back(roi);
	    ROIsaved_map[roi] = 1;
	    
	    while(temp_rois.size()){
	      SignalROI *temp_roi = temp_rois.back();
	      temp_rois.pop_back();
	      // save all its neighbour into a temporary holder
	      if (front_rois.find(temp_roi)!=front_rois.end()){
		for (auto it1 = front_rois[temp_roi].begin();it1!=front_rois[temp_roi].end();it1++){
		  if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		    temp_rois.push_back(*it1);
		    ROIsaved_map[*it1] = 1;
		  }
		}
	      }
	      if (back_rois.find(temp_roi)!=back_rois.end()){
		for (auto it1 = back_rois[temp_roi].begin();it1!=back_rois[temp_roi].end();it1++){
		  if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		    temp_rois.push_back(*it1);
		    ROIsaved_map[*it1] = 1;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    //remove the bad ones ...
    //int counter2 = 0;
    for (size_t i=0;i!=rois_u_loose.size();i++){
      SignalROISelection to_be_removed;
      for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
	SignalROI *roi = *it;
	if (ROIsaved_map.find(roi) == ROIsaved_map.end()){
	  //	counter2 ++;
	  to_be_removed.push_back(roi);
	  //it = rois_u_loose.at(i).erase(it);
	  // check contained map
	  if (contained_rois.find(roi)!= contained_rois.end()){
	    std::cout << "Wrong! " << std::endl;
	  }
	  // check front map
	  if (front_rois.find(roi)!=front_rois.end()){
	    for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++){
	      auto it2 = find(back_rois[*it1].begin(),back_rois[*it1].end(),roi);
	      back_rois[*it1].erase(it2);
	    }
	    front_rois.erase(roi);
	  }
	  // check back map
	  if (back_rois.find(roi)!=back_rois.end()){
	    for (auto it1 = back_rois[roi].begin(); it1!=back_rois[roi].end(); it1++){
	      auto it2 = find(front_rois[*it1].begin(),front_rois[*it1].end(),roi);
	      front_rois[*it1].erase(it2);
	    }
	    back_rois.erase(roi);
	  }
	}
      }
      
      for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
	auto it1 = find(rois_u_loose.at(i).begin(), rois_u_loose.at(i).end(),*it);
	rois_u_loose.at(i).erase(it1);
      }
    }
  }else if (plane==1){
    
    // int counter = 0;
    for (size_t i=0;i!=rois_v_loose.size();i++){
      //  counter += rois_v_loose.at(i).size();
      for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end();it++){
	SignalROI *roi = *it;
	if (ROIsaved_map.find(roi)==ROIsaved_map.end()){
	  if (contained_rois.find(roi) != contained_rois.end()){
	    // contain good stuff
	    SignalROISelection temp_rois;
	    temp_rois.push_back(roi);
	    ROIsaved_map[roi] = 1;
	    
	    while(temp_rois.size()){
	      SignalROI *temp_roi = temp_rois.back();
	      temp_rois.pop_back();
	      // save all its neighbour into a temporary holder
	      if (front_rois.find(temp_roi)!=front_rois.end()){
		for (auto it1 = front_rois[temp_roi].begin();it1!=front_rois[temp_roi].end();it1++){
		  if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		    temp_rois.push_back(*it1);
		    ROIsaved_map[*it1] = 1;
		  }
		}
	      }
	      if (back_rois.find(temp_roi)!=back_rois.end()){
		for (auto it1 = back_rois[temp_roi].begin();it1!=back_rois[temp_roi].end();it1++){
		  if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		    temp_rois.push_back(*it1);
		    ROIsaved_map[*it1] = 1;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    //remove the bad ones ...
    //int counter2 = 0;
    for (size_t i=0;i!=rois_v_loose.size();i++){
      SignalROISelection to_be_removed;
      for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end();it++){
	SignalROI *roi = *it;
	if (ROIsaved_map.find(roi) == ROIsaved_map.end()){
	  //	counter2 ++;
	  to_be_removed.push_back(roi);
	  //it = rois_v_loose.at(i).erase(it);
	  // check contained map
	  if (contained_rois.find(roi)!= contained_rois.end()){
	    std::cout << "Wrong! " << std::endl;
	  }
	  // check front map
	  if (front_rois.find(roi)!=front_rois.end()){
	    for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++){
	      auto it2 = find(back_rois[*it1].begin(),back_rois[*it1].end(),roi);
	      back_rois[*it1].erase(it2);
	    }
	    front_rois.erase(roi);
	  }
	  // check back map
	  if (back_rois.find(roi)!=back_rois.end()){
	    for (auto it1 = back_rois[roi].begin(); it1!=back_rois[roi].end(); it1++){
	      auto it2 = find(front_rois[*it1].begin(),front_rois[*it1].end(),roi);
	      front_rois[*it1].erase(it2);
	    }
	    back_rois.erase(roi);
	  }
	}
      }
      
      for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
	auto it1 = find(rois_v_loose.at(i).begin(), rois_v_loose.at(i).end(),*it);
	rois_v_loose.at(i).erase(it1);
      }
    }
  }



  // int counter1 = 0;
  // for (int i=0;i!=rois_v_loose.size();i++){
  //   counter1+=rois_v_loose.at(i).size();
  // }
  
  // std::cout << counter << " " << ROIsaved_map.size() << " " << counter1 << " " << counter2 << std::endl;
}

void ROI_refinement::generate_merge_ROIs(int plane){
  // find tight ROIs not contained by the loose ROIs
  if (plane==0){
    for (int i = 0;i!=nwire_u;i++){
      std::map<SignalROI*,int> covered_tight_rois;
      for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
	SignalROI *roi = *it;
	if (contained_rois.find(roi) != contained_rois.end()){
	  for (auto it1 = contained_rois[roi].begin(); it1!= contained_rois[roi].end(); it1++){
	    if (covered_tight_rois.find(*it1)==covered_tight_rois.end()){
	      covered_tight_rois[*it1]  =1;
	    }
	  }
	}
      }
      SignalROISelection saved_rois;
      for (auto it = rois_u_tight.at(i).begin();it!=rois_u_tight.at(i).end();it++){
	SignalROI *roi = *it;
	if (covered_tight_rois.find(roi) == covered_tight_rois.end()){
	  saved_rois.push_back(roi);
	}
      }
      // if (i == 1212)
      //   std::cout << saved_rois.size() << " " << saved_rois.at(0)->get_start_bin() << " " << saved_rois.at(0)->get_end_bin() << std::endl;
      
      for (auto it = saved_rois.begin(); it!=saved_rois.end();it++){
	SignalROI *roi = *it;
	// Duplicate them 
	SignalROI *loose_roi = new SignalROI(roi);
	
	rois_u_loose.at(i).push_back(loose_roi);
	
	// update all the maps     
	// contained
	SignalROISelection temp_rois;
	temp_rois.push_back(roi);
	contained_rois[loose_roi] = temp_rois;
	// front map loose ROI
	if (i < nwire_u-1){
	  for (auto it1 = rois_u_loose.at(i+1).begin(); it1!=rois_u_loose.at(i+1).end(); it1++){
	    SignalROI *next_roi = *it1;
	    
	    if (loose_roi->overlap(next_roi)){
	      if (back_rois.find(next_roi) == back_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(loose_roi);
		back_rois[next_roi] = temp_rois;
	      }else{
		back_rois[next_roi].push_back(loose_roi);
	      }
	      
	      if (front_rois.find(loose_roi) == front_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(next_roi);
		front_rois[loose_roi] = temp_rois;
	      }else{
		front_rois[loose_roi].push_back(next_roi);
	      }
	    }
	    
	  }
	}
	// back map loose ROI
	if (i > 0){
	  for (auto it1 = rois_u_loose.at(i-1).begin(); it1!=rois_u_loose.at(i-1).end(); it1++){
	    SignalROI *prev_roi = *it1;
	    if (loose_roi->overlap(prev_roi)){
	      if (front_rois.find(prev_roi) == front_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(loose_roi);
		front_rois[prev_roi] = temp_rois;
	      }else{
		front_rois[prev_roi].push_back(loose_roi);
	      }
	      if (back_rois.find(loose_roi) == back_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(prev_roi);
		back_rois[loose_roi] = temp_rois;
	      }else{
		back_rois[loose_roi].push_back(prev_roi);
	      }
	    }
	  }
	}
      }
    }
  }else if (plane==1){
    for (int i = 0;i!=nwire_v;i++){
      std::map<SignalROI*,int> covered_tight_rois;
      for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
	SignalROI *roi = *it;
	if (contained_rois.find(roi) != contained_rois.end()){
	  for (auto it1 = contained_rois[roi].begin(); it1!= contained_rois[roi].end(); it1++){
	    if (covered_tight_rois.find(*it1)==covered_tight_rois.end()){
	      covered_tight_rois[*it1]  =1;
	    }
	  }
	}
      }
      SignalROISelection saved_rois;
      for (auto it = rois_v_tight.at(i).begin();it!=rois_v_tight.at(i).end();it++){
	SignalROI *roi = *it;
	if (covered_tight_rois.find(roi) == covered_tight_rois.end()){
	  saved_rois.push_back(roi);
	}
      }
      //if (i == 3885-2400)
      //  std::cout << saved_rois.size() << std::endl;
      //   std::cout << saved_rois.size() << " " << saved_rois.at(0)->get_start_bin() << " " << saved_rois.at(0)->get_end_bin() << std::endl;
      
      for (auto it = saved_rois.begin(); it!=saved_rois.end();it++){
	SignalROI *roi = *it;
	// Duplicate them 
	SignalROI *loose_roi = new SignalROI(roi);
	
	rois_v_loose.at(i).push_back(loose_roi);
	
	// update all the maps     
	// contained
	SignalROISelection temp_rois;
	temp_rois.push_back(roi);
	contained_rois[loose_roi] = temp_rois;
	// front map loose ROI
	if (i < nwire_v-1){
	  for (auto it1 = rois_v_loose.at(i+1).begin(); it1!=rois_v_loose.at(i+1).end(); it1++){
	    SignalROI *next_roi = *it1;
	    
	    if (loose_roi->overlap(next_roi)){
	      if (back_rois.find(next_roi) == back_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(loose_roi);
		back_rois[next_roi] = temp_rois;
	      }else{
		back_rois[next_roi].push_back(loose_roi);
	      }
	      
	      if (front_rois.find(loose_roi) == front_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(next_roi);
		front_rois[loose_roi] = temp_rois;
	      }else{
		front_rois[loose_roi].push_back(next_roi);
	      }
	    }
	    
	  }
	}
	// back map loose ROI
	if (i > 0){
	  for (auto it1 = rois_v_loose.at(i-1).begin(); it1!=rois_v_loose.at(i-1).end(); it1++){
	    SignalROI *prev_roi = *it1;
	    if (loose_roi->overlap(prev_roi)){
	      if (front_rois.find(prev_roi) == front_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(loose_roi);
		front_rois[prev_roi] = temp_rois;
	      }else{
		front_rois[prev_roi].push_back(loose_roi);
	      }
	      if (back_rois.find(loose_roi) == back_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(prev_roi);
		back_rois[loose_roi] = temp_rois;
	      }else{
		back_rois[loose_roi].push_back(prev_roi);
	      }
	    }
	  }
	}
      }
    }
  }
}

void ROI_refinement::CheckROIs(int plane,ROI_formation& roi_form){

  if (plane==0){
    std::vector<float>& rms_u = roi_form.get_uplane_rms();
    
    // for (int i=0;i!=rois_u_loose.size();i++){
    //   for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
    //     SignalROI *roi = *it;
    //     int chid = roi->get_chid();
    //     if (chid != i) 
    // 	std::cout << roi << std::endl;
    
    //     if (front_rois.find(roi)!=front_rois.end()){
    //  	for (auto it1 = front_rois[roi].begin();it1!=front_rois[roi].end();it1++){
    // 	  SignalROI *roi1 = *it1;
    // 	  int chid1 = roi1->get_chid();
    // 	  if (chid1!=i+1)
    // 	    std::cout << roi1 << " " << chid << " " << chid1 << std::endl;
    // 	}
    //     }
    //   }
    // }
    
    //loop over u loose
    for (size_t i=0;i!=rois_u_loose.size();i++){
      for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
	SignalROI *roi = *it;
	int chid = roi->get_chid();
	float th;
	th = th_factor*rms_u.at(chid);
	
	if (front_rois.find(roi)!=front_rois.end()){
	  SignalROISelection temp_rois;
	  for (auto it1 = front_rois[roi].begin();it1!=front_rois[roi].end();it1++){
	    SignalROI *roi1 = *it1;
	    int chid1 = roi1->get_chid();
	    //std::cout << "F " << i << " " << rois_u_loose.size() << " " << roi1 << " " << chid << " " << chid1 << std::endl;
	    float th1;
	    th1 = th_factor*rms_u.at(chid1);
	    if (roi->overlap(roi1,th,th1)){
	    }else{
	      temp_rois.push_back(roi1);
	      //unlink(roi,roi1);
	    }
	  }
	  for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	    unlink(roi,*it2);
	  }
	}

	if (back_rois.find(roi)!=back_rois.end()){
	  SignalROISelection temp_rois;
	  for (auto it1 = back_rois[roi].begin();it1!=back_rois[roi].end();it1++){
	    SignalROI *roi1 = *it1;
	    int chid1 = roi1->get_chid();
	    //std::cout << "B " << roi1 << " " << chid << " " << chid1 << std::endl;
	    float th1;
	    th1 = th_factor*rms_u.at(chid1);
	    if (roi->overlap(roi1,th,th1)){
	    }else{
	      temp_rois.push_back(roi1);
	      //unlink(roi,roi1);
	    }
	  }
	  for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	    unlink(roi,*it2);
	  }
	}
	
      }
    }
  }else if (plane==1){
    std::vector<float>& rms_v = roi_form.get_vplane_rms();
    //loop over v loose
    for (size_t i=0;i!=rois_v_loose.size();i++){
      for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end();it++){
	SignalROI *roi = *it;
	int chid = roi->get_chid()-nwire_u;
	float th;
	th = th_factor*rms_v.at(chid);
	if (front_rois.find(roi)!=front_rois.end()){
	  SignalROISelection temp_rois;
	  for (auto it1 = front_rois[roi].begin();it1!=front_rois[roi].end();it1++){
	    SignalROI *roi1 = *it1;
	    int chid1 = roi1->get_chid()-nwire_u;
	    float th1;
	    th1 = th_factor*rms_v.at(chid1);
	    if (roi->overlap(roi1,th,th1)){
	    }else{
	      temp_rois.push_back(roi1);
	      // unlink(roi,roi1);
	    }
	  }
	  for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	    unlink(roi,*it2);
	  }
	}
	
	if (back_rois.find(roi)!=back_rois.end()){
	  SignalROISelection temp_rois;
	  for (auto it1 = back_rois[roi].begin();it1!=back_rois[roi].end();it1++){
	    SignalROI *roi1 = *it1;
	    int chid1 = roi1->get_chid()-nwire_u;
	    float th1;
	    th1 = th_factor*rms_v.at(chid1);
	    if (roi->overlap(roi1,th,th1)){
	    }else{
	      temp_rois.push_back(roi1);
	      // unlink(roi,roi1);
	    }
	  }
	  for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	    unlink(roi,*it2);
	  }
	}
      }
    }
  }
}

void ROI_refinement::CleanUpCollectionROIs(){
  // deal with tight ROIs, 
  // scan with all the tight ROIs to look for peaks above certain threshold, put in a temporary set
  float threshold = fake_signal_low_th; //electrons, about 1/2 of MIP per tick ...
  std::set<SignalROI*> Good_ROIs;
  for (int i=0;i!=nwire_w;i++){
    for (auto it = rois_w_tight.at(i).begin();it!=rois_w_tight.at(i).end();it++){
      SignalROI* roi = *it;
      if (roi->get_above_threshold(threshold).size()!=0)
	Good_ROIs.insert(roi);
    }
  }
  // for a particular ROI if it is not in, or it is not connected with one in the temporary map, then remove it
  std::list<SignalROI*> Bad_ROIs;
  for (int i=0;i!=nwire_w;i++){
    for (auto it = rois_w_tight.at(i).begin();it!=rois_w_tight.at(i).end();it++){
      SignalROI* roi = *it;
      
      if (Good_ROIs.find(roi)!=Good_ROIs.end()) continue;
      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection next_rois = front_rois[roi];
	int flag_qx = 0;
	for (size_t i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection next_rois = back_rois[roi];
	int flag_qx = 0;
	for (size_t i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      Bad_ROIs.push_back(roi);
    }
  }
  
  // remove the ROI and then update the map
  
  for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    SignalROI* roi = *it;
    int chid = roi->get_chid()-nwire_u-nwire_v;
    //std::cout << chid << std::endl;
    if (front_rois.find(roi)!=front_rois.end()){
      SignalROISelection next_rois = front_rois[roi];
      for (size_t i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      front_rois.erase(roi);
    }
    
    if (back_rois.find(roi)!=back_rois.end()){
      SignalROISelection next_rois = back_rois[roi];
      for (size_t i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      back_rois.erase(roi);
    }
    auto it1 = find(rois_w_tight.at(chid).begin(), rois_w_tight.at(chid).end(),roi);
    if (it1 != rois_w_tight.at(chid).end())
      rois_w_tight.at(chid).erase(it1);
    
    delete roi;
  }
  
}

void ROI_refinement::CleanUpInductionROIs(int plane){
   // deal with loose ROIs
  // focus on the isolated ones first
  float threshold = fake_signal_high_th;
  std::list<SignalROI*> Bad_ROIs;
  if (plane==0){
    for (int i=0;i!=nwire_u;i++){
      for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
	SignalROI* roi = *it;
	if (front_rois.find(roi)==front_rois.end() && back_rois.find(roi)==back_rois.end()){
	  if (roi->get_above_threshold(threshold).size()==0)
	    Bad_ROIs.push_back(roi);
	}
      }
    }
    for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
      SignalROI* roi = *it;
      int chid = roi->get_chid();
      auto it1 = find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(),roi);
      if (it1 != rois_u_loose.at(chid).end())
	rois_u_loose.at(chid).erase(it1);
      delete roi;
    }
    Bad_ROIs.clear();
  }else if (plane==1){
    for (int i=0;i!=nwire_v;i++){
      for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
	SignalROI* roi = *it;
	if (front_rois.find(roi)==front_rois.end() && back_rois.find(roi)==back_rois.end()){
	  if (roi->get_above_threshold(threshold).size()==0)
	    Bad_ROIs.push_back(roi);
	}
      }
    }
    for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
      SignalROI* roi = *it;
      int chid = roi->get_chid()-nwire_u;
      auto it1 = find(rois_v_loose.at(chid).begin(), rois_v_loose.at(chid).end(),roi);
      if (it1 != rois_v_loose.at(chid).end())
	rois_v_loose.at(chid).erase(it1);
      delete roi;
    }
  }


  threshold = fake_signal_low_th;
  std::set<SignalROI*> Good_ROIs;
  if (plane==0){
    for (int i=0;i!=nwire_u;i++){
      for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
	SignalROI* roi = *it;
	if (roi->get_above_threshold(threshold).size()!=0)
	  Good_ROIs.insert(roi);
      }
    }
    Bad_ROIs.clear();
    for (int i=0;i!=nwire_u;i++){
      for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
	SignalROI* roi = *it;
	
	if (Good_ROIs.find(roi)!=Good_ROIs.end()) continue;
	if (front_rois.find(roi)!=front_rois.end()){
	  SignalROISelection next_rois = front_rois[roi];
	  int flag_qx = 0;
	  for (size_t i=0;i!=next_rois.size();i++){
	    SignalROI* roi1 = next_rois.at(i);
	    if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	      flag_qx = 1;
	      continue;
	    }
	  }
	  if (flag_qx == 1) continue;
	}
	
	if (back_rois.find(roi)!=back_rois.end()){
	  SignalROISelection next_rois = back_rois[roi];
	  int flag_qx = 0;
	  for (size_t i=0;i!=next_rois.size();i++){
	    SignalROI* roi1 = next_rois.at(i);
	    if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	      flag_qx = 1;
	      continue;
	    }
	  }
	  if (flag_qx == 1) continue;
	}
	
	Bad_ROIs.push_back(roi);
      }
    }
    for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
      SignalROI* roi = *it;
      int chid = roi->get_chid();
      //std::cout << chid << std::endl;
      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection next_rois = front_rois[roi];
	for (size_t i=0;i!=next_rois.size();i++){
	  //unlink the current roi
	  unlink(roi,next_rois.at(i));
	}
	front_rois.erase(roi);
      }
      
      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection next_rois = back_rois[roi];
	for (size_t i=0;i!=next_rois.size();i++){
	  //unlink the current roi
	  unlink(roi,next_rois.at(i));
	}
	back_rois.erase(roi);
      }
      auto it1 = find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(),roi);
      if (it1 != rois_u_loose.at(chid).end())
	rois_u_loose.at(chid).erase(it1);
      
      delete roi;
    }
  }else if (plane==1){
    
    Good_ROIs.clear();
    for (int i=0;i!=nwire_v;i++){
      for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
	SignalROI* roi = *it;
	if (roi->get_above_threshold(threshold).size()!=0)
	  Good_ROIs.insert(roi);
      }
    }
    Bad_ROIs.clear();
    for (int i=0;i!=nwire_v;i++){
      for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
	SignalROI* roi = *it;
	
	if (Good_ROIs.find(roi)!=Good_ROIs.end()) continue;
	if (front_rois.find(roi)!=front_rois.end()){
	  SignalROISelection next_rois = front_rois[roi];
	  int flag_qx = 0;
	  for (size_t i=0;i!=next_rois.size();i++){
	    SignalROI* roi1 = next_rois.at(i);
	    if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	      flag_qx = 1;
	      continue;
	    }
	  }
	  if (flag_qx == 1) continue;
	}
	
	if (back_rois.find(roi)!=back_rois.end()){
	  SignalROISelection next_rois = back_rois[roi];
	  int flag_qx = 0;
	  for (size_t i=0;i!=next_rois.size();i++){
	    SignalROI* roi1 = next_rois.at(i);
	    if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	      flag_qx = 1;
	      continue;
	    }
	  }
	  if (flag_qx == 1) continue;
	}
	
	Bad_ROIs.push_back(roi);
      }
    }
    for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
      SignalROI* roi = *it;
      int chid = roi->get_chid()-nwire_u;
      //std::cout << chid << std::endl;
      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection next_rois = front_rois[roi];
	for (size_t i=0;i!=next_rois.size();i++){
	  //unlink the current roi
	  unlink(roi,next_rois.at(i));
	}
	front_rois.erase(roi);
      }
      
      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection next_rois = back_rois[roi];
	for (size_t i=0;i!=next_rois.size();i++){
	  //unlink the current roi
	  unlink(roi,next_rois.at(i));
	}
	back_rois.erase(roi);
      }
      auto it1 = find(rois_v_loose.at(chid).begin(), rois_v_loose.at(chid).end(),roi);
      if (it1 != rois_v_loose.at(chid).end())
	rois_v_loose.at(chid).erase(it1);
      
      delete roi;
    }
  }
}

void ROI_refinement::ShrinkROI(SignalROI *roi, ROI_formation& roi_form){
  // get tight ROI as a inner boundary
  // get the nearby ROIs with threshold as some sort of boundary 
  int start_bin = roi->get_start_bin();
  int end_bin = roi->get_end_bin();
  if (start_bin <0 || end_bin <0) return;
  
  int chid = roi->get_chid();
  int plane = roi->get_plane();
  std::vector<float>& contents = roi->get_contents();
  
  float threshold1=0;
  if (plane==0){
    threshold1 = roi_form.get_uplane_rms().at(chid) * 3.0;
  }else if (plane==1){
    threshold1 = roi_form.get_vplane_rms().at(chid-nwire_u)*3.0;
  }
  
  int channel_save = 1240;
  int print_flag = 0;

  // std::cout << "check tight ROIs " << std::endl;
  // use to save contents
  Waveform::realseq_t temp_signal(end_bin-start_bin+1,0);
  // TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);
  
  // check tight ROIs
  if (contained_rois.find(roi)!=contained_rois.end()){
    for (auto it = contained_rois[roi].begin();it!=contained_rois[roi].end();it++){
      SignalROI *tight_roi = *it;
      int start_bin1 = tight_roi->get_start_bin();
      int end_bin1 = tight_roi->get_end_bin();
      
      if (chid == channel_save && print_flag)
   	std::cout << "Tight "  " " << start_bin1 << " " << end_bin1 << std::endl;

      for (int i=start_bin1;i<=end_bin1;i++){
   	if (i-start_bin >=0 && i-start_bin <int(temp_signal.size())){
	  // 	  htemp->SetBinContent(i-start_bin+1,1);
	  temp_signal.at(i-start_bin) = 1;
   	}
      }
    }
  }

  // std::cout << "check front ROIs " << std::endl;

  //check front ROIs
  if (front_rois.find(roi)!=front_rois.end()){
    for (auto it=front_rois[roi].begin();it!=front_rois[roi].end();it++){
      SignalROI *next_roi = *it;
      int start_bin1 = next_roi->get_start_bin();
      int chid1 = next_roi->get_chid();
      int plane1 = next_roi->get_plane();
      float threshold=0;
      if (plane1==0){
   	threshold = roi_form.get_uplane_rms().at(chid1) * th_factor;
      }else if (plane1==1){
   	threshold = roi_form.get_vplane_rms().at(chid1-nwire_u) * th_factor;
      }
      std::vector<std::pair<int,int>> contents_above_threshold = next_roi->get_above_threshold(threshold);
      for (size_t i=0;i!=contents_above_threshold.size();i++){
   	if (chid == channel_save && print_flag)
   	  std::cout << "Front " << chid1 << " " << start_bin1 + contents_above_threshold.at(i).first << " " << start_bin1 + contents_above_threshold.at(i).second << std::endl;

   	for (int j=contents_above_threshold.at(i).first;j<=contents_above_threshold.at(i).second;j++){
   	  if (j+start_bin1-start_bin >=0 && j+start_bin1-start_bin < int(temp_signal.size())){
   	    if (contents.at(j+start_bin1-start_bin) > threshold1)
	      temp_signal.at(j+start_bin1-start_bin) = 1;
	    //	      htemp->SetBinContent(j+start_bin1-start_bin+1,1);
   	  }
   	}
      }
    }
  }
  
  //std::cout << "check back ROIs " << std::endl;

  //check back ROIs
  if (back_rois.find(roi)!=back_rois.end()){
    for (auto it=back_rois[roi].begin();it!=back_rois[roi].end();it++){
      SignalROI *prev_roi = *it;
      int start_bin1 = prev_roi->get_start_bin();
      int chid1 = prev_roi->get_chid();
      int plane1 = prev_roi->get_plane();
      float threshold = 0;
      if (plane1==0){
   	threshold = roi_form.get_uplane_rms().at(chid1) * th_factor;
      }else if (plane1==1){
   	threshold = roi_form.get_vplane_rms().at(chid1-nwire_u)* th_factor;
      }
      std::vector<std::pair<int,int>> contents_above_threshold = prev_roi->get_above_threshold(threshold);
      for (size_t i=0;i!=contents_above_threshold.size();i++){
	if (chid == channel_save && print_flag)
	  std::cout << "Back " << chid1 << " " << start_bin1 + contents_above_threshold.at(i).first << " " << start_bin1 + contents_above_threshold.at(i).second << std::endl;

   	for (int j=contents_above_threshold.at(i).first;j<=contents_above_threshold.at(i).second;j++){
   	  if (j+start_bin1-start_bin >=0 && j+start_bin1-start_bin <int(temp_signal.size())){
   	    if (contents.at(j+start_bin1-start_bin) > threshold1)
	      temp_signal.at(j+start_bin1-start_bin) = 1;
	    //	      htemp->SetBinContent(j+start_bin1-start_bin+1,1);
   	  }
   	}
      }
    }
  }

  //std::cout << "check contents " << std::endl;

  // // consider the 1/2 of the peak as threshold;
  // float max = 0;
  // for (int i=0;i!=contents.size();i++){
  //   if (contents.at(i) > max)
  //     max = contents.at(i);
  // }
  // for (int i=0;i!=contents.size();i++){
  //   // if (contents.at(i) > max/2. && contents.at(i) > threshold1*2 ) htemp->SetBinContent(i+1,1);
  // }
  
  // get the first bin, and last bin, add pad
  int pad = 5;
  int new_start_bin=start_bin;
  int new_end_bin=end_bin;
  for (size_t i=0;i!= temp_signal.size(); i++){
    if (temp_signal.at(i) >0){
      new_start_bin = i + start_bin;
      break;
    }
  }
  for (int i = int(temp_signal.size())-1;i>=0;i--){
    if (temp_signal.at(i) > 0){
      new_end_bin = i + start_bin;
      break;
    }
  }
  new_start_bin -= pad;
  new_end_bin += pad;
  if (new_start_bin < start_bin) new_start_bin = start_bin;
  if (new_end_bin > end_bin) new_end_bin = end_bin;
  
  if (chid == channel_save && print_flag)
    std::cout << "check contents " << " " << start_bin << " " << end_bin << " " << new_start_bin << " " << new_end_bin << std::endl;
  

  // create a new ROI
  Waveform::realseq_t signal(end_bin+1);
  // TH1F *h1 = new TH1F("h1","h1",end_bin+1,0,end_bin+1);
  for (int i=new_start_bin; i<=new_end_bin;i++){
    signal.at(i) = contents.at(i-start_bin);
    //   h1->SetBinContent(i+1,contents.at(i-start_bin));
  }
  
  SignalROISelection new_rois;
  if (new_start_bin >=0 && new_end_bin > new_start_bin){
    SignalROI *new_roi = new SignalROI(plane,chid,new_start_bin,new_end_bin,signal);
    new_rois.push_back(new_roi);
  }

  // std::cout << "update maps " << std::endl;
  
  // update the list 
  if (chid < nwire_u){
    auto it = std::find(rois_u_loose.at(chid).begin(),rois_u_loose.at(chid).end(),roi);
    rois_u_loose.at(chid).erase(it);
    for (size_t i=0;i!=new_rois.size();i++){
      rois_u_loose.at(chid).push_back(new_rois.at(i));
    }
  }else if (chid < nwire_u+nwire_v){
    auto it = std::find(rois_v_loose.at(chid-nwire_u).begin(),rois_v_loose.at(chid-nwire_u).end(),roi);
    rois_v_loose.at(chid-nwire_u).erase(it);
    for (size_t i=0;i!=new_rois.size();i++){
      rois_v_loose.at(chid-nwire_u).push_back(new_rois.at(i));
    }
  }
  
  // update all the maps 
  // update front map
  if (front_rois.find(roi)!=front_rois.end()){
    SignalROISelection next_rois = front_rois[roi];
    for (size_t i=0;i!=next_rois.size();i++){
      //unlink the current roi
      unlink(roi,next_rois.at(i));
      //loop new rois and link them
      for (size_t j=0;j!=new_rois.size();j++){
	if (new_rois.at(j)->overlap(next_rois.at(i)))
	  link(new_rois.at(j),next_rois.at(i));
      }
    }
    front_rois.erase(roi);
  }
  // update back map
  if (back_rois.find(roi)!=back_rois.end()){
    SignalROISelection prev_rois = back_rois[roi];
    for (size_t i=0;i!=prev_rois.size();i++){
      // unlink the current roi
      unlink(prev_rois.at(i),roi);
      // loop new rois and link them
      for (size_t j=0;j!=new_rois.size();j++){
	if (new_rois.at(j)->overlap(prev_rois.at(i)))
	  link(prev_rois.at(i),new_rois.at(j));
      }
    }
    back_rois.erase(roi);
  }
  
  // update contained map 
  if (contained_rois.find(roi)!=contained_rois.end()){
    SignalROISelection tight_rois = contained_rois[roi];
    for (size_t i=0;i!=tight_rois.size();i++){
      for (size_t j=0;j!=new_rois.size();j++){
	if (new_rois.at(j)->overlap(tight_rois.at(i))){
	  if (contained_rois.find(new_rois.at(j)) == contained_rois.end()){
	    SignalROISelection temp_rois;
	    temp_rois.push_back(tight_rois.at(i));
	    contained_rois[new_rois.at(j)] = temp_rois;
	  }else{
	    contained_rois[new_rois.at(j)].push_back(tight_rois.at(i));
	  }
	}
      }
    }
    contained_rois.erase(roi);
  }
  
  // delete the old ROI
  delete roi;

  // delete htemp;
  // delete h1;
}

void ROI_refinement::ShrinkROIs(int plane, ROI_formation& roi_form){
  // collect all ROIs
  SignalROISelection all_rois;
  if (plane==0){
    for (size_t i=0;i!=rois_u_loose.size();i++){
      for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end(); it++){
	all_rois.push_back(*it);
      }
    }
  }else if (plane==1){
    for (size_t i=0;i!=rois_v_loose.size();i++){
      for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end(); it++){
	all_rois.push_back(*it);
      }
    }
  }
  for (size_t i=0;i!=all_rois.size();i++){
    ShrinkROI(all_rois.at(i),roi_form);
  }
}

void ROI_refinement::BreakROI(SignalROI *roi, float rms){

}
void ROI_refinement::BreakROI1(SignalROI *roi){
  
}
void ROI_refinement::BreakROIs(int plane, ROI_formation& roi_form){
  SignalROISelection all_rois;
  if (plane==0){
    std::vector<float>& rms_u = roi_form.get_uplane_rms();
    for (size_t i=0;i!=rois_u_loose.size();i++){
      for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end(); it++){
	BreakROI(*it,rms_u.at(i));
	all_rois.push_back(*it);
	
      }
    }
  }else if (plane==1){
    std::vector<float>& rms_v = roi_form.get_vplane_rms();
    for (size_t i=0;i!=rois_v_loose.size();i++){
      for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end(); it++){
	BreakROI(*it,rms_v.at(i));
	all_rois.push_back(*it);
      }
    }
  }
  
  for (size_t i=0;i!=all_rois.size();i++){
    // if (all_rois.at(i)->get_chid()==1151){
    //   std::cout << all_rois.at(i)->get_chid() << " " << all_rois.at(i)->get_start_bin() << " " << all_rois.at(i)->get_end_bin() << std::endl;
    //   for (int j=0;j!=all_rois.at(i)->get_contents().size();j++){
    // 	std::cout << j << " " << all_rois.at(i)->get_contents().at(j) << std::endl;  
    //   }
    // }
    
    BreakROI1(all_rois.at(i));
  }

  
}


void ROI_refinement::refine_data(int plane, ROI_formation& roi_form){
  CleanUpROIs(plane);
  generate_merge_ROIs(plane);
  
}
