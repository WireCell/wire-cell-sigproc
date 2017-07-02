#include "ROI_refinement.h"
#include <iostream>

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_refinement::ROI_refinement(Waveform::ChannelMaskMap& cmm,int nwire_u, int nwire_v, int nwire_w, float th_factor)
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
  , th_factor(th_factor)
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

void ROI_refinement::refine_data(int plane){
  CleanUpROIs(plane);
  generate_merge_ROIs(plane);
}
