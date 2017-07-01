
#include "ROI_formation.h"

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_formation::ROI_formation(int nwire_u, int nwire_v, int nwire_w, float th_factor_ind, float th_factor_col, int pad, float asy, int nbins)
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
  , th_factor_ind(th_factor_ind)
  , th_factor_col(th_factor_col)
  , pad(pad)
  , asy(asy)
  , nbins(nbins)
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


void ROI_formation::find_ROI_by_decon_itself(int plane, const Array::array_xxf& r_data){
  
}

void ROI_formation::extend_ROI_self(){
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

void ROI_formation::create_ROI_connect_info(){
  
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



