#include "ROI_refinement.h"
#include <iostream>

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_refinement::ROI_refinement(Waveform::ChannelMaskMap& cmm,int nwire_u, int nwire_v, int nwire_w)
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
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
  
}
