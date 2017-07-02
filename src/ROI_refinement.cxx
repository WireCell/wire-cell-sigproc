#include "ROI_refinement.h"
#include <iostream>

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_refinement::ROI_refinement(int nwire_u, int nwire_v, int nwire_w)
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
