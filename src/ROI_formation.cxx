
#include "ROI_formation.h"

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_formation::ROI_formation(int nwire_u, int nwire_v, int nwire_w)
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
{
  
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
