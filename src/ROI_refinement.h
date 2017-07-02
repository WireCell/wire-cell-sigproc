#ifndef WIRECELLSIGPROC_ROIREFINEMENT
#define WIRECELLSIGPROC_ROIREFINEMENT

#include "SignalROI.h"
#include "ROI_formation.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Waveform.h"

#include <vector>
#include <map>

namespace WireCell{
  namespace SigProc{
    class ROI_refinement{
    public:
      ROI_refinement(Waveform::ChannelMaskMap& cmm,int nwire_u, int nwire_v, int nwire_w);
      ~ROI_refinement();

      void Clear();

      // initialize the ROIs
      void load_data(int plane, const Array::array_xxf& r_data, ROI_formation& roi_form);
      
      SignalROIChList& get_u_rois(){return rois_u_loose;};
      SignalROIChList& get_v_rois(){return rois_v_loose;};
      SignalROIChList& get_w_rois(){return rois_w_tight;};
      
    private:
      int nwire_u;
      int nwire_v;
      int nwire_w;
      
      void unlink(SignalROI *prev_roi, SignalROI *next_roi);
      void link(SignalROI *prev_roi, SignalROI *next_roi);

      
      std::map<int,std::vector<std::pair<int,int>>> bad_ch_map;
      
      
      SignalROIChList rois_u_tight;
      SignalROIChList rois_v_tight;
      SignalROIChList rois_w_tight;
    
      SignalROIChList rois_u_loose;
      SignalROIChList rois_v_loose;
   
    
      SignalROIMap front_rois;
      SignalROIMap back_rois;
      SignalROIMap contained_rois;
      
    };
  }
}
#endif
