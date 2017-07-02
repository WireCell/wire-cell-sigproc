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
      ROI_refinement(Waveform::ChannelMaskMap& cmm,int nwire_u, int nwire_v, int nwire_w, float th_factor = 3.0, float fake_signal_low_th = 1200, float fake_signal_high_th = 1500, int pad = 5, int break_roi_loop = 2, float th_peak = 3.0, float sep_peak = 6.0, float low_peak_sep_threshold_pre = 1200);
      ~ROI_refinement();

      void Clear();

      // initialize the ROIs
      void load_data(int plane, const Array::array_xxf& r_data, ROI_formation& roi_form);
      void refine_data(int plane, ROI_formation& roi_form);
      
      SignalROIChList& get_u_rois(){return rois_u_loose;};
      SignalROIChList& get_v_rois(){return rois_v_loose;};
      SignalROIChList& get_w_rois(){return rois_w_tight;};
      
    private:
      int nwire_u;
      int nwire_v;
      int nwire_w;

      float th_factor;
      float fake_signal_low_th;
      float fake_signal_high_th;
      int pad;
      int break_roi_loop;
      float th_peak;
      float sep_peak;
      float low_peak_sep_threshold_pre;
      
      void unlink(SignalROI *prev_roi, SignalROI *next_roi);
      void link(SignalROI *prev_roi, SignalROI *next_roi);
      void CleanUpROIs(int plane);
      void generate_merge_ROIs(int plane);
      void CheckROIs(int plane, ROI_formation& roi_form);

      void CleanUpCollectionROIs();
      void CleanUpInductionROIs(int plane);
      void ShrinkROIs(int plane, ROI_formation& roi_form);
      void ShrinkROI(SignalROI *roi, ROI_formation& roi_form);

      void BreakROIs(int plane, ROI_formation& roi_form);
      void BreakROI(SignalROI *roi, float rms);
      void BreakROI1(SignalROI *roi);
      
      
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
