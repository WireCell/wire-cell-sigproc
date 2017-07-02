

#ifndef WIRECELLSIGPROC_ROIFORMATION
#define WIRECELLSIGPROC_ROIFORMATION

#include "WireCellUtil/Array.h"
#include "WireCellUtil/Waveform.h"

#include <vector>
#include <map>

namespace WireCell{
  namespace SigProc{
    class ROI_formation{
    public:
      ROI_formation(Waveform::ChannelMaskMap& cmm,int nwire_u, int nwire_v, int nwire_w, float th_factor_ind = 3, float th_factor_col = 5, int pad = 5, float asy = 0.1, int nbins = 9594);
      ~ROI_formation();

      void Clear();

      void find_ROI_by_decon_itself(int plane, const Array::array_xxf& r_data);
      void find_ROI_by_decon_itself(int plane, const Array::array_xxf& r_data, const Array::array_xxf& r_data_tight);
      void extend_ROI_self(int plane);
      void create_ROI_connect_info(int plane);
      
      std::vector<std::pair<int,int>>& get_self_rois(int chid) {
	if (chid < nwire_u){
	  return self_rois_u.at(chid);
	}else if (chid < nwire_u + nwire_v){
	  return self_rois_v.at(chid - nwire_u);
	}else{
	  return self_rois_w.at(chid - nwire_u - nwire_v);
	}
      }
      
      std::vector<std::pair<int,int>>& get_loose_rois(int chid) {
	if (chid < nwire_u){
	  return loose_rois_u.at(chid);
	}else if (chid < nwire_u + nwire_v){
	  return loose_rois_v.at(chid - nwire_u);
	}else{
	  return loose_rois_w.at(chid - nwire_u - nwire_v);
	}
      }
      
      std::vector <float>& get_uplane_rms(){return uplane_rms;};
      std::vector <float>& get_vplane_rms(){return vplane_rms;};
      std::vector <float>& get_wplane_rms(){return wplane_rms;};
      
      
    private:
      double cal_RMS(Waveform::realseq_t signal);
      
      int nwire_u, nwire_v, nwire_w;
      float th_factor_ind, th_factor_col;
      int pad;
      float asy;
      int nbins;

      std::map<int,std::vector<std::pair<int,int>>> bad_ch_map;
      
      std::vector<std::vector<std::pair<int,int>>> self_rois_u; // tight ROIs
      std::vector<std::vector<std::pair<int,int>>> self_rois_v; // tight ROIs
      std::vector<std::vector<std::pair<int,int>>> self_rois_w; // tight ROIs
      
      std::vector<std::vector<std::pair<int,int>>> loose_rois_u; // tight ROIs
      std::vector<std::vector<std::pair<int,int>>> loose_rois_v; // tight ROIs
      std::vector<std::vector<std::pair<int,int>>> loose_rois_w; // tight ROIs
      
      std::vector<float> uplane_rms; // calibrated field response
      std::vector<float> vplane_rms; // calibrated field response
      std::vector<float> wplane_rms; // calibrated field response
    };
  }
}

#endif
