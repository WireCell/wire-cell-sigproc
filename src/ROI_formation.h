

#ifndef WIRECELLSIGPROC_ROIFORMATION
#define WIRECELLSIGPROC_ROIFORMATION


#include <vector>


namespace WireCell{
  namespace SigProc{
    class ROI_formation{
    public:
      ROI_formation(int nwire_u, int nwire_v, int nwire_w);
      ~ROI_formation();

      void Clear();
      
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
      int nwire_u, nwire_v, nwire_w;
      
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
