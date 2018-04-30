#ifndef WIRECELLSIGPROC_OMNIBUSSIGPROC
#define WIRECELLSIGPROC_OMNIBUSSIGPROC

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Array.h"

namespace WireCell {
  namespace SigProc {
    class OmnibusSigProc : public WireCell::IFrameFilter, public WireCell::IConfigurable {
    public:
      OmnibusSigProc(const std::string anode_tn = "AnodePlane",
                     double fine_time_offset = 0.0 * units::microsecond,
                     double coarse_time_offset = -8.0 * units::microsecond,
                     double gain = 14.0 * units::mV/units::fC,
                     double shaping_time = 2.0 * units::microsecond,
                     double inter_gain = 1.2,
                     double ADC_mV = 4096/(2000.*units::mV),
                     bool flag_ch_corr = true,
                     float th_factor_ind = 3,
                     float th_factor_col = 5,
                     int pad = 5,
                     float asy = 0.1,
                     int rebin =6,
                     double l_factor=3.5,
                     double l_max_th=10000,
                     double l_factor1=0.7,
                     int l_short_length = 3,
                     double r_th_factor = 3.0,
                     double r_fake_signal_low_th = 500,
                     double r_fake_signal_high_th = 1000,
                     int r_pad = 5,
                     int r_break_roi_loop = 2,
                     double r_th_peak = 3.0,
                     double r_sep_peak=6.0,
                     double r_low_peak_sep_threshold_pre = 1200,
                     int r_max_npeaks = 200,
                     double r_sigma = 2.0,
                     double r_th_percent = 0.1,
                     int charge_ch_offset = 10000 );
      virtual ~OmnibusSigProc();
      
      virtual bool operator()(const input_pointer& in, output_pointer& out);
      
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;
      
    private:
      
      // convert data into Eigen Matrix
      void load_data(const input_pointer& in, int plane);

      // deconvolution
      void decon_2D_init(int plane); // main decon code 
      void decon_2D_ROI_refine(int plane);
      void decon_2D_tightROI(int plane);
      void decon_2D_tighterROI(int plane); 
      void decon_2D_looseROI(int plane);
      void decon_2D_hits(int plane);
      void decon_2D_charge(int plane);
      
      // save data into the out frame and collect the indices
      void save_data(ITrace::vector& itraces, IFrame::trace_list_t& indices, int plane);

      // initialize the overall response function ...
      void init_overall_response(IFrame::pointer frame);

      void restore_baseline(WireCell::Array::array_xxf& arr);
      
      
      // Anode plane for geometry
      std::string m_anode_tn;
      IAnodePlane::pointer m_anode;
      
      // Overall time offset
      double m_fine_time_offset; // must be positive, between 0-0.5 us, shift the response function to earlier time --> shift the deconvoluted signal to a later time
      double m_coarse_time_offset; // additional coarse time shift ...
      double m_intrinsic_time_offset;
      int m_wire_shift[3];
      
      // bins
      double m_period;
      int m_nticks;
      
      // gain, shaping time, other applification factors
      double m_gain, m_shaping_time;
      double m_inter_gain, m_ADC_mV;
      bool m_flag_ch_corr;

      // some parameters for ROI creating
      float m_th_factor_ind;
      float m_th_factor_col;
      int m_pad;
      float m_asy ;
      int m_rebin;
      double m_l_factor;
      double m_l_max_th;
      double m_l_factor1;
      int m_l_short_length;


       // ROI_refinement
      double m_r_th_factor;
      double m_r_fake_signal_low_th;
      double m_r_fake_signal_high_th;
      int m_r_pad;
      int m_r_break_roi_loop;
      double m_r_th_peak;
      double m_r_sep_peak;
      double m_r_low_peak_sep_threshold_pre;
      int m_r_max_npeaks;
      double m_r_sigma;
      double m_r_th_percent;

      // channel offset
      int m_charge_ch_offset;
      
      // Some global data useful
      int nwire_u, nwire_v, nwire_w;
      Waveform::ChannelMaskMap cmm;
      std::map<int,int> ch_plane_map;

      // data after decon steps before final ifft ...
      Array::array_xxf m_r_data; // evil
      Array::array_xxc m_c_data; // evil
      
      //average overall responses
      std::vector<Waveform::realseq_t> overall_resp[3];
      
    };
  }
}


#endif
// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
