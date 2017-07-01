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
      OmnibusSigProc(const std::string anode_tn = "AnodePlane", double fine_time_offset = 0.2 * units::microsecond, double coarse_time_offset = -8.5 * units::microsecond, double period = 0.5*units::microsecond, int nticks = 9594, double gain = 14.0 * units::mV/units::fC, double shaping_time = 2.0 * units::microsecond, double inter_gain = 1.2, double ADC_mV = 4096/2000., bool flag_ch_corr = false);
      virtual ~OmnibusSigProc();
      
      virtual bool operator()(const input_pointer& in, output_pointer& out);
      
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;
      
    private:
      
      // convert data into Eigen Matrix
      void load_data(const input_pointer& in, int plane);

      // deconvolution
      void decon_2D_init(int plane); // main decon code 
      void decon_2D_tightROI(int plane);
      void decon_2D_tighterROI(int plane); 
      void decon_2D_looseROI(int plane);
      void decon_2D_hits(int plane);
      void decon_2D_charge(int plane);
      
      // save data into the out frame
      void save_data(ITrace::vector& itraces, int plane, int total_offset=0);

      // initialize the overall response function ...
      void init_overall_response();

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
      
      // Some global data useful
      int nwire_u, nwire_v, nwire_w;
      Waveform::ChannelMaskMap cmm;
      std::map<int,int> ch_plane_map;

      // data after decon steps before final ifft ...
      Array::array_xxf r_data;
      Array::array_xxc c_data;
      
      //average overall responses
      std::vector<Waveform::realseq_t> overall_resp[3];
      
    };
  }
}


#endif
