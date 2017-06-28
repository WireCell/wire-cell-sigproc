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
      OmnibusSigProc(const std::string anode_tn = "AnodePlane", double fine_time_offset = 0, double coarse_time_offset = 0, double period = 0.5*units::microsecond, int nticks = 9594);
      virtual ~OmnibusSigProc();
      
      virtual bool operator()(const input_pointer& in, output_pointer& out);
      
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;
      
    private:
      
      // convert data into Eigen Matrix
      void load_data(const input_pointer& in, int plane);

      // deconvolution
      void decon_2D(int plane);
      
      // save data into the out frame
      void save_data(ITrace::vector& itraces, int plane, int total_offset=0);
      

      
      
      // Anode plane for geometry
      std::string m_anode_tn;
      IAnodePlane::pointer m_anode;
      
      // Overall time offset
      double m_fine_time_offset;
      double m_coarse_time_offset;

      // bins
      double m_period;
      int m_nticks;

      
      // Some global data useful
      int nwire_u, nwire_v, nwire_w;
      Waveform::ChannelMaskMap cmm;
      std::map<int,int> ch_plane_map;

      // data after decon steps before final ifft ...
      Array::array_xxf r_data;
      Array::array_xxc c_data;
      
    };
  }
}


#endif
