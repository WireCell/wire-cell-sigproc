#ifndef WIRECELLSIGPROC_OMNIBUSSIGPROC
#define WIRECELLSIGPROC_OMNIBUSSIGPROC

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"

namespace WireCell {
  namespace SigProc {
    class OmnibusSigProc : public WireCell::IFrameFilter, public WireCell::IConfigurable {
    public:
      OmnibusSigProc(const std::string anode_tn = "AnodePlane", double fine_time_offset = 0, double coarse_time_offset = 0);
      virtual ~OmnibusSigProc();
      
      virtual bool operator()(const input_pointer& in, output_pointer& out);
      
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;
      
    private:
      
      // Anode plane for geometry
      std::string m_anode_tn;
      IAnodePlane::pointer m_anode;
      
      // Overall time offset
      double m_fine_time_offset;
      double m_coarse_time_offset;
      
      
    
      
    };
  }
}


#endif
