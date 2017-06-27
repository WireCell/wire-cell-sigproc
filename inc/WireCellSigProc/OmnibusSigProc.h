#ifndef WIRECELLSIGPROC_OMNIBUSSIGPROC
#define WIRECELLSIGPROC_OMNIBUSSIGPROC

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"

namespace WireCell {
  namespace SigProc {
    class OmnibusSigProc : public WireCell::IFrameFilter, public WireCell::IConfigurable {
    public:
      OmnibusSigProc(const std::string anode_tn = "AnodePlane");
      virtual ~OmnibusSigProc();
      
      virtual bool operator()(const input_pointer& in, output_pointer& out);
      
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;
      
    private:
      // Anode plane for geometry
      std::string m_anode_tn;
      IAnodePlane::pointer m_anode;
      
      // Overall time offset

      
      
      // various software filters ...
      
      // one HF for ROI,  two LF for ROI finding

      // one HF for charge

      // one HF for hit

      // two HF filters for wire dimension 
      
    };
  }
}


#endif
