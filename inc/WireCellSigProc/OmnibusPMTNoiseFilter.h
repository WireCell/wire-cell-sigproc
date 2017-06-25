#ifndef WIRECELLSIGPROC_OMNIBUSPMTNOISEFILTER
#define WIRECELLSIGPROC_OMNIBUSPMTNOISEFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/Waveform.h"

namespace WireCell {
  namespace SigProc {
    
    class OmnibusPMTNoiseFilter : public WireCell::IFrameFilter, public WireCell::IConfigurable {
    public:
      typedef std::vector< std::vector<int> > grouped_channels_t;
      
      /// Create an OmnibusPMTNoiseFilter.
      OmnibusPMTNoiseFilter();
      virtual ~OmnibusPMTNoiseFilter();
      
      /// IFrameFilter interface.
      virtual bool operator()(const input_pointer& in, output_pointer& out);
      
      /// IConfigurable interface.
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;
      
      
      /// Explicitly inject required services
      
      
    private:
      
      
    };
  }
}

#endif
