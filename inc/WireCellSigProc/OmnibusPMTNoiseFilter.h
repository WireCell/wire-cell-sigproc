#ifndef WIRECELLSIGPROC_OMNIBUSPMTNOISEFILTER
#define WIRECELLSIGPROC_OMNIBUSPMTNOISEFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Waveform.h"

namespace WireCell {
  namespace SigProc {
    
    class OmnibusPMTNoiseFilter : public WireCell::IFrameFilter, public WireCell::IConfigurable {
    public:
      typedef std::vector< std::vector<int> > grouped_channels_t;
      
      /// Create an OmnibusPMTNoiseFilter.
      OmnibusPMTNoiseFilter(const std::string anode_tn = "AnodePlane", int pad_window = 5, int min_window_length = 4, int threshold = 5, float rms_threshold = 0.5);
      virtual ~OmnibusPMTNoiseFilter();
      
      /// IFrameFilter interface.
      virtual bool operator()(const input_pointer& in, output_pointer& out);
      
      /// IConfigurable interface.
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;
      
      
      /// Explicitly inject required services
      void RemovePMTSignalCollection(Waveform::realseq_t& signal,double rms, int ch);
      void IDPMTSignalInduction(Waveform::realseq_t& signal, double rms, int ch);
      void RemovePMTSignalInduction(Waveform::realseq_t& signal, int start_bin, int end_bin);
      
    private:
      std::string m_anode_tn;
      IAnodePlane::pointer m_anode;

      int m_pad_window;
      int m_min_window_length;
      int m_threshold;
      float m_rms_threshold;
      
    };
  }
}

#endif
