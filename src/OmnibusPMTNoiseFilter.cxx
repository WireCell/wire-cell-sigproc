#include "WireCellSigProc/OmnibusPMTNoiseFilter.h"

#include "WireCellSigProc/Derivations.h"


#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"


using namespace WireCell;

using namespace WireCell::SigProc;


OmnibusPMTNoiseFilter::OmnibusPMTNoiseFilter(const std::string anode_tn, int pad_window, int min_window_length, int threshold, float rms_threshold )
  : m_anode_tn(anode_tn)
  , m_pad_window(pad_window)
  , m_min_window_length(min_window_length)
  , m_threshold(threshold)
  , m_rms_threshold(rms_threshold)
{
  configure(default_configuration());
}
OmnibusPMTNoiseFilter::~OmnibusPMTNoiseFilter()
{
}

void OmnibusPMTNoiseFilter::configure(const WireCell::Configuration& config)
{
  m_anode_tn = get(config, "anode", m_anode_tn);
  m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
  if (!m_anode) {
    THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
  }

  m_pad_window = get(config,"pad_window", m_pad_window);
  m_min_window_length = get(config,"min_window_length",m_min_window_length);
  m_threshold = get(config,"threshold",m_threshold);
  m_rms_threshold = get(config,"rms_threshold",m_rms_threshold);
}
WireCell::Configuration OmnibusPMTNoiseFilter::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["pad_window"] = m_pad_window;
    cfg["min_window_length"] = m_min_window_length;
    cfg["threshold"] = m_threshold;
    cfg["rms_threshold"]= m_rms_threshold;
    return cfg;
}

bool OmnibusPMTNoiseFilter::operator()(const input_pointer& in, output_pointer& out)
{
  std::map<int, Waveform::realseq_t> bychan_coll;
  std::map<int, Waveform::realseq_t> bychan_indu;
  std::map<int, Waveform::realseq_t> bychan_indv;
  std::map<int,double> by_chan_rms;
  // go through all channels and calculate RMS as well as categorize them
  auto traces = in->traces();
  for (auto trace : *traces.get()) {
    int ch = trace->channel();

    auto wpid = m_anode->resolve(ch);      
    const int iplane = wpid.index();

    Waveform::realseq_t signal=trace->charge();
    std::pair<double,double> results = Derivations::CalcRMS(signal);
    by_chan_rms[ch] = results.second;

    //std::cout << iplane << " " << ch << " " << results.second << std::endl;
    
    if (iplane ==0){
      bychan_indu[ch] = trace->charge();
    }else if (iplane == 1){
      bychan_indv[ch] = trace->charge();
    }else{
      bychan_coll[ch] = trace->charge();
    }
    
  }
  
  // Remove PMT signal from Collection
  for (auto cs : bychan_coll) {
    // ignore dead channels ... 
    if (by_chan_rms[cs.first]>m_rms_threshold)
      RemovePMTSignalCollection(cs.second,by_chan_rms[cs.first],cs.first);
  }
  // ID PMT signal in induction signal

  // Remove PMT signal from Induction ...

  //load results ...

  ITrace::vector itraces;
  for (auto cs : bychan_indu) {
    // fixme: that tbin though
    SimpleTrace *trace = new SimpleTrace(cs.first, 0, cs.second);
    itraces.push_back(ITrace::pointer(trace));
  }
  for (auto cs : bychan_indv) {
    // fixme: that tbin though
    SimpleTrace *trace = new SimpleTrace(cs.first, 0, cs.second);
    itraces.push_back(ITrace::pointer(trace));
  }
  for (auto cs : bychan_coll) {
    // fixme: that tbin though
    SimpleTrace *trace = new SimpleTrace(cs.first, 0, cs.second);
    itraces.push_back(ITrace::pointer(trace));
  }
  
  SimpleFrame* sframe = new SimpleFrame(in->ident(), in->time(), itraces, in->tick(), in->masks());
  out = IFrame::pointer(sframe);
  
  return true;
}


void OmnibusPMTNoiseFilter::RemovePMTSignalCollection(Waveform::realseq_t& signal,double rms, int ch){
  
  int flag_start = 0;
  int start_bin=0;
  int end_bin=0;
  int peak_bin=0;

  for (int i=0;i!=int(signal.size());i++){
    float content = signal.at(i);

    if (flag_start ==0){
      if (content < -m_threshold * rms){
	start_bin = i;
	flag_start = 1;
      }
    }else{
      if (content >= -m_threshold * rms){
	end_bin = i-1;
	if (end_bin > start_bin + m_min_window_length){
	  float min = signal.at(start_bin);
	  peak_bin = start_bin;
	  for (int j=start_bin+1;j!=end_bin;j++){
	    if (signal.at(j) < min)
	      peak_bin = j;
	  }

	  
	  
	}
      }
    }
  }
  
}


void OmnibusPMTNoiseFilter::IDPMTSignalInduction(Waveform::realseq_t& signal, double rms, int ch){
  
}

void OmnibusPMTNoiseFilter::RemovePMTSignalInduction(Waveform::realseq_t& signal, int start_bin, int end_bin){
  
}
