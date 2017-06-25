#include "WireCellSigProc/OmnibusPMTNoiseFilter.h"

#include "WireCellSigProc/Derivations.h"


#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "PMTNoiseROI.h"

using namespace WireCell;

using namespace WireCell::SigProc;


std::vector<PMTNoiseROI*> PMT_ROIs;

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
  for (auto cs : bychan_indu) {
    // ignore dead channels ... 
    if (by_chan_rms[cs.first]>m_rms_threshold)
      IDPMTSignalInduction(cs.second,by_chan_rms[cs.first],cs.first,0);
  }
  for (auto cs : bychan_indv) {
    // ignore dead channels ... 
    if (by_chan_rms[cs.first]>m_rms_threshold)
      IDPMTSignalInduction(cs.second,by_chan_rms[cs.first],cs.first,1);
  }
  
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

	  PMTNoiseROI *ROI = new PMTNoiseROI(start_bin,end_bin,peak_bin,ch,signal.at(peak_bin));

	  //std::cout << start_bin << " " << end_bin << std::endl;
	  if (PMT_ROIs.size()==0){
	    PMT_ROIs.push_back(ROI);
	  }else{
	    bool flag_merge = false;
	    for (int i=0;i!=int(PMT_ROIs.size());i++){
	      flag_merge = PMT_ROIs.at(i)->merge_ROI(*ROI);
	      if (flag_merge){
		delete ROI;
		break;
	      }
	    }
	    if (!flag_merge){
	      PMT_ROIs.push_back(ROI);
	    }
	  }
	  // std::cout << "h " << PMT_ROIs.size() << std::endl;
	  
	  flag_start = 0;


	   //adaptive baseline
	    float start_content =signal.at(start_bin);
	    for (int j=start_bin;j>=start_bin - m_pad_window;j--){
	      if (j<0) continue;
	      if (fabs(signal.at(j)) < fabs(start_content)){
		start_bin = j;
		start_content = signal.at(start_bin);
	      }
	    }
	    float end_content = signal.at(end_bin);
	    for (int j=end_bin; j<=end_bin+m_pad_window;j++){
	      if (j>=int(signal.size())) continue;
	      if (fabs(signal.at(j)) < fabs(end_content)){
		end_bin = j;
		end_content = signal.at(end_bin);
	      }
	    }
	    
	    for (int j=start_bin;j<=end_bin;j++){
	      float content = start_content + (end_content - start_content) * (j - start_bin) / (end_bin - start_bin*1.0);
	      signal.at(j)=content;
	    }
	  
	  
	}
      }
    }
  }
  
}


void OmnibusPMTNoiseFilter::IDPMTSignalInduction(Waveform::realseq_t& signal, double rms, int ch, int plane){
  //std::cout <<  PMT_ROIs.size() << std::endl;

  for (int i=0;i!=int(PMT_ROIs.size());i++){
    PMTNoiseROI* ROI= PMT_ROIs.at(i);
    for (int j=0;j!= int(ROI->get_peaks().size());j++){
      int peak = ROI->get_peaks().at(j);
      int peak_m1 = peak - 1; if (peak_m1 <0) peak_m1 = 0;
      int peak_m2 = peak - 2; if (peak_m2 <0) peak_m2 = 0;
      int peak_m3 = peak - 3; if (peak_m3 <0) peak_m3 = 0;
      int peak_p1 = peak + 1; if (peak_p1 >= int(signal.size())) peak_p1 = int(signal.size())-1;
      int peak_p2 = peak + 2; if (peak_p2 >= int(signal.size())) peak_p2 = int(signal.size())-1;
      int peak_p3 = peak + 3; if (peak_p3 >= int(signal.size())) peak_p3 = int(signal.size())-1;
      
      if (fabs(signal.at(peak))> m_threshold * rms && 
	  fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) > fabs(signal.at(peak_m1)) + fabs(signal.at(peak_m2)) + fabs(signal.at(peak)) &&
	  fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) > fabs(signal.at(peak_p1)) + fabs(signal.at(peak_p2)) + fabs(signal.at(peak)) && 
	  fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) > fabs(signal.at(peak_m2)) + fabs(signal.at(peak_m3)) + fabs(signal.at(peak_m1)) && 
	  fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) > fabs(signal.at(peak_p2)) + fabs(signal.at(peak_p3)) + fabs(signal.at(peak_p1)) ){
	
	//	    fabs(signal.at(peak))>= fabs(signal.at(peak_m2)) && 
	// fabs(signal.at(peak))>= fabs(signal.at(peak_p1)) && 
	// fabs(signal.at(peak))>= fabs(signal.at(peak_p2)) && 
	// fabs(signal.at(peak))>= fabs(signal.at(peak_p3)) && 
	// fabs(signal.at(peak))>= fabs(signal.at(peak_m3)) 
	// ){
	if (plane == 0 ){
	  ROI->insert_uwires(ch,fabs(signal.at(peak)));
	  break;
	}else{
	  ROI->insert_vwires(ch,fabs(signal.at(peak)));
	  break;
	}
      }
    }
  }
}

void OmnibusPMTNoiseFilter::RemovePMTSignalInduction(Waveform::realseq_t& signal, int start_bin, int end_bin){
  
}
