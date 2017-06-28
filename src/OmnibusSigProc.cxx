#include "WireCellSigProc/OmnibusSigProc.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusSigProc::OmnibusSigProc(const std::string anode_tn, double fine_time_offset, double coarse_time_offset, double period, int nticks)
  : m_anode_tn (anode_tn)
  , m_fine_time_offset(fine_time_offset)
  , m_coarse_time_offset(coarse_time_offset)
  , m_period(period)
  , m_nticks(nticks)
{
  configure(default_configuration());
  // get wires for each plane

  nwire_u = 0;
  nwire_v = 0;
  nwire_w = 0;
  for (int i=0;i!=int(m_anode->channels().size());i++){
    int ch = m_anode->channels().at(i);
    auto wpid = m_anode->resolve(ch);      
    int iplane = wpid.index();
    ch_plane_map[ch] = iplane;

    if (iplane==0){
      nwire_u ++;
    }else if (iplane==1){
      nwire_v ++;
    }else if (iplane==2){
      nwire_w ++;
    }
    
  }
  //std::cout << m_anode->channels().size() << " " << nwire_u << " " << nwire_v << " " << nwire_w << std::endl;
  
}

OmnibusSigProc::~OmnibusSigProc()
{
}

void OmnibusSigProc::configure(const WireCell::Configuration& config)
{
  m_fine_time_offset = get(config,"ftoffset",m_fine_time_offset);
  m_coarse_time_offset = get(config,"ctoffset",m_coarse_time_offset);
  m_anode_tn = get(config, "anode", m_anode_tn);
  m_nticks = get(config,"nticks",m_nticks);
  m_period = get(config,"period",m_period);
  
  m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
  if (!m_anode) {
    THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
  }
}

WireCell::Configuration OmnibusSigProc::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = m_anode_tn;
  cfg["ftoffset"] = m_fine_time_offset;
  cfg["ctoffset"] = m_coarse_time_offset;
  cfg["nticks"] = m_nticks;
  cfg["period"] = m_period;
  
  return cfg;
  
}

void OmnibusSigProc::load_data(const input_pointer& in, int plane){
  if (plane ==0){
    r_data = Array::array_xxf::Zero(nwire_u,m_nticks);
  }else if (plane==1){
    r_data = Array::array_xxf::Zero(nwire_v,m_nticks);
  }else if (plane==2){
    r_data = Array::array_xxf::Zero(nwire_w,m_nticks);
  }
  
  auto traces = in->traces();
  int offset=0;
  if (plane==0){
    offset = 0;
  }else if (plane==1){
    offset = nwire_u;
  }else if (plane==2){
    offset = nwire_u+nwire_v;
  }


  //fix me, this mapping needs to be fixed ... 
  for (auto trace : *traces.get()) {
    int tbin = trace->tbin();
    int ch = trace->channel();
    auto charges = trace->charge();


    if (plane==ch_plane_map[ch]){
      int counter = 0;
      for (auto q : charges) {
	r_data(ch - offset, tbin + counter) = q;
	counter ++;
      }

      //ensure dead channels are indeed dead ...
      if (cmm["bad"].find(ch)!=cmm["bad"].end()){
	for (size_t ind = 0; ind < cmm["bad"][ch].size();++ind){
 	  for (int i=cmm["bad"][ch][ind].first; i!=cmm["bad"][ch][ind].second; i++){
 	    r_data(ch-offset,i) = 0;
 	  }
 	}
      }
      
    }
  }
  
  // std::cout << r_data(14,2000) << " " << r_data(1000,2000) << " " << r_data(1000,3000) << std::endl;
  
}


void OmnibusSigProc::save_data(ITrace::vector& itraces, int plane, int total_offset){

  int offset=0;
  int nwire;
  if (plane==0){
    offset = total_offset;
    nwire = nwire_u;
  }else if (plane==1){
    offset = nwire_u + total_offset;
    nwire = nwire_v;
  }else if (plane==2){
    offset = nwire_u + nwire_v + total_offset;
    nwire = nwire_w;
  }
  
  
  for (int i=0;i!=nwire;i++){
    ITrace::ChargeSequence charge;
    for (int j=0;j!=m_nticks;j++){
      charge.push_back(r_data(i,j));
    }
    SimpleTrace *trace = new SimpleTrace(i+offset, 0, charge);
    itraces.push_back(ITrace::pointer(trace));
  }
  
}



void OmnibusSigProc::decon_2D(int plane){

  // first round of FFT on time
  c_data = Array::dft_rc(r_data,0);

  // now apply the ch-by-ch response ...
  // to be added
  
  //second round of FFT on wire
  c_data = Array::dft_cc(c_data,1);

  //make ratio to the response and apply wire filter
  // to be added
  
  //do the first round of inverse FFT on wire
  c_data = Array::idft_cc(c_data,1);
  
  //do the second round of inverse FFT on wire
  c_data = Array::idft_cr(c_data,0);
}




bool OmnibusSigProc::operator()(const input_pointer& in, output_pointer& out)
{
  cmm = in->masks();
  ITrace::vector itraces;

  for (int i=0;i!=3;i++){
    // load data into EIGEN matrices ...
    load_data(in,i);
    // initial decon ... 
    decon_2D(i);
    // Form ROIs
    
    // Refine ROIs

    // merge results ...
    
    // Get results
    save_data(itraces,i);
  }
  
  SimpleFrame* sframe = new SimpleFrame(in->ident(), in->time(), itraces, in->tick(), cmm);
  out = IFrame::pointer(sframe);
  
  return true;
}
