#include "WireCellSigProc/OmnibusSigProc.h"

#include "WireCellUtil/NamedFactory.h"

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

void OmnibusSigProc::load_data(const input_pointer& in){
  u_data = Array::array_xxf::Zero(nwire_u,m_nticks);
  v_data = Array::array_xxf::Zero(nwire_v,m_nticks);
  w_data = Array::array_xxf::Zero(nwire_w,m_nticks);
  
  auto traces = in->traces();
  Array::array_xxf* temp=&u_data;
  int offset = 0;


  //fix me, this mapping needs to be fixed ... 
  for (auto trace : *traces.get()) {
    int tbin = trace->tbin();
    int ch = trace->channel();
    auto charges = trace->charge();

    int plane = ch_plane_map[ch];

    //std::cout << nwire_u << " " << m_nticks << " " << charges.size() << " " << tbin << " " << ch << " " << plane << std::endl;
    
    if (plane == 0){
      temp = &u_data;
      offset = 0;
    }else if (plane==1){
      temp = &v_data;
      offset = nwire_u;
    }else if (plane==2){
      temp = &w_data;
      offset = nwire_u+nwire_v;
    }

    int counter = 0;
    for (auto q : charges) {
      (*temp)(ch - offset, tbin + counter) = q;
      counter ++;
    }
  }

  //ensure dead channels are indeed dead ...
  // zero out the bad channels ...
  for (auto const& it: cmm) {
    if (it.first == "bad"){
      for (auto const &it1 : it.second){
	int chid = it1.first;
	int plane = ch_plane_map[chid];
	if (plane == 0){
	  temp = &u_data;
	  offset = 0;
	}else if (plane==1){
	  temp = &v_data;
	  offset = nwire_u;
	}else if (plane==2){
	  temp = &w_data;
	  offset = nwire_u+nwire_v;
	}

	for (size_t ind = 0; ind < it1.second.size();++ind){
	  for (int i=it1.second[ind].first; i!= it1.second[ind].second; i++){
	    (*temp)(chid-offset,i) = 0;
	  }
	}
	
	
      }
    }

  
  
    //std::cout << u_data(14,2000) << " " << v_data(1000,2000) << " " << w_data(1000,2000) << std::endl;
  
  }

}

void OmnibusSigProc::do_time_fft(){
  uc_data = Array::dft_rc(u_data,0);
  vc_data = Array::dft_rc(v_data,0);
  wc_data = Array::dft_rc(w_data,0);

  // now apply the ch-by-ch response ... 
}

void OmnibusSigProc::do_wire_fft(){
  uc_data = Array::dft_cc(uc_data,1);
  vc_data = Array::dft_cc(vc_data,1);
  wc_data = Array::dft_cc(wc_data,1);

  //std::cout << uc_data1(1000,3000) << std::endl;
}

void OmnibusSigProc::do_wire_inv_fft(){
  uc_data = Array::idft_cc(uc_data,1);
  vc_data = Array::idft_cc(vc_data,1);
  wc_data = Array::idft_cc(wc_data,1);
}

void OmnibusSigProc::do_time_inv_fft(){
  u_data = Array::idft_cr(uc_data,0);
  v_data = Array::idft_cr(vc_data,0);
  w_data = Array::idft_cr(wc_data,0);
}

bool OmnibusSigProc::operator()(const input_pointer& in, output_pointer& out)
{
  cmm = in->masks();

  // load data into EIGEN matrices ...
  load_data(in);

  //std::cout << u_data(1000,3000) << std::endl;

  // Deconvolute
  //do the first time FFT and can correct for the ch-by-ch response
  do_time_fft();

  //std::cout << uc_data(1000,3000) << std::endl;
  //do the second time wire FFT
  do_wire_fft();
  //std::cout << uc_data(1000,3000) << std::endl;
  
  //do the first time inverse FFT in wire
  do_wire_inv_fft();

  //std::cout << uc_data(1000,3000) << std::endl;
  
  //do the second time inverse FFT in time with different filters ...
  do_time_inv_fft();
  
  //std::cout << u_data(1000,3000) << std::endl;
  
  // Form ROIs

  // Refine ROIs

  // Get results
  
  return true;
}
