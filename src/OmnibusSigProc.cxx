#include "WireCellSigProc/OmnibusSigProc.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellIface/IFieldResponse.h"

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusSigProc::OmnibusSigProc(const std::string anode_tn, double fine_time_offset, double coarse_time_offset, double period, int nticks, double gain , double shaping_time , double inter_gain , double ADC_mV )
  : m_anode_tn (anode_tn)
  , m_fine_time_offset(fine_time_offset)
  , m_coarse_time_offset(coarse_time_offset)
  , m_period(period)
  , m_nticks(nticks)
  , m_gain(gain)
  , m_shaping_time(shaping_time)
  , m_inter_gain(inter_gain)
  , m_ADC_mV(ADC_mV)
{
  configure(default_configuration());
  // get wires for each plane

 
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

  m_gain = get(config,"gain",m_gain);
  m_shaping_time = get(config,"shaping_time",m_shaping_time);
  m_inter_gain = get(config,"inter_gain", m_inter_gain);
  m_ADC_mV = get(config,"ADC_mV", m_ADC_mV);

  
  m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
  if (!m_anode) {
    THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
  }


  // give the configuration, trying to understand how many wires are there for each plane
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
  
}

WireCell::Configuration OmnibusSigProc::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = m_anode_tn;
  cfg["ftoffset"] = m_fine_time_offset;
  cfg["ctoffset"] = m_coarse_time_offset;
  cfg["nticks"] = m_nticks;
  cfg["period"] = m_period;
  
  cfg["gain"] = m_gain;
  cfg["shaping_time"] = m_shaping_time;
  cfg["inter_gain"] = m_inter_gain;
  cfg["ADC_mV"] = m_ADC_mV;
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


void OmnibusSigProc::init_overall_response(){
  auto ifr = Factory::find<IFieldResponse>("FieldResponse");
  // Get full, "fine-grained" field responses defined at impact
  // positions.
  Response::Schema::FieldResponse fr = ifr->field_response();
  
  // Make a new data set which is the average FR
  Response::Schema::FieldResponse fravg = Response::wire_region_average(fr);

  WireCell::Waveform::compseq_t elec;
  WireCell::Binning tbins(Response::as_array(fravg.planes[0]).cols(), 0, Response::as_array(fravg.planes[0]).cols() * fravg.period);
  Response::ColdElec ce(m_gain, m_shaping_time);
  auto ewave = ce.generate(tbins);
  Waveform::scale(ewave, m_inter_gain * m_ADC_mV);
  elec = Waveform::dft(ewave);

  std::complex<float> fine_period(fravg.period,0);
  
  Waveform::realseq_t wfs(m_nticks);
  Waveform::realseq_t ctbins(m_nticks);
  for (int i=0;i!=m_nticks;i++){
    ctbins.at(i) = i * m_period;
  }

  int fine_nticks = Response::as_array(fravg.planes[0]).cols();
  
  Waveform::realseq_t ftbins(fine_nticks);
  for (int i=0;i!=fine_nticks;i++){
    ftbins.at(i) = i * fravg.period;
  }
  

  // Convert each average FR to a 2D array
  for (int ind=0; ind<3; ++ind) {
    auto arr = Response::as_array(fravg.planes[ind]);

    // do FFT for response ... 
    Array::array_xxc c_data = Array::dft_rc(arr,0);
    int nrows = c_data.rows();
    int ncols = c_data.cols();

    for (int irow = 0; irow < nrows; ++irow){
      for (int icol = 0; icol < ncols; ++ icol){
	c_data(irow,icol) = c_data(irow,icol) * elec.at(icol) * fine_period;
      }
    }
    
    arr = Array::idft_cr(c_data,0);

	
    // figure out how to do fine ... shift (good ...) 
    auto arr1 = arr.block(0,0,nrows,100);
    arr.block(0,0,nrows,ncols-100) = arr.block(0,100,nrows,ncols-100);
    arr.block(0,ncols-100,nrows,100) = arr1;
    
	
	
    // redigitize ... 
    for (int irow = 0; irow < nrows; ++ irow){
      // gtemp = new TGraph();
      
      int fcount = 1;
      for (int i=0;i!=m_nticks;i++){
	double ctime = ctbins.at(i);
	
	if (fcount < fine_nticks)
	  while(ctime > ftbins.at(fcount)){
	    fcount ++;
	    if (fcount >= fine_nticks) break;
	  }
	
	    
	if(fcount < fine_nticks){
	  wfs.at(i) = (ctime - ftbins.at(fcount-1)) / m_period * arr(irow,fcount-1) + (ftbins.at(fcount)-ctime)/m_period * arr(irow,fcount);
	}else{
	  wfs.at(i) = 0;
	}
      }
    } // loop inside wire ... 

  }// loop over plane
  
  
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
