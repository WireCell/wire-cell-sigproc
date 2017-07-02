#include "WireCellSigProc/OmnibusSigProc.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IFilterWaveform.h"
#include "WireCellIface/IChannelResponse.h"

#include "ROI_formation.h"
#include "ROI_refinement.h"

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusSigProc::OmnibusSigProc(const std::string anode_tn, double fine_time_offset, double coarse_time_offset, double period, int nticks, double gain , double shaping_time , double inter_gain , double ADC_mV, bool flag_ch_corr, float th_factor_ind, float th_factor_col, int pad, float asy, int rebin, double l_factor, double l_max_th, double l_factor1, int l_short_length )
  : m_anode_tn (anode_tn)
  , m_fine_time_offset(fine_time_offset)
  , m_coarse_time_offset(coarse_time_offset)
  , m_period(period)
  , m_nticks(nticks)
  , m_gain(gain)
  , m_shaping_time(shaping_time)
  , m_inter_gain(inter_gain)
  , m_ADC_mV(ADC_mV)
  , m_flag_ch_corr(flag_ch_corr)
  , m_th_factor_ind(th_factor_ind)
  , m_th_factor_col(th_factor_col)
  , m_pad(pad)
  , m_asy(asy)
  , m_rebin(rebin)
  , m_l_factor(l_factor)
  , m_l_max_th(l_max_th)
  , m_l_factor1(l_factor1)
  , m_l_short_length(l_short_length)
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

  m_flag_ch_corr = get(config,"ch_corr",m_flag_ch_corr);

  m_th_factor_ind = get(config,"troi_ind_th_factor",m_th_factor_ind);
  m_th_factor_col = get(config,"troi_col_th_factor",m_th_factor_col);
  m_pad = get(config,"troi_pad",m_pad);
  m_asy = get(config,"troi_asy",m_asy);
  m_rebin = get(config,"lroi_rebin",m_rebin);
  m_l_factor = get(config,"lroi_th_factor",m_l_factor);
  m_l_max_th = get(config,"lroi_max_th",m_l_max_th);
  m_l_factor1 = get(config,"lori_th_factor1",m_l_factor1);
  m_l_short_length = get(config,"lroi_short_length",m_l_short_length);
  
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

  cfg["ch_corr"] = m_flag_ch_corr;

  cfg["troi_ind_th_factor"] = m_th_factor_ind;
  cfg["troi_col_th_factor"] = m_th_factor_col;
  cfg["troi_pad"] = m_pad;
  cfg["troi_asy"] = m_asy;
  cfg["lroi_rebin"] = m_rebin; 
  cfg["lroi_th_factor"] = m_l_factor;
  cfg["lroi_max_th"] = m_l_max_th;
  cfg["lori_th_factor1"] = m_l_factor1;
  cfg["lroi_short_length"] = m_l_short_length; 

      
      
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
  int offset1=0;
  int nwire = 0;
  if (plane==0){
    offset = total_offset;
    offset1 = 0;
    nwire = nwire_u;
  }else if (plane==1){
    offset = nwire_u + total_offset;
    offset1 = nwire_u;
    nwire = nwire_v;
  }else if (plane==2){
    offset = nwire_u + nwire_v + total_offset;
    offset1 = nwire_u + nwire_v;
    nwire = nwire_w;
  }
  
  
  for (int i=0;i!=nwire;i++){
    ITrace::ChargeSequence charge(m_nticks);
    for (int j=0;j!=m_nticks;j++){
      charge.at(j) = r_data(i,j);
    }
    // correct the dead channels ... 
    if (cmm["bad"].find(i+offset1)!=cmm["bad"].end()){
      for (size_t k=0;k!=cmm["bad"][i+offset1].size();k++){
	for (int j=cmm["bad"][i+offset1].at(k).first; j!=cmm["bad"][i+offset1].at(k).second;j++){
	  charge.at(j)=0;
	}
      }
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

  // clear the overall response
  for (int i=0;i!=3;i++){
    overall_resp[i].clear();
  }

  m_intrinsic_time_offset = fr.origin/fr.speed;

  
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
    int fine_time_shift = m_fine_time_offset / fravg.period;
    if (fine_time_shift>0){
      Array::array_xxf arr1(nrows,ncols - fine_time_shift);
      arr1 = arr.block(0,0,nrows,ncols - fine_time_shift);
      Array::array_xxf arr2(nrows,fine_time_shift);
      arr2 = arr.block(0,ncols-fine_time_shift,nrows,fine_time_shift);
      arr.block(0,0,nrows,fine_time_shift) = arr2;
      arr.block(0,fine_time_shift,nrows,ncols-fine_time_shift) = arr1;
      
      // Array::array_xxf arr1(nrows,fine_time_shift);
      // arr1 = arr.block(0,0,nrows,fine_time_shift);
      // Array::array_xxf arr2(nrows,ncols-fine_time_shift);
      // arr2 = arr.block(0,fine_time_shift,nrows,ncols-fine_time_shift);
      // arr.block(0,0,nrows,ncols-fine_time_shift) = arr2;
      // arr.block(0,ncols-fine_time_shift,nrows,fine_time_shift) = arr1;
    }
	
	
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
	  wfs.at(i) = ((ctime - ftbins.at(fcount-1)) /fravg.period * arr(irow,fcount-1) + (ftbins.at(fcount)-ctime)/fravg.period * arr(irow,fcount)) / units::mV / (-1);
	}else{
	  wfs.at(i) = 0;
	}
      }
      
      overall_resp[ind].push_back(wfs);
      //wfs.clear();
    } // loop inside wire ...
    // calculated the wire shift ...     
    m_wire_shift[ind] = (int(overall_resp[ind].size())-1)/2;

    //    std::cout << /m_period << std::endl;
  }//  loop over plane

  
  
}

void OmnibusSigProc::restore_baseline(Array::array_xxf& arr){
  
  for (int i=0;i!=arr.rows();i++){
    Waveform::realseq_t signal(arr.cols());
    int ncount = 0;
    for (int j=0;j!=arr.cols();j++){
      if (arr(i,j)!=0){
	signal.at(ncount) = arr(i,j);
	ncount ++;
      }
    }
    signal.resize(ncount);
    float baseline = WireCell::Waveform::median_binned(signal);
    // std::cout << i << " " << baseline << " " << signal.size() << " ";
    Waveform::realseq_t temp_signal(arr.cols());
    ncount = 0;
    for (size_t j =0; j!=signal.size();j++){
      if (fabs(signal.at(j)-baseline) < 500){
	temp_signal.at(ncount) = signal.at(j);
	ncount ++;
      }
    }
    temp_signal.resize(ncount);
    
    baseline = WireCell::Waveform::median_binned(temp_signal);
    //std::cout << baseline << " " << temp_signal.size() << std::endl;
    
    for (int j=0;j!=arr.cols();j++){
      if (arr(i,j)!=0)
	arr(i,j) -= baseline;
    }
  }
}


void OmnibusSigProc::decon_2D_init(int plane){

  // data part ... 
  // first round of FFT on time
  c_data = Array::dft_rc(r_data,0);

  
  // now apply the ch-by-ch response ...
  //  bool flag_ch_corr = false;
  if (m_flag_ch_corr){
    auto cr = Factory::find<IChannelResponse>("PerChannelResponse");
    auto bins = cr->channel_response_binning();
    assert(bins.binsize()==m_period);
    //starndard electronics response ... 
    WireCell::Binning tbins(m_nticks, 0, m_nticks*m_period);
    Response::ColdElec ce(m_gain, m_shaping_time);
    auto ewave = ce.generate(tbins);
    //ch-by-ch electronics response
    int offset = 0;
    int nwire = 0;
    if (plane==0){
      offset = 0;
      nwire = nwire_u;
    }else if (plane==1){
      offset = nwire_u;
      nwire = nwire_v;
    }else{
      offset = nwire_u + nwire_v;
      nwire = nwire_w;
    }
    std::vector<Waveform::realseq_t> ch_wfs;
    ch_wfs.resize(nwire);
    for (int ich=0;ich!=nwire;ich++){
      ch_wfs.at(ich).resize(m_nticks,0);
      auto& resp = cr->channel_response(ich+offset);
      for (size_t i=0;i!=resp.size();i++){
	ch_wfs.at(ich).at(i) = resp.at(i);
      }
      //std::cout << ich << " " << resp.size() << std::endl;
    }
    
    WireCell::Waveform::compseq_t elec = Waveform::dft(ewave);
    for (int irow = 0; irow != c_data.rows(); irow++){
      WireCell::Waveform::compseq_t ch_elec = Waveform::dft(ch_wfs.at(irow));
      for (int icol = 0; icol != c_data.cols(); icol++){
	if (abs(ch_elec.at(icol))!=0){
	  c_data(irow,icol) *= elec.at(icol) / ch_elec.at(icol);
	}else{
	  c_data(irow,icol) = 0;
	}
      }
    }
  }


  
  
  //second round of FFT on wire
  c_data = Array::dft_cc(c_data,1);
  
  //response part ...
  Array::array_xxf r_resp = Array::array_xxf::Zero(r_data.rows(),m_nticks);
  for (size_t i=0;i!=overall_resp[plane].size();i++){
    for (int j=0;j!=m_nticks;j++){
      r_resp(i,j) = overall_resp[plane].at(i).at(j);
    }
  }
  
  // do first round FFT on the resposne on time
  Array::array_xxc c_resp = Array::dft_rc(r_resp,0);
  // do second round FFT on the response on wire
  c_resp = Array::dft_cc(c_resp,1);

  
  //make ratio to the response and apply wire filter
  c_data = c_data/c_resp;
  //std::cout << "Apply Wire Filter! " << std::endl;
  // apply software filter on wire
  Waveform::realseq_t wire_filter_wf;
  if (plane ==0 || plane == 1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wire_ind");
    wire_filter_wf = ncr1->filter_waveform();
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wire_col");
    wire_filter_wf = ncr1->filter_waveform();
  }
  for (int irow=0; irow<c_data.rows(); ++irow) {
    for (int icol=0; icol<c_data.cols(); ++icol) {
      float val = abs(c_data(irow,icol));
      if (std::isnan(val)) {
	c_data(irow,icol) = -0.0;
      }
      if (std::isinf(val)) {
	c_data(irow,icol) = 0.0;
      }
      c_data(irow,icol) *= wire_filter_wf.at(irow);
    }
  }
  
  //do the first round of inverse FFT on wire
  c_data = Array::idft_cc(c_data,1);

  // do the second round of inverse FFT on time
  r_data = Array::idft_cr(c_data,0);

  // do the shift in wire 
  int nrows = r_data.rows();
  int ncols = r_data.cols();
  {    // std::cout << nrows << " " << ncols << " " << m_wire_shift[plane] << std::endl;
    Array::array_xxf arr1(m_wire_shift[plane], ncols) ;
    arr1 = r_data.block(nrows-m_wire_shift[plane] , 0 , m_wire_shift[plane], ncols);
    Array::array_xxf arr2(nrows-m_wire_shift[plane],ncols);
    arr2 = r_data.block(0,0,nrows-m_wire_shift[plane],ncols);
    r_data.block(0,0,m_wire_shift[plane],ncols) = arr1;
    r_data.block(m_wire_shift[plane],0,nrows-m_wire_shift[plane],ncols) = arr2;
  }
  
  //do the shift in time
  int time_shift = (m_coarse_time_offset + m_intrinsic_time_offset)/m_period;
  if (time_shift > 0){
    Array::array_xxf arr1(nrows,ncols - time_shift);
    arr1 = r_data.block(0,0,nrows,ncols - time_shift);
    Array::array_xxf arr2(nrows,time_shift);
    arr2 = r_data.block(0,ncols-time_shift,nrows,time_shift);
    r_data.block(0,0,nrows,time_shift) = arr2;
    r_data.block(0,time_shift,nrows,ncols-time_shift) = arr1;
  }
  c_data = Array::dft_rc(r_data,0);
  
  
  
}

void OmnibusSigProc::decon_2D_ROI_refine(int plane){
   // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_U");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_V");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_W");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }

  Array::array_xxc c_data_afterfilter(r_data.rows(),r_data.cols());
  for (int irow=0; irow<c_data.rows(); ++irow) {
    for (int icol=0; icol<c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(r_data);
}

void OmnibusSigProc::decon_2D_tightROI(int plane){
  // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_U");
    roi_hf_filter_wf = ncr1->filter_waveform();
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
    auto temp_filter = ncr2->filter_waveform();
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_V");
    roi_hf_filter_wf = ncr1->filter_waveform();
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
    auto temp_filter = ncr2->filter_waveform();
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_W");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }

  Array::array_xxc c_data_afterfilter(r_data.rows(),r_data.cols());
  for (int irow=0; irow<c_data.rows(); ++irow) {
    for (int icol=0; icol<c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(r_data);
}
 
void OmnibusSigProc::decon_2D_tighterROI(int plane){
  // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_U");
    roi_hf_filter_wf = ncr1->filter_waveform();
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tighter_lf");
    auto temp_filter = ncr2->filter_waveform();
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_V");
    roi_hf_filter_wf = ncr1->filter_waveform();
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tighter_lf");
    auto temp_filter = ncr2->filter_waveform();
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_W");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }

  Array::array_xxc c_data_afterfilter(r_data.rows(),r_data.cols());
  for (int irow=0; irow<c_data.rows(); ++irow) {
    for (int icol=0; icol<c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(r_data);
}
 

void OmnibusSigProc::decon_2D_looseROI(int plane){
  if (plane == 2) return;
   // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  Waveform::realseq_t roi_hf_filter_wf1;
  Waveform::realseq_t roi_hf_filter_wf2;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_U");
    roi_hf_filter_wf = ncr1->filter_waveform();
    roi_hf_filter_wf1 = roi_hf_filter_wf;
    {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_loose_lf");
      auto temp_filter = ncr2->filter_waveform();
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf.at(i) *= temp_filter.at(i);
      }
    }
    {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
      auto temp_filter = ncr2->filter_waveform();
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf1.at(i) *= temp_filter.at(i);
      }
    }
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_V");
    roi_hf_filter_wf = ncr1->filter_waveform();
    roi_hf_filter_wf1 = roi_hf_filter_wf;
    {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_loose_lf");
      auto temp_filter = ncr2->filter_waveform();
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf.at(i) *= temp_filter.at(i);
      }
    }
     {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
      auto temp_filter = ncr2->filter_waveform();
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf1.at(i) *= temp_filter.at(i);
      }
    }
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_W");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }

  Array::array_xxc c_data_afterfilter(r_data.rows(),r_data.cols());
  for (int irow=0; irow<c_data.rows(); ++irow) {
    if (plane == 0 || plane == 1){
      if (cmm["lf_noisy"].find(nwire_u*plane+irow)==cmm["lf_noisy"].end()){
	roi_hf_filter_wf2 = roi_hf_filter_wf;
      }else{
	roi_hf_filter_wf2 = roi_hf_filter_wf1;
      }
      //Waveform::realseq_t roi_hf_filter_wf;
    }else{
      roi_hf_filter_wf2 = roi_hf_filter_wf;
    }
    
    for (int icol=0; icol<c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = c_data(irow,icol) * roi_hf_filter_wf2.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(r_data);
}

void OmnibusSigProc::decon_2D_hits(int plane){
   // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_wide_U");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_wide_V");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_wide_W");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }

  Array::array_xxc c_data_afterfilter(r_data.rows(),r_data.cols());
  for (int irow=0; irow<c_data.rows(); ++irow) {
    for (int icol=0; icol<c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  r_data = Array::idft_cr(c_data_afterfilter,0);
  if (plane==2)
    restore_baseline(r_data);
}

void OmnibusSigProc::decon_2D_charge(int plane){
  // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Gaus_wide");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Gaus_wide");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Gaus_wide");
    roi_hf_filter_wf = ncr1->filter_waveform();
  }

  Array::array_xxc c_data_afterfilter(r_data.rows(),r_data.cols());
  for (int irow=0; irow<c_data.rows(); ++irow) {
    for (int icol=0; icol<c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  r_data = Array::idft_cr(c_data_afterfilter,0);
  if (plane==2)
    restore_baseline(r_data);
}


bool OmnibusSigProc::operator()(const input_pointer& in, output_pointer& out)
{
  cmm = in->masks();
  ITrace::vector itraces;

  // initialize the overall response function ... 
  init_overall_response();

  // create a class for ROIs ... 
  ROI_formation roi_form(cmm,nwire_u, nwire_v, nwire_w, m_nticks, m_th_factor_ind, m_th_factor_col, m_pad, m_asy, m_rebin, m_l_factor, m_l_max_th, m_l_factor1, m_l_short_length);
  ROI_refinement roi_refine(nwire_u,nwire_v,nwire_w);
  
  for (int i=0;i!=3;i++){
    // load data into EIGEN matrices ...
    load_data(in,i);
    // initial decon ... 
    decon_2D_init(i);

    
    // Form tight ROIs
    if (i!=2){ // induction wire planes
      decon_2D_tighterROI(i);
      Array::array_xxf r_data_tight(r_data.rows(),r_data.cols());
      r_data_tight = r_data;
      decon_2D_tightROI(i);
      roi_form.find_ROI_by_decon_itself(i,r_data,r_data_tight);
    }else{ // collection wire planes
      decon_2D_tightROI(i);
      roi_form.find_ROI_by_decon_itself(i,r_data);
    }
    
    // Form loose ROIs
    if (i!=2){
      decon_2D_looseROI(i);
      roi_form.find_ROI_loose(i,r_data);
      decon_2D_ROI_refine(i);
    }
    // Refine ROIs
    roi_refine.load_data(i,r_data,roi_form);


    // merge results ...
    decon_2D_hits(i);
    decon_2D_charge(i);
    
    
    // Get results
    save_data(itraces,i);
  }

  // put threshold into cmm ...
  std::vector<float>& u_rms = roi_form.get_uplane_rms();
  std::vector<float>& v_rms = roi_form.get_vplane_rms();
  std::vector<float>& w_rms = roi_form.get_wplane_rms();
  Waveform::ChannelMasks temp;
  cmm["threshold"] = temp;
  for ( size_t i=0;i!=u_rms.size();i++){
    cmm["threshold"][i].push_back(std::make_pair(int(u_rms.at(i)*10), 10));
  }
  for ( size_t i=0;i!=v_rms.size();i++){
    cmm["threshold"][i+nwire_u].push_back(std::make_pair(int(v_rms.at(i)*10), 10));
  }
  for ( size_t i=0;i!=w_rms.size();i++){
    cmm["threshold"][i+nwire_u+nwire_v].push_back(std::make_pair(int(w_rms.at(i)*10), 10));
  }
  
  
  
  SimpleFrame* sframe = new SimpleFrame(in->ident(), in->time(), itraces, in->tick(), cmm);
  out = IFrame::pointer(sframe);
  
  return true;
}
