#include "WireCellSigProc/OmnibusSigProc.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IFilterWaveform.h"
#include "WireCellIface/IChannelResponse.h"

#include "ROI_formation.h"
#include "ROI_refinement.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(OmnibusSigProc, WireCell::SigProc::OmnibusSigProc,
                 WireCell::IFrameFilter, WireCell::IConfigurable);


using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusSigProc::OmnibusSigProc(const std::string& anode_tn,
                               const std::string& per_chan_resp_tn,
                               const std::string& field_response,
                               double fine_time_offset,
                               double coarse_time_offset,
                               double gain, 
                               double shaping_time,
                               double inter_gain , 
                               double ADC_mV,
                               float th_factor_ind,
                               float th_factor_col,
                               int pad,
                               float asy,
                               int rebin,
                               double l_factor,
                               double l_max_th,
                               double l_factor1,
                               int l_short_length,
                               double r_th_factor ,
                               double r_fake_signal_low_th ,
                               double r_fake_signal_high_th,
                               int r_pad ,
                               int r_break_roi_loop ,
                               double r_th_peak ,
                               double r_sep_peak,
                               double r_low_peak_sep_threshold_pre ,
                               int r_max_npeaks ,
                               double r_sigma ,
                               double r_th_percent ,
                               int charge_ch_offset  )
  : m_anode_tn (anode_tn)
  , m_per_chan_resp(per_chan_resp_tn)
  , m_field_response(field_response)
  , m_fine_time_offset(fine_time_offset)
  , m_coarse_time_offset(coarse_time_offset)
  , m_period(0)
  , m_nticks(0)
  , m_gain(gain)
  , m_shaping_time(shaping_time)
  , m_inter_gain(inter_gain)
  , m_ADC_mV(ADC_mV)
  , m_th_factor_ind(th_factor_ind)
  , m_th_factor_col(th_factor_col)
  , m_pad(pad)
  , m_asy(asy)
  , m_rebin(rebin)
  , m_l_factor(l_factor)
  , m_l_max_th(l_max_th)
  , m_l_factor1(l_factor1)
  , m_l_short_length(l_short_length)
  , m_r_th_factor(r_th_factor)
  , m_r_fake_signal_low_th(r_fake_signal_low_th)
  , m_r_fake_signal_high_th(r_fake_signal_high_th)
  , m_r_pad(r_pad)
  , m_r_break_roi_loop(r_break_roi_loop)
  , m_r_th_peak(r_th_peak)
  , m_r_sep_peak(r_sep_peak)
  , m_r_low_peak_sep_threshold_pre(r_low_peak_sep_threshold_pre)
  , m_r_max_npeaks(r_max_npeaks)
  , m_r_sigma(r_sigma)
  , m_r_th_percent(r_th_percent)
  , m_charge_ch_offset(charge_ch_offset)
{
  // get wires for each plane

 
  //std::cout << m_anode->channels().size() << " " << nwire_u << " " << nwire_v << " " << nwire_w << std::endl;
  
}

OmnibusSigProc::~OmnibusSigProc()
{
}

std::string WireCell::SigProc::OmnibusSigProc::OspChan::str() const
{
  std::stringstream ss;
  ss<<"OspChan<c:"<<channel<<",w:"<<wire<<",p:"<<plane<<",i:"<<ident<<">";
  return ss.str();
}


void OmnibusSigProc::configure(const WireCell::Configuration& config)
{
  m_fine_time_offset = get(config,"ftoffset",m_fine_time_offset);
  m_coarse_time_offset = get(config,"ctoffset",m_coarse_time_offset);
  m_anode_tn = get(config, "anode", m_anode_tn);
  m_nticks = get(config,"nticks",m_nticks);
  m_period = get(config,"period",m_period);

  m_gain = get(config,"gain",m_gain);
  m_shaping_time = get(config,"shaping",m_shaping_time);
  m_inter_gain = get(config,"postgain", m_inter_gain);
  m_ADC_mV = get(config,"ADC_mV", m_ADC_mV);

  m_per_chan_resp = get(config, "per_chan_resp", m_per_chan_resp);
  m_field_response = get(config, "field_response", m_field_response);

  m_th_factor_ind = get(config,"troi_ind_th_factor",m_th_factor_ind);
  m_th_factor_col = get(config,"troi_col_th_factor",m_th_factor_col);
  m_pad = get(config,"troi_pad",m_pad);
  m_asy = get(config,"troi_asy",m_asy);
  m_rebin = get(config,"lroi_rebin",m_rebin);
  m_l_factor = get(config,"lroi_th_factor",m_l_factor);
  m_l_max_th = get(config,"lroi_max_th",m_l_max_th);
  m_l_factor1 = get(config,"lori_th_factor1",m_l_factor1);
  m_l_short_length = get(config,"lroi_short_length",m_l_short_length);


  m_r_th_factor = get(config,"r_th_factor",m_r_th_factor);
  m_r_fake_signal_low_th = get(config,"r_fake_signal_low_th",m_r_fake_signal_low_th);
  m_r_fake_signal_high_th = get(config,"r_fake_signal_high_th",m_r_fake_signal_high_th);
  m_r_pad = get(config,"r_pad",m_r_pad);
  m_r_break_roi_loop = get(config,"r_break_roi_loop",m_r_break_roi_loop);
  m_r_th_peak = get(config,"r_th_peak",m_r_th_peak);
  m_r_sep_peak = get(config,"r_sep_peak",m_r_sep_peak);
  m_r_low_peak_sep_threshold_pre = get(config,"r_low_peak_sep_threshold_pre",m_r_low_peak_sep_threshold_pre);
  m_r_max_npeaks = get(config,"r_max_npeaks",m_r_max_npeaks);
  m_r_sigma = get(config,"r_sigma",m_r_sigma);
  m_r_th_percent = get(config,"r_th_percent",m_r_th_percent);

  m_charge_ch_offset = get(config,"charge_ch_offset",m_charge_ch_offset);
  
  // this throws if not found
  m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

  // Build up the channel map.  The OSP channel must run contiguously
  // first up the U, then V, then W "wires".  Ie, face-major order,
  // but we have plane-major order so make a temporary collection.
  IChannel::vector plane_channels[3];
  for (auto face : m_anode->faces()) {
    if (!face) { // A null face means one sided AnodePlane.  
      continue;  // Can be "back" or "front" face.
    }
    for (auto plane: face->planes()) {
      int plane_index = plane->planeid().index();
      auto& pchans = plane_channels[plane_index];
      // These IChannel vectors are ordered in same order as wire-in-plane.
      const auto& ichans = plane->channels();
      // Append
      pchans.reserve(pchans.size() + ichans.size());
      pchans.insert(pchans.end(), ichans.begin(), ichans.end());
      std::cerr << "a" << m_anode->ident() << " "
                << "f" << face->ident() << " "
                << "p" << plane->ident() << " "
                << "c" << ichans.front()->ident() << " -> c" << ichans.back()->ident() << ", "
                << "i" << ichans.front()->index() << " -> i" << ichans.back()->index() << "\n";
    }
  }

  int osp_channel_number = 0;
  for (int iplane = 0; iplane < 3; ++iplane) {
    m_nwires[iplane] = plane_channels[iplane].size();
    int osp_wire_number = 0;
    for (auto ichan : plane_channels[iplane]) {
      const int wct_chan_ident = ichan->ident();
      OspChan och(osp_channel_number, osp_wire_number, iplane, wct_chan_ident);
      m_channel_map[wct_chan_ident] = och; // we could save some space by storing
      m_channel_range[iplane].push_back(och);// wct ident here instead of a whole och.
      ++osp_wire_number;
      ++osp_channel_number;
    }
    // std::cerr << iplane << ": och:["
    //           << m_channel_range[iplane].front().str()
    //           << " -> "
    //           << m_channel_range[iplane].back().str()
    //           << std::endl;
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

  cfg["per_chan_resp"] = m_per_chan_resp;
  cfg["field_response"] = m_field_response;

  cfg["troi_ind_th_factor"] = m_th_factor_ind;
  cfg["troi_col_th_factor"] = m_th_factor_col;
  cfg["troi_pad"] = m_pad;
  cfg["troi_asy"] = m_asy;
  cfg["lroi_rebin"] = m_rebin; 
  cfg["lroi_th_factor"] = m_l_factor;
  cfg["lroi_max_th"] = m_l_max_th;
  cfg["lori_th_factor1"] = m_l_factor1;
  cfg["lroi_short_length"] = m_l_short_length; 

  cfg["r_th_factor"] = m_r_th_factor;
  cfg["r_fake_signal_low_th"] = m_r_fake_signal_low_th;
  cfg["r_fake_signal_high_th"] = m_r_fake_signal_high_th;
  cfg["r_pad"] = m_r_pad;
  cfg["r_break_roi_loop"] = m_r_break_roi_loop;
  cfg["r_th_peak"] = m_r_th_peak;
  cfg["r_sep_peak"] = m_r_sep_peak;
  cfg["r_low_peak_sep_threshold_pre"] = m_r_low_peak_sep_threshold_pre;
  cfg["r_max_npeaks"] = m_r_max_npeaks;
  cfg["r_sigma"] = m_r_sigma;
  cfg["r_th_precent"] = m_r_th_percent;
      
  cfg["charge_ch_offset"] = m_charge_ch_offset;
  
  return cfg;
  
}

void OmnibusSigProc::load_data(const input_pointer& in, int plane){
  m_r_data = Array::array_xxf::Zero(m_nwires[plane],m_nticks);

  auto traces = in->traces();

  auto& bad = m_cmm["bad"];
  int nbad = 0;

  for (auto trace : *traces.get()) {
    int wct_channel_ident = trace->channel();
    OspChan och = m_channel_map[wct_channel_ident];
    if (plane != och.plane) {
      continue;         // we'll catch it in another call to load_data
    }

    // fixme: this code uses tbin() but other places in this file will barf if tbin!=0.
    int tbin = trace->tbin();
    auto const& charges = trace->charge();
    const int ntbins = std::min((int)charges.size(), m_nticks);
    for (int qind = 0; qind < ntbins; ++qind) {
      m_r_data(och.wire, tbin + qind) = charges[qind];
    }

    //ensure dead channels are indeed dead ...
    auto const& badch = bad.find(och.channel);
    if (badch == bad.end()) {
      continue;
    }

    auto const& binranges = badch->second;
    for (auto const& br : binranges) {
      ++nbad;
      for (int i = br.first; i != br.second; ++i) {
        m_r_data(och.wire, i) = 0;
      }
      //std::cerr << plane << " " << ch << ": [" << br.first << "," << br.second << "]\n";
    }

  }
  std::cerr << "OmnibusSigProc: plane index: " << plane << " configured with " << nbad << " bad regions\n";
  
}


// Note, this MUST fill itraces in OSP channel/wire order
void OmnibusSigProc::save_data(ITrace::vector& itraces, IFrame::trace_list_t& indices, int plane)
{
  double qtot = 0.0;
  for (auto och : m_channel_range[plane]) { // ordered by osp channel
    
    ITrace::ChargeSequence charge(m_nticks);
    for (int itick=0;itick!=m_nticks;itick++){
      const float q = m_r_data(och.wire, itick);
      charge.at(itick) = q > 0.0 ? q : 0.0;
    }

    // zero the bad channels ... 
    {
      auto& bad = m_cmm["bad"];
      auto badit = bad.find(och.channel);
      if (badit != bad.end()) {
        for (auto bad : badit->second) {
          for (int itick=bad.first; itick < bad.second; ++itick) {
            charge.at(itick) = 0.0;
          }
        }
      }
    }

    // debug
    for (int j=0;j!=m_nticks;j++){
      qtot += charge.at(j);
    }

    SimpleTrace *trace = new SimpleTrace(och.ident, 0, charge);
    const size_t trace_index = itraces.size();
    indices.push_back(trace_index);
    itraces.push_back(ITrace::pointer(trace));
  }
  const int nadded = indices.back() - indices.front() + 1;
  std::cerr << "OmnibusSigProc::save_data plane index: " << plane << " Qtot=" << qtot
            << " added " << nadded << " traces to total " << indices.size()
            << " indices:[" << indices.front() << "," << indices.back() << "]\n";
  
}


void OmnibusSigProc::init_overall_response(IFrame::pointer frame)
{
  m_period = frame->tick();
  {
    std::vector<int> tbins;
    for (auto trace : *frame->traces()) {
      const int tbin = trace->tbin();
      const int nbins = trace->charge().size();
      tbins.push_back(tbin);
      tbins.push_back(tbin+nbins);
    }
    auto mme = std::minmax_element(tbins.begin(), tbins.end());
    int tbinmin = *mme.first;
    int tbinmax = *mme.second;
    m_nticks = tbinmax-tbinmin;
    std::cerr <<"OmnibusSigProc: nticks=" << m_nticks << " tbinmin="<<tbinmin << " tbinmax="<<tbinmax<<std::endl;
  }

  auto ifr = Factory::find_tn<IFieldResponse>(m_field_response);
  // Get full, "fine-grained" field responses defined at impact
  // positions.
  Response::Schema::FieldResponse fr = ifr->field_response();
  
  // Make a new data set which is the average FR
  Response::Schema::FieldResponse fravg = Response::wire_region_average(fr);

  WireCell::Waveform::compseq_t elec;
  WireCell::Binning tbins(Response::as_array(fravg.planes[0]).cols(), 0, Response::as_array(fravg.planes[0]).cols() * fravg.period);
  Response::ColdElec ce(m_gain, m_shaping_time);
  auto ewave = ce.generate(tbins);
  Waveform::scale(ewave, m_inter_gain * m_ADC_mV  * (-1));
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
  for (int iplane=0; iplane<3; ++iplane) {
    auto arr = Response::as_array(fravg.planes[iplane]);

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
	  wfs.at(i) = ((ctime - ftbins.at(fcount-1)) /fravg.period * arr(irow,fcount-1) + (ftbins.at(fcount)-ctime)/fravg.period * arr(irow,fcount)) ;// / (-1);
	}else{
	  wfs.at(i) = 0;
	}
      }
      
      overall_resp[iplane].push_back(wfs);
      //wfs.clear();
    } // loop inside wire ...
    // calculated the wire shift ...     
    m_wire_shift[iplane] = (int(overall_resp[iplane].size())-1)/2;

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
  m_c_data = Array::dft_rc(m_r_data,0);

  
  // now apply the ch-by-ch response ...
  if (! m_per_chan_resp.empty()) {
    std::cerr<<"OmnibusSigProc: CH-BY-CH ELECTRONICS RESPONSE CORRECTION\n";
    auto cr = Factory::find_tn<IChannelResponse>(m_per_chan_resp);
    auto cr_bins = cr->channel_response_binning();
    assert(cr_bins.binsize()==m_period);
    //starndard electronics response ... 
    // WireCell::Binning tbins(m_nticks, 0-m_period/2., m_nticks*m_period-m_period/2.);
    // Response::ColdElec ce(m_gain, m_shaping_time);

    // temporary hack ...
    //float scaling = 1./(1e-9*0.5/1.13312);
    //WireCell::Binning tbins(m_nticks, (-5-0.5)*m_period, (m_nticks-5-0.5)*m_period-m_period);
    //Response::ColdElec ce(m_gain*scaling, m_shaping_time);
    //// this is moved into wirecell.sigproc.main production of
    //// microboone-channel-responses-v1.json.bz2
    WireCell::Binning tbins(m_nticks, cr_bins.min(), cr_bins.min() + m_nticks*m_period);
    Response::ColdElec ce(m_gain, m_shaping_time);
    
    const auto ewave = ce.generate(tbins);
    const WireCell::Waveform::compseq_t elec = Waveform::dft(ewave);

    for (auto och : m_channel_range[plane]) {
      const auto& ch_resp = cr->channel_response(och.ident);
      const WireCell::Waveform::compseq_t ch_elec = Waveform::dft(ch_resp);

      const int irow = och.wire;
      for (int icol = 0; icol != m_c_data.cols(); icol++){
        const auto four = ch_elec.at(icol);
	if (std::abs(four) != 0){
	  m_c_data(irow,icol) *= elec.at(icol) / four;
	}else{
	  m_c_data(irow,icol) = 0;
	}
      }
    }
  }


  
  
  //second round of FFT on wire
  m_c_data = Array::dft_cc(m_c_data,1);
  
  //response part ...
  Array::array_xxf r_resp = Array::array_xxf::Zero(m_r_data.rows(),m_nticks);
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
  m_c_data = m_c_data/c_resp;
  //std::cout << "Apply Wire Filter! " << std::endl;
  // apply software filter on wire
  const std::vector<std::string> filter_names{"Wire_ind", "Wire_ind", "Wire_col"};
  Waveform::realseq_t wire_filter_wf;
  auto ncr1 = Factory::find<IFilterWaveform>("HfFilter", filter_names[plane]);
  wire_filter_wf = ncr1->filter_waveform(m_c_data.rows());
  for (int irow=0; irow<m_c_data.rows(); ++irow) {
    for (int icol=0; icol<m_c_data.cols(); ++icol) {
      float val = abs(m_c_data(irow,icol));
      if (std::isnan(val)) {
	m_c_data(irow,icol) = -0.0;
      }
      if (std::isinf(val)) {
	m_c_data(irow,icol) = 0.0;
      }
      m_c_data(irow,icol) *= wire_filter_wf.at(irow);
    }
  }
  
  //do the first round of inverse FFT on wire
  m_c_data = Array::idft_cc(m_c_data,1);

  // do the second round of inverse FFT on time
  m_r_data = Array::idft_cr(m_c_data,0);

  // do the shift in wire 
  const int nrows = m_r_data.rows();
  const int ncols = m_r_data.cols();
  {    // std::cout << nrows << " " << ncols << " " << m_wire_shift[plane] << std::endl;
    Array::array_xxf arr1(m_wire_shift[plane], ncols) ;
    arr1 = m_r_data.block(nrows-m_wire_shift[plane] , 0 , m_wire_shift[plane], ncols);
    Array::array_xxf arr2(nrows-m_wire_shift[plane],ncols);
    arr2 = m_r_data.block(0,0,nrows-m_wire_shift[plane],ncols);
    m_r_data.block(0,0,m_wire_shift[plane],ncols) = arr1;
    m_r_data.block(m_wire_shift[plane],0,nrows-m_wire_shift[plane],ncols) = arr2;
  }
  
  //do the shift in time
  int time_shift = (m_coarse_time_offset + m_intrinsic_time_offset)/m_period;
  if (time_shift > 0){
    Array::array_xxf arr1(nrows,ncols - time_shift);
    arr1 = m_r_data.block(0,0,nrows,ncols - time_shift);
    Array::array_xxf arr2(nrows,time_shift);
    arr2 = m_r_data.block(0,ncols-time_shift,nrows,time_shift);
    m_r_data.block(0,0,nrows,time_shift) = arr2;
    m_r_data.block(0,time_shift,nrows,ncols-time_shift) = arr1;
  }
  m_c_data = Array::dft_rc(m_r_data,0);
}


void OmnibusSigProc::decon_2D_ROI_refine(int plane){
   // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  const std::vector<std::string> filter_names{"Wiener_tight_U", "Wiener_tight_V", "Wiener_tight_W"};
  Waveform::realseq_t roi_hf_filter_wf;

  auto ncr1 = Factory::find<IFilterWaveform>("HfFilter", filter_names[plane]);
  roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());

  Array::array_xxc c_data_afterfilter(m_r_data.rows(),m_r_data.cols());
  for (int irow=0; irow<m_c_data.rows(); ++irow) {
    for (int icol=0; icol<m_c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = m_c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  m_r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(m_r_data);
}


void OmnibusSigProc::decon_2D_tightROI(int plane){
  // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_U");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
    auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_V");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
    auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_W");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }

  Array::array_xxc c_data_afterfilter(m_r_data.rows(),m_r_data.cols());
  for (int irow=0; irow<m_c_data.rows(); ++irow) {
    for (int icol=0; icol<m_c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = m_c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  m_r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(m_r_data);
}
 
// same as above but with "tight" -> "tighter" for ROI filterss.
void OmnibusSigProc::decon_2D_tighterROI(int plane){
  // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_U");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tighter_lf");
    auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_V");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
    auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tighter_lf");
    auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
    for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
      roi_hf_filter_wf.at(i) *= temp_filter.at(i);
    }
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_W");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }

  Array::array_xxc c_data_afterfilter(m_r_data.rows(),m_r_data.cols());
  for (int irow=0; irow<m_c_data.rows(); ++irow) {
    for (int icol=0; icol<m_c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = m_c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  m_r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(m_r_data);
}
 

void OmnibusSigProc::decon_2D_looseROI(int plane){
  if (plane == 2) {
    return;                     // don't filter colleciton
  }

   // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  Waveform::realseq_t roi_hf_filter_wf1;
  Waveform::realseq_t roi_hf_filter_wf2;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_U");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
    roi_hf_filter_wf1 = roi_hf_filter_wf;
    {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_loose_lf");
      auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf.at(i) *= temp_filter.at(i);
      }
    }
    {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
      auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf1.at(i) *= temp_filter.at(i);
      }
    }
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_V");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
    roi_hf_filter_wf1 = roi_hf_filter_wf;
    {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_loose_lf");
      auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf.at(i) *= temp_filter.at(i);
      }
    }
     {
      auto ncr2 = Factory::find<IFilterWaveform>("LfFilter","ROI_tight_lf");
      auto temp_filter = ncr2->filter_waveform(m_c_data.cols());
      for(size_t i=0;i!=roi_hf_filter_wf.size();i++){
	roi_hf_filter_wf1.at(i) *= temp_filter.at(i);
      }
    }
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_tight_W");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }

  const int n_lfn_nn = 2;
  const int n_bad_nn = plane ? 1 : 2;

  Array::array_xxc c_data_afterfilter(m_r_data.rows(),m_r_data.cols());
  for (auto och : m_channel_range[plane]) {
    const int irow = och.wire;

    roi_hf_filter_wf2 = roi_hf_filter_wf;
    if (masked_neighbors("bad", och, n_bad_nn) or
        masked_neighbors("lf_noisy", och, n_lfn_nn))
    {
      roi_hf_filter_wf2 = roi_hf_filter_wf1;
    }
    
    for (int icol=0; icol<m_c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = m_c_data(irow,icol) * roi_hf_filter_wf2.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  m_r_data = Array::idft_cr(c_data_afterfilter,0);
  restore_baseline(m_r_data);
}


// return true if any channels w/in +/- nnn, inclusive, of the channel has the mask.
bool OmnibusSigProc::masked_neighbors(const std::string& cmname, OspChan& ochan, int nnn)
{
  // take care of boundary cases
  int lo_wire = ochan.wire - nnn;
  int lo_chan = ochan.channel - nnn;
  while (lo_wire < 0) {
    ++lo_wire;
    ++lo_chan;
  }
  const int nwires = m_nwires[ochan.plane];
  int hi_wire = ochan.wire + nnn;
  int hi_chan = ochan.channel + nnn;
  while (hi_wire >= nwires) {
    --hi_wire;
    --hi_chan;
  }
  if (hi_chan < lo_chan) {      // how?  bogus inputs?
    return false;              
  }

  auto& cm = m_cmm[cmname];
  for (int och = lo_chan; och <= hi_chan; ++och) {
    if (cm.find(och) != cm.end()) {
      return true;
    }
  }
  return false;
}

void OmnibusSigProc::decon_2D_hits(int plane){
   // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_wide_U");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_wide_V");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Wiener_wide_W");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }

  Array::array_xxc c_data_afterfilter(m_r_data.rows(),m_r_data.cols());
  for (int irow=0; irow<m_c_data.rows(); ++irow) {
    for (int icol=0; icol<m_c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = m_c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  m_r_data = Array::idft_cr(c_data_afterfilter,0);
  if (plane==2) {
    restore_baseline(m_r_data);
  }
}

void OmnibusSigProc::decon_2D_charge(int plane){
  // apply software filter on time
  //std::cout << "Apply Time Filter! " << std::endl;
  Waveform::realseq_t roi_hf_filter_wf;
  if (plane ==0){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Gaus_wide");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }else if (plane==1){
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Gaus_wide");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }else{
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter","Gaus_wide");
    roi_hf_filter_wf = ncr1->filter_waveform(m_c_data.cols());
  }

  Array::array_xxc c_data_afterfilter(m_r_data.rows(),m_r_data.cols());
  for (int irow=0; irow<m_c_data.rows(); ++irow) {
    for (int icol=0; icol<m_c_data.cols(); ++icol) {
      c_data_afterfilter(irow,icol) = m_c_data(irow,icol) * roi_hf_filter_wf.at(icol);
    }
  }
  
  //do the second round of inverse FFT on wire
  m_r_data = Array::idft_cr(c_data_afterfilter,0);
  if (plane==2) {
    restore_baseline(m_r_data);
  }
}


bool OmnibusSigProc::operator()(const input_pointer& in, output_pointer& out)
{
  if (!in) {
    // A null input indicates end of stream and is to let us flush
    // any data we may have buffered.  Since we do not buffer,
    // just return.
    out = nullptr;
    return true;
  }

  // Convert to OSP cmm indexed by OSB sequential channels, NOT WCT channel ID.
  m_cmm.clear();
  // double emap: name -> channel -> pair<int,int>
  for (auto cm : in->masks()) {
    const std::string name = cm.first;
    for (auto m: cm.second) {
      const int wct_channel_ident = m.first;
      const OspChan& och = m_channel_map[wct_channel_ident];
      m_cmm[name][och.channel] = m.second;
      //std::cerr << wct_channel_ident << " " << och.str() << std::endl;
    }
  }

  ITrace::vector itraces;
  IFrame::trace_list_t wiener_traces, gauss_traces, perframe_traces[3];

  // initialize the overall response function ... 
  init_overall_response(in);

  // create a class for ROIs ... 
  ROI_formation roi_form(m_cmm, m_nwires[0], m_nwires[1], m_nwires[2], m_nticks, m_th_factor_ind, m_th_factor_col, m_pad, m_asy, m_rebin, m_l_factor, m_l_max_th, m_l_factor1, m_l_short_length);
  ROI_refinement roi_refine(m_cmm, m_nwires[0], m_nwires[1], m_nwires[2],m_r_th_factor,m_r_fake_signal_low_th,m_r_fake_signal_high_th,m_r_pad,m_r_break_roi_loop,m_r_th_peak,m_r_sep_peak,m_r_low_peak_sep_threshold_pre,m_r_max_npeaks,m_r_sigma,m_r_th_percent);//

  
  for (int iplane = 0; iplane != 3; ++iplane){
    // load data into EIGEN matrices ...
    load_data(in, iplane);

    // initial decon ... 
    decon_2D_init(iplane);


    //std::cout << "Form tight ROIs " << i  << std::endl;
    // Form tight ROIs
    if (iplane != 2){ // induction wire planes
      decon_2D_tighterROI(iplane);
      Array::array_xxf r_data_tight(m_r_data.rows(), m_r_data.cols());
      r_data_tight = m_r_data;
      decon_2D_tightROI(iplane);
      roi_form.find_ROI_by_decon_itself(iplane, m_r_data, r_data_tight);
    }else{ // collection wire planes
      decon_2D_tightROI(iplane);
      roi_form.find_ROI_by_decon_itself(iplane, m_r_data);
    }

    
    // Form loose ROIs
    if (iplane != 2){
      decon_2D_looseROI(iplane);

      roi_form.find_ROI_loose(iplane,m_r_data);
      decon_2D_ROI_refine(iplane);
    }

    // Refine ROIs
    roi_refine.load_data(iplane, m_r_data, roi_form);
    roi_refine.refine_data(iplane, roi_form);

    // merge results ...
    decon_2D_hits(iplane);
    roi_refine.apply_roi(iplane, m_r_data);
    save_data(itraces, perframe_traces[iplane], iplane);
    wiener_traces.insert(wiener_traces.end(), perframe_traces[iplane].begin(), perframe_traces[iplane].end());

    decon_2D_charge(iplane);
    roi_refine.apply_roi(iplane, m_r_data);
    save_data(itraces, gauss_traces, iplane);

    m_c_data.resize(0,0); // clear memory
    m_r_data.resize(0,0); // clear memory
  }


  SimpleFrame* sframe = new SimpleFrame(in->ident(), in->time(), itraces, in->tick(), m_cmm);
  sframe->tag_frame("sigproc");

  // this assumes save_data produces itraces in OSP channel order
  std::vector<float> perplane_thresholds[3] = {
    roi_form.get_uplane_rms(),
    roi_form.get_vplane_rms(),
    roi_form.get_wplane_rms()
  };

  IFrame::trace_summary_t threshold;
  for (int iplane=0; iplane<3; ++iplane) {
    for (float val : perplane_thresholds[iplane]) {
      threshold.push_back((double)val);
    }
  }

  sframe->tag_traces("wiener", wiener_traces);
  sframe->tag_traces("threshold", wiener_traces, threshold);
  sframe->tag_traces("gauss", gauss_traces);

  std::cerr << "OmnibusSigProc: produce " << itraces.size() << " traces\n"
	    << "\t" << wiener_traces.size() << " wiener\n"
	    << "\t" << gauss_traces.size() << " gauss\n";

  out = IFrame::pointer(sframe);
  
  return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
