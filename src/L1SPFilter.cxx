#include "WireCellSigProc/L1SPFilter.h"
#include "WireCellIface/FrameTools.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellIface/SimpleFrame.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/IFieldResponse.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

#include <numeric>
#include <iostream>

WIRECELL_FACTORY(L1SPFilter, WireCell::SigProc::L1SPFilter,
                 WireCell::IFrameFilter, WireCell::IConfigurable);

using namespace WireCell;
using namespace WireCell::SigProc;


L1SPFilter::L1SPFilter(double gain, 
		       double shaping_time,
		       double inter_gain , 
		       double ADC_mV)
  : m_gain(gain)
  , m_shaping_time(shaping_time)
  , m_inter_gain(inter_gain)
  , m_ADC_mV(ADC_mV)
{
}

L1SPFilter::~L1SPFilter()
{
}

WireCell::Configuration L1SPFilter::default_configuration() const
{
    Configuration cfg;

    /// An array holding a waveform to use as the "smearing" filter.
    cfg["filter"] = Json::arrayValue;

    /// The tag identifying traces which represent "raw" (not
    /// deconvolved) ADC values.
    cfg["adctag"] = "raw";

    /// The tag identifying traces which represent "signal" processed
    /// (deconvolved) waveforms.
    cfg["sigtag"] = "gauss";

    /// The tag to place on the output waveforms
    cfg["outtag"] = "l1sp";

    // 4 sigma for raw waveform ROI identification
    cfg["raw_ROI_th_nsigma"] = 4;
    // 10 ADC for upper limit on ADC ... 
    cfg["raw_ROI_th_adclimit"] = 10;
    // global offset 
    cfg["overall_time_offset"] = 0;
    // need 3 us offset for collection plane relative to the induction plane ...
    cfg["collect_time_offset"] = 3.0;

    // ROI padding ticks ...
    cfg["roi_pad"] = 20;

    // L1 fit parameters ...
    cfg["adc_l1_threshold"] = 6;
    cfg["adc_sum_threshold"] = 160;
    cfg["adc_sum_rescaling"] = 90.;
    cfg["adc_sum_rescaling_limit"] = 50.;
    cfg["l1_seg_length"] = 120;
    cfg["l1_scaling_factor"] = 500;
    cfg["l1_lambda"] = 5;
    cfg["l1_epsilon"] = 0.05;
    cfg["l1_niteration"] = 100000;
    cfg["l1_decon_limit"] = 50; // 50 electrons

    cfg["l1_resp_scale"] = 0.5;
    cfg["l1_col_scale"] = 1.15;
    cfg["l1_ind_scale"] = 0.5;

    
    cfg["gain"] = m_gain;
    cfg["shaping_time"] = m_shaping_time;
    cfg["inter_gain"] = m_inter_gain;
    cfg["ADC_mV"] = m_ADC_mV;
    
    return cfg;
}

void L1SPFilter::configure(const WireCell::Configuration& cfg)
{
    m_cfg = cfg;

    m_gain = get(cfg,"gain",m_gain);
    m_shaping_time = get(cfg,"shaping_time",m_shaping_time);
    m_inter_gain = get(cfg,"inter_gain", m_inter_gain);
    m_ADC_mV = get(cfg,"ADC_mV", m_ADC_mV);

    
}

bool L1SPFilter::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        return true;            // eos
    }

    std::string adctag = get<std::string>(m_cfg, "adctag");
    std::string sigtag = get<std::string>(m_cfg, "sigtag");
    std::string outtag = get<std::string>(m_cfg, "outtag");
    std::vector<float> signal = get< std::vector<float> >(m_cfg, "filter");

    double raw_ROI_th_nsigma = get(m_cfg,"raw_ROI_th_nsigma",raw_ROI_th_nsigma);
    double raw_ROI_th_adclimit = get(m_cfg,"raw_ROI_th_adclimit",raw_ROI_th_adclimit);
    double overall_time_offset = get(m_cfg,"overall_time_offset",overall_time_offset);
    double collect_time_offset = get(m_cfg,"collect_time_offset",collect_time_offset);
    int roi_pad = 0;
    roi_pad = get(m_cfg,"roi_pad",roi_pad);

    double adc_l1_threshold = get(m_cfg,"adc_l1_threshold",adc_l1_threshold);
    double adc_sum_threshold= get(m_cfg,"adc_sum_threshold",adc_sum_threshold);
    double adc_sum_rescaling= get(m_cfg,"adc_sum_rescaling",adc_sum_rescaling);
    double adc_sum_rescaling_limit= get(m_cfg,"adc_sum_rescaling_limit",adc_sum_rescaling_limit);
    double l1_seg_length= get(m_cfg,"l1_seg_length",l1_seg_length);
    double l1_scaling_factor= get(m_cfg,"l1_scaling_factor",l1_scaling_factor);
    double l1_lambda= get(m_cfg,"l1_lambda",l1_lambda);
    double l1_epsilon= get(m_cfg,"l1_epsilon",l1_epsilon);
    double l1_niteration= get(m_cfg,"l1_niteration",l1_niteration);
    double l1_decon_limit= get(m_cfg,"l1_decon_limit",l1_decon_limit);
    
    double l1_resp_scale = get(m_cfg,"l1_resp_scale",l1_resp_scale);
    double l1_col_scale = get(m_cfg,"l1_col_scale",l1_col_scale);
    double l1_ind_scale = get(m_cfg,"l1_ind_scale",l1_ind_scale);
    
    std::cout << "Xin: " << raw_ROI_th_nsigma << " " << raw_ROI_th_adclimit << " " << overall_time_offset << " " << collect_time_offset << " " << roi_pad << " " << adc_l1_threshold << " " << adc_sum_threshold << " " << adc_sum_rescaling << " " << adc_sum_rescaling_limit << " " << l1_seg_length << " " << l1_scaling_factor << " " << l1_lambda << " " << l1_epsilon << " " << l1_niteration << " " << l1_decon_limit << " " << l1_resp_scale << " " << l1_col_scale << " " << l1_ind_scale << std::endl;

    
    // // get field response ... 
    // auto ifr = Factory::find<IFieldResponse>("FieldResponse");
    // Response::Schema::FieldResponse fr = ifr->field_response();
    // // Make a new data set which is the average FR
    // Response::Schema::FieldResponse fravg = Response::wire_region_average(fr);

    // get electronics response
    WireCell::Waveform::compseq_t elec;
    Response::ColdElec ce(m_gain, m_shaping_time);
    // auto ewave = ce.generate(tbins);
    // Waveform::scale(ewave, m_inter_gain * m_ADC_mV);

    // get the overall response function ...
    // how to use it ... 
    
    auto adctraces = FrameTools::tagged_traces(in, adctag);
    auto sigtraces = FrameTools::tagged_traces(in, sigtag);


    /// here, use the ADC and signal traces to do L1SP
    ///  put result in out_traces
    ITrace::vector out_traces;
    
    std::map<int,std::set<int>> init_map;
    // do ROI from the decon signal
    for (auto trace : sigtraces) {
      int ch = trace->channel();
      int tbin = trace->tbin();
      auto const& charges = trace->charge();
      const int ntbins = charges.size();
      std::set<int> time_ticks;

      for (int qi = 0; qi < ntbins; qi++){
	if (charges[qi]!=0){
	  time_ticks.insert(tbin+qi);
	}
      }
      
      init_map[ch] = time_ticks;
      // if (time_ticks.size()>0){
      // 	std::cout << ch << " " << time_ticks.size() << std::endl;
      // }
    }
    
    // do ROI from the raw signal
    int ntot_ticks=0;
    
    for (auto trace : adctraces) {
      int ch = trace->channel();
      int tbin = trace->tbin();
      auto const& charges = trace->charge();
      const int ntbins = charges.size();
      std::set<int>& time_ticks = init_map[ch];

      if (ntot_ticks < ntbins)
	ntot_ticks = ntbins;
      
      double mean = Waveform::percentile(charges,0.5);
      double mean_p1sig = Waveform::percentile(charges,0.5+0.34);
      double mean_n1sig = Waveform::percentile(charges,0.5-0.34);
      double cut = raw_ROI_th_nsigma * sqrt((pow(mean_p1sig-mean,2)+pow(mean_n1sig-mean,2))/2.);
      if (cut < raw_ROI_th_adclimit) cut = raw_ROI_th_adclimit;

      for (int qi = 0; qi < ntbins; qi++){
	if (fabs(charges[qi])>cut){
	  time_ticks.insert(tbin+qi);
	}
      }
      // if (time_ticks.size()>0){
      // 	std::cout << ch << " " << time_ticks.size() << std::endl;
      // }
    }

    
    // create ROIs ... 
    std::map<int, std::vector<std::pair<int,int>>> map_ch_rois;
    
    for (auto it = init_map.begin(); it!=init_map.end(); it++){
      int wire_index = it->first;
      std::set<int>& time_slices_set = it->second;
      if (time_slices_set.size()==0) continue;
      std::vector<int> time_slices;
      std::copy(time_slices_set.begin(), time_slices_set.end(), std::back_inserter(time_slices));
      
      std::vector<std::pair<int,int>> rois;
      std::vector<std::pair<int,int>> rois_save;
      
      rois.push_back(std::make_pair(time_slices.front(),time_slices.front()));
      for (size_t i=1; i<time_slices.size();i++){
	if (time_slices.at(i) - rois.back().second <= roi_pad*2){
	  rois.back().second = time_slices.at(i);
	}else{
	  rois.push_back(std::make_pair(time_slices.at(i),time_slices.at(i)));
	}
      }
      
      // extend the rois to both side according to the bin content
      for (auto it = rois.begin(); it!= rois.end();  it++){
	int start_bin = it->first;
	int end_bin = it->second;
	start_bin = start_bin - roi_pad;
	end_bin = end_bin + roi_pad;
	if (start_bin <0) start_bin = 0;
	if (end_bin>ntot_ticks-1) end_bin = ntot_ticks-1;
	it->first = start_bin;
	it->second = end_bin;
      }
      
      for (auto it = rois.begin(); it!= rois.end();  it++){
	if (rois_save.size()==0){
	  rois_save.push_back(*it);
	}else if (it->first <= rois_save.back().second){
	  rois_save.back().second = it->second;
	}else{
	  rois_save.push_back(*it);
	}
      }

      if (rois_save.size()>0)
	map_ch_rois[wire_index] = rois_save;
      // for (auto it = rois_save.begin(); it!=rois_save.end(); it++){
      // std::cout << wire_index << " " << it->first << " " << it->second +1 << std::endl;
      // }
    }
    
    
    // prepare for the output signal ...
    
    for (auto trace : sigtraces) {
      auto newtrace = std::make_shared<SimpleTrace>(trace->channel(), trace->tbin(), trace->charge());
      // How to access the sigtraces together ???
      if (map_ch_rois.find(trace->channel()) != map_ch_rois.end()){
	std::vector<std::pair<int,int>>& rois_save = map_ch_rois[trace->channel()];
	for (auto it = rois_save.begin(); it!=rois_save.end(); it++){
	  for (int time_tick = it->first; time_tick!=it->second+1; time_tick++){
	    // temporary hack to reset the data ... 
	    newtrace->charge().at(time_tick-trace->tbin())=0;
	  }
	}
      }

      
      // std::cout << trace->channel() << std::endl;
      out_traces.push_back(newtrace);
    }



    std::cerr << "L1SPFilter: frame: " << in->ident() << " "
              << adctag << "[" << adctraces.size() << "] + "
              << sigtag << "[" << sigtraces.size() << "] --> "
              << outtag << "[" << out_traces.size() << "]\n";


    
    // Finally, we save the traces to an output frame with tags.

    IFrame::trace_list_t tl(out_traces.size());
    std::iota(tl.begin(), tl.end(), 0);
    
    auto sf = new SimpleFrame(in->ident(), in->time(), out_traces, in->tick());
    sf->tag_traces(outtag, tl);
    out = IFrame::pointer(sf);
    
    return true;
}


