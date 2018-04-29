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
    
    std::cout << "Xin: " << raw_ROI_th_nsigma << " " << raw_ROI_th_adclimit << " " << overall_time_offset << " " << collect_time_offset << " " << roi_pad << " " << adc_l1_threshold << " " << adc_sum_threshold << " " << adc_sum_rescaling << " " << adc_sum_rescaling_limit << " " << l1_seg_length << " " << l1_scaling_factor << " " << l1_lambda << " " << l1_epsilon << " " << l1_niteration << " " << l1_decon_limit << std::endl;

    
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
    
    auto adctraces = FrameTools::tagged_traces(in, adctag);
    auto sigtraces = FrameTools::tagged_traces(in, sigtag);


    /// here, use the ADC and signal traces to do L1SP
    ///  put result in out_traces
    ITrace::vector out_traces;

    // do ROI from the raw signal
    
    // do ROI from the decon signal

    // merge ROIs ... 
    
    // prepare for the output signal ...
    
    for (auto trace : sigtraces) {
      auto newtrace = std::make_shared<SimpleTrace>(trace->channel(), trace->tbin(), trace->charge());
      // How to access the sigtraces together ???
      
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


