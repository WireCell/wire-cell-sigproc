#include "WireCellSigProc/L1SPFilter.h"
#include "WireCellIface/FrameTools.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellIface/SimpleFrame.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

#include <numeric>
#include <iostream>

WIRECELL_FACTORY(L1SPFilter, WireCell::SigProc::L1SPFilter,
                 WireCell::IFrameFilter, WireCell::IConfigurable);

using namespace WireCell;
using namespace WireCell::SigProc;


L1SPFilter::L1SPFilter()
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
    return cfg;
}

void L1SPFilter::configure(const WireCell::Configuration& cfg)
{
    m_cfg = cfg;
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
    
    std::cout << "Xin: " << raw_ROI_th_nsigma << " " << raw_ROI_th_adclimit << " " << overall_time_offset << " " << collect_time_offset << " " << roi_pad << std::endl;
    
    auto adctraces = FrameTools::tagged_traces(in, adctag);
    auto sigtraces = FrameTools::tagged_traces(in, sigtag);


    /// here, use the ADC and signal traces to do L1SP
    ///  put result in out_traces
    ITrace::vector out_traces;

    /// Here is a dummy L1SP filter.  All it does is create new traces
    /// from the input signals.
    for (auto trace : sigtraces) {
        auto newtrace = std::make_shared<SimpleTrace>(trace->channel(), trace->tbin(), trace->charge());
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


