#include "WireCellSigProc/FrameMerger.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/FrameTools.h"
#include "WireCellIface/ITrace.h"

#include <string>

WIRECELL_FACTORY(FrameMerger, WireCell::SigProc::FrameMerger,
                 WireCell::IFrameJoiner, WireCell::IConfigurable)

using namespace WireCell;


Configuration SigProc::FrameMerger::default_configuration() const
{
    Configuration cfg;

    // If tags are given then the merge is applied on a per-tag basis
    // and only given tags will be found in the output frame.  If no
    // tags are given then the merge is applied to all traces
    // regardless of any tags and not tags will be found in the output
    // frame.
    cfg["tags"] = Json::arrayValue;

    // Rule says what to do when a trace in frame 2 (tag set) is found
    // of the same channel of a trace in frame 1 (tag set) is found.
    // - replace :: drop frame 2's trace, replace it with frame 1's
    // - include :: simply include 2's trace, keeping also frame 1's as-is
    // - tbd :: more may be added (eg, sum the two)
    cfg["rule"] = "replace";

    return cfg;
}

void SigProc::FrameMerger::configure(const Configuration& cfg)
{
    // fixme: this is maybe uncessessarily slow to stash the cfg.
    m_cfg = cfg;
}

bool SigProc::FrameMerger::operator()(const input_tuple_type& intup,
                                      output_pointer& out)
{
    out = nullptr;

    auto one = std::get<0>(intup);
    auto two = std::get<1>(intup);
    if (!one or !two) {
        std::cerr << "FrameMerger: EOS\n";
        return true;
    }
    std::cerr << "FrameMerger: see frame: "<<one->ident()<<"\n";

    std::vector<ITrace::vector> tracesv1, tracesv2;
    const int ntags = m_cfg["tags"].size();
    if (!ntags) {
        tracesv1.push_back(FrameTools::untagged_traces(one));
        tracesv2.push_back(FrameTools::untagged_traces(two));
    }
    else {
        for (int ind=0; ind<ntags; ++ind) {
            std::string tag = m_cfg["tags"][ind].asString();
            tracesv1.push_back(FrameTools::tagged_traces(one, tag));
            tracesv2.push_back(FrameTools::tagged_traces(two, tag));
        }
    }
        
    ITrace::vector out_traces;
    std::vector<IFrame::trace_list_t> tagged_trace_indices;
    
    auto rule = get<std::string>(m_cfg, "rule"); 
    if (rule == "replace") {

        for (size_t ind=0; ind<tracesv1.size(); ++ind) {
            auto& traces1 = tracesv1[ind];
            auto& traces2 = tracesv2[ind];
            
            std::unordered_map<int, ITrace::pointer> ch2tr;
            for (auto trace : traces2) {
                ch2tr[trace->channel()] = trace;
            }
            for (auto trace : traces1) { // now replace any from frame 1
                ch2tr[trace->channel()] = trace;
            }

            IFrame::trace_list_t tl;
            for (auto chtr : ch2tr) {
                tl.push_back(out_traces.size());
                out_traces.push_back(chtr.second);
            }
            tagged_trace_indices.push_back(tl);
        }            
    }
    if (rule == "include") {
        for (size_t ind=0; ind<tracesv1.size(); ++ind) {
            auto& traces1 = tracesv1[ind];
            auto& traces2 = tracesv2[ind];

            IFrame::trace_list_t tl;
            for (size_t trind=0; trind < traces1.size(); ++trind) {
                tl.push_back(out_traces.size());
                out_traces.push_back(traces1[trind]);
            }
            for (size_t trind=0; trind < traces2.size(); ++trind) {
                tl.push_back(out_traces.size());
                out_traces.push_back(traces2[trind]);
            }
            tagged_trace_indices.push_back(tl);
        }
    }

    auto sf = new SimpleFrame(two->ident(), two->time(), out_traces, two->tick());
    if (ntags) {
        for (int ind=0; ind<ntags; ++ind) {
            std::string tag = m_cfg["tags"][ind].asString();
            sf->tag_traces(tag, tagged_trace_indices[ind]);
        }
    }
    out = IFrame::pointer(sf);
    return true;
}


SigProc::FrameMerger::FrameMerger()
{
}

SigProc::FrameMerger::~FrameMerger()
{
}
