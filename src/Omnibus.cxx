#include "WireCellSigProc/Omnibus.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ExecMon.h"

#include "FrameUtils.h"
using wct::sigproc::dump_frame;

#include <iostream>

WIRECELL_FACTORY(Omnibus, WireCell::SigProc::Omnibus, WireCell::IApplication, WireCell::IConfigurable);

using namespace WireCell;

SigProc::Omnibus::Omnibus()
{
}

SigProc::Omnibus::~Omnibus()
{
}


WireCell::Configuration SigProc::Omnibus::default_configuration() const
{
    Configuration cfg;
    cfg["source"] = "";         // required
    cfg["filters"] = Json::arrayValue;
    cfg["sink"] = "";           // optional
    return cfg;
}

void SigProc::Omnibus::configure(const WireCell::Configuration& cfg)
{
    m_input = Factory::find_tn<IFrameSource>(cfg["source"].asString());
    std::string sink_tn = cfg["sink"].asString();
    if (sink_tn.empty()) {
        std::cerr << "Omnibus has no data sink, will produce no output\n";
        m_output = nullptr;
    }
    else {
        m_output = Factory::find_tn<IFrameSink>(cfg["sink"].asString());
    }
    m_filters.clear();
    auto jffl = cfg["filters"];
    for (auto jff : jffl) {
        auto ff = Factory::find_tn<IFrameFilter>(jff.asString());
        m_filters.push_back(ff);
    }
}


void SigProc::Omnibus::execute()
{
    ExecMon em("omnibus starts");

    IFrame::pointer frame;
    if (!(*m_input)(frame)) {
        std::cerr << "Omnibus: failed to get input frame\n";
        THROW(RuntimeError() << errmsg{"Omnibus: failed to get input frame"});
    }
    if (!frame) {
        std::cerr << "Omnibus: got null frame, forwarding, assuming we have reached EOS\n";
    }
    else {
        std::cerr << "Omnibus: got input frame with " << frame->traces()->size() << " traces\n";
        dump_frame(frame);
    }

    em("sourced frame");

    for (auto ff : m_filters) {
        IFrame::pointer nextframe;
        if (!(*ff)(frame, nextframe)) {
            std::cerr << "Failed to filter frame\n"; // fixme, give more info
            THROW(RuntimeError() << errmsg{"failed to filter frame"});
        }
        if (!nextframe && !frame) {
            continue;           // processing EOS
        }
        if (!nextframe) {
            std::cerr << "Omnibus: filter returned a null frame\n";
            THROW(RuntimeError() << errmsg{"filter returned a null frame"});
        }
        em("filtered frame");

        frame = nextframe;
        nextframe = nullptr;
        if (frame) {
            std::cerr << "Omnibus: filtered frame now with " << frame->traces()->size() << " traces\n";
            dump_frame(frame);
        }
    }

    if (m_output) {
        if (!(*m_output)(frame)) {
            std::cerr << "Omnibus: failed to send output frame\n";
            THROW(RuntimeError() << errmsg{"failed to send output frame"});
        }
    }
    em("sunk frame");

    std::cerr << em.summary() << std::endl;
}
