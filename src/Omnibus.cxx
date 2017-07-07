#include "WireCellSigProc/Omnibus.h"

#include "WireCellUtil/NamedFactory.h"

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
    cfg["source"] = "";
    cfg["filters"] = Json::arrayValue;
    cfg["sink"] = "";
    return cfg;
}

void SigProc::Omnibus::configure(const WireCell::Configuration& cfg)
{
    m_input = Factory::find_tn<IFrameSource>(cfg["source"].asString());
    m_output = Factory::find_tn<IFrameSink>(cfg["sink"].asString());
    m_filters.clear();
    auto jffl = cfg["filters"];
    for (auto jff : jffl) {
        auto ff = Factory::find_tn<IFrameFilter>(jff.asString());
        m_filters.push_back(ff);
    }
}

void SigProc::Omnibus::execute()
{
    IFrame::pointer frame;
    if (!(*m_input)(frame)) {
        std::cerr << "Failed to get input frame\n";
        return;
    }

    for (auto ff : m_filters) {
        IFrame::pointer nextframe;
        if (!(*ff)(frame, nextframe)) {
            std::cerr << "Failed to filter frame\n"; // fixme, give more info
            return;
        }
        frame = nextframe;
        nextframe = nullptr;
    }

    if (!(*m_output)(frame)) {
        std::cerr << "Failed to get output frame\n";
        return;
    }
}
