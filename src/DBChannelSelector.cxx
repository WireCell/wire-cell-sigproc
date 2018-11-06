#include "WireCellSigProc/DBChannelSelector.h"
#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/FrameTools.h"
#include "WireCellIface/IChannelNoiseDatabase.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(DBChannelSelector, WireCell::SigProc::DBChannelSelector,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SigProc;

DBChannelSelector::DBChannelSelector()
{
}

DBChannelSelector::~DBChannelSelector()
{
}


WireCell::Configuration DBChannelSelector::default_configuration() const
{
    Configuration cfg;

    //must supply
    //database: wclsMiscfgChannelDB, OmniChannelNoiseDB, wclsChannelNoiseDB
    cfg["channelDB"] = "";

    return cfg;
}

void DBChannelSelector::configure(const WireCell::Configuration& cfg)
{
    // channels from database
    auto jcndb = cfg["channelDB"];
    if (!jcndb.empty()) {
        auto db = Factory::find_tn<IChannelNoiseDatabase>(jcndb.asString());
        std::cerr <<"DBChannelSelector: using channel database object: "
            << " \"" << jcndb.asString() << "\"\n";
        std::vector<int> miscfg_channels = db->miscfg_channels(); 
        std::cout << "miscfg channel size: " << miscfg_channels.size() <<"\n";
        for(size_t i=0; i< miscfg_channels.size(); i++) {
        std::cout << "miscfg channel: " << miscfg_channels.at(i) <<"\n";
            m_channels.insert(miscfg_channels.at(i));
        }
    }
    else{
        THROW(ValueError() << errmsg{"DBChannelSelector: no database configured"});
    }
}
