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

    cfg["type"] = "misconfigured";

    //must supply
    //database: wclsMiscfgChannelDB, OmniChannelNoiseDB, wclsChannelNoiseDB
    cfg["channelDB"] = "";

    return cfg;
}

void DBChannelSelector::configure(const WireCell::Configuration& cfg)
{
    std::string type = get<std::string>(cfg, "type", "misconfigured");

    // channels from database
    auto jcndb = cfg["channelDB"];
    if (!jcndb.empty()) {
        auto db = Factory::find_tn<IChannelNoiseDatabase>(jcndb.asString());
        std::cerr <<"DBChannelSelector: using channel database object: "
            << " \"" << jcndb.asString() << "\" type: \"" << type << "\"\n";
        if( type == "bad" ) ChannelSelector::set_channels(db->bad_channels());
        if( type == "misconfigured" ) ChannelSelector::set_channels(db->miscfg_channels()); 
        //std::vector<int> miscfg_channels = db->miscfg_channels(); 
        //std::cout << "miscfg channel size: " << miscfg_channels.size() <<"\n";
        //for(size_t i=0; i< miscfg_channels.size(); i++) {
        //std::cout << "miscfg channel: " << miscfg_channels.at(i) <<"\n";
        //}
    }
    else{
        THROW(ValueError() << errmsg{"DBChannelSelector: no database configured"});
    }
}
