#include "WireCellSigProc/OmnibusSigProc.h"

#include "WireCellUtil/NamedFactory.h"

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusSigProc::OmnibusSigProc(const std::string anode_tn)
  : m_anode_tn (anode_tn)
{
}

OmnibusSigProc::~OmnibusSigProc()
{
}

void OmnibusSigProc::configure(const WireCell::Configuration& config)
{
  m_anode_tn = get(config, "anode", m_anode_tn);
  m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
  if (!m_anode) {
    THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
  }
}

WireCell::Configuration OmnibusSigProc::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = m_anode_tn; 
  return cfg;
  
}


bool OmnibusSigProc::operator()(const input_pointer& in, output_pointer& out)
{
  return true;
}
