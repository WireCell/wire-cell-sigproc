#include "WireCellSigProc/OmnibusPMTNoiseFilter.h"

#include "WireCellSigProc/Diagnostics.h"

using namespace WireCell;

using namespace WireCell::SigProc;


OmnibusPMTNoiseFilter::OmnibusPMTNoiseFilter()
{
    configure(default_configuration());
}
OmnibusPMTNoiseFilter::~OmnibusPMTNoiseFilter()
{
}

void OmnibusPMTNoiseFilter::configure(const WireCell::Configuration& config)
{
    

    
}
WireCell::Configuration OmnibusPMTNoiseFilter::default_configuration() const
{
    Configuration cfg;
    return cfg;
}

bool OmnibusPMTNoiseFilter::operator()(const input_pointer& in, output_pointer& out)
{
  return true;
}
