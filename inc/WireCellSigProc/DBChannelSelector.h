/**
 Database Channel Selector:
 select channels from database
 */

#ifndef WIRECELLSIGPROC_DBCHANNELSELECTOR
#define WIRECELLSIGPROC_DBCHANNELSELECTOR

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSigProc/ChannelSelector.h"

#include <string>
#include <vector>

namespace WireCell {
    namespace SigProc {

        class DBChannelSelector : public ChannelSelector {
        public:

            DBChannelSelector();
            virtual ~DBChannelSelector();

	    /// IConfigurable interface.
	    void configure(const WireCell::Configuration& config);
	    WireCell::Configuration default_configuration() const;
            
        };
    }
}

#endif
