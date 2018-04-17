// This "merges" two frames together by taking all traces from frame 1
// on adding them to frame 2 on a per-tag basis.  If any traces in
// frame 2 of the same channel exist then the rule is applied.

#ifndef WIRECELL_SIGPROC_FRAMEMERGER
#define WIRECELL_SIGPROC_FRAMEMERGER

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameJoiner.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {
    namespace SigProc {
        class FrameMerger : public IFrameJoiner, public IConfigurable {
        public:
            FrameMerger();
            virtual ~FrameMerger();
            
            // IJoinNode
            virtual bool operator()(const input_tuple_type& intup,
                                    output_pointer& out);
            
            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
        private:
            Configuration m_cfg;

        };
    }
}

#endif
