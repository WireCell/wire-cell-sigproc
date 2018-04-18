/** This component applies "compressed sensing" influenced signal
 * processing based on an L1 norm minimzation which fits both a
 * unipolar collection and a bipolar induction response to regions
 * channels in shorted regions known to have a mix.
 */
#ifndef WIRECELLSIGPROC_L1SPFILTER
#define WIRECELLSIGPROC_L1SPFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace SigProc {

        class L1SPFilter : public IFrameFilter,
                           public IConfigurable
        {
        public:
            L1SPFilter();
            virtual ~L1SPFilter();

	    /// IFrameFilter interface.
	    virtual bool operator()(const input_pointer& in, output_pointer& out);

	    /// IConfigurable interface.
	    virtual void configure(const WireCell::Configuration& config);
	    virtual WireCell::Configuration default_configuration() const;

        private:
            Configuration m_cfg;
        };
    }
}

#endif
