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
            L1SPFilter(double gain = 14.0 * units::mV/units::fC, 
		       double shaping = 2.0 * units::microsecond,
		       double postgain = 1.2, 
		       double ADC_mV = 4096/(2000.*units::mV),
		       double fine_time_offset = 0.0 * units::microsecond,
		       double coarse_time_offset = -8.0 * units::microsecond);
            virtual ~L1SPFilter();

	    /// IFrameFilter interface.
	    virtual bool operator()(const input_pointer& in, output_pointer& out);

	    /// IConfigurable interface.
	    virtual void configure(const WireCell::Configuration& config);
	    virtual WireCell::Configuration default_configuration() const;

        private:
            Configuration m_cfg;

	    double m_gain;
	    double m_shaping;
	    double m_postgain;
	    double m_ADC_mV;
	    double m_fine_time_offset;
	    double m_coarse_time_offset;
	    
        };
    }
}

#endif
