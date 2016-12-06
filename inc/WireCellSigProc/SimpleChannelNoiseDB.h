#ifndef WIRECELLSIGPROC_SIMPLECHANNELNOISEDB
#define WIRECELLSIGPROC_SIMPLECHANNELNOISEDB

#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Units.h"

#include <vector>
#include <tuple>
#include <unordered_map>
#include <memory>

namespace WireCell {
    namespace SigProc {

	class SimpleChannelNoiseDB : public WireCell::IChannelNoiseDatabase , public WireCell::IConfigurable {
	public:

	    /// Create a simple channel noise DB for digitized waveforms
	    /// with the given size and number of samples.  Default is for
	    /// microboone.
	    SimpleChannelNoiseDB(double tick=0.5*units::us, int nsamples=9600);
	    virtual ~SimpleChannelNoiseDB();

	    // IConfigurable
	    virtual void configure(const WireCell::Configuration& config);
	    virtual WireCell::Configuration default_configuration() const;

	    // IChannelNoiseDatabase
	    virtual int number_samples() const { return m_nsamples; }
	    virtual double sample_time() const { return m_tick; }
	    virtual double nominal_baseline(int channel) const;
	    virtual double gain_correction(int channel) const;
	    virtual const filter_t& rcrc(int channel) const;
	    virtual const filter_t& config(int channel) const;
	    virtual const filter_t& noise(int channel) const;
	    virtual std::vector<channel_group_t> coherent_channels() const {
		return m_channel_groups;
	    }
	    virtual channel_group_t bad_channels() const {
		return m_bad_channels;
	    }

	    // concrete helper methods

	    /// Set the size and number of samples of a channel's
	    /// waveform, default is for microboone.
	    ///
	    /// Warning: calling this will reset any settings for
	    /// gains+shaping and rcrc as they depend on knowing the
	    /// sampling.
	    void set_sampling(double tick=0.5*units::us, int nsamples=9600);
	
	    /// Set nominal baseline in units of ADC (eg uB is -2048 for U/V, -400 for W)
	    void set_nominal_baseline(const std::vector<int>& channels, double baseline);

	    /// Set gain/shaping corrections for cnofig_correction.  Gains
	    /// are assumed to be in mV/fC.  Shaping times should be in
	    /// the system of units.  Defaults are microboone (but you
	    /// need to give channels).
	    void set_gains_shapings(const std::vector<int>& channels,
				    double from_gain_mVfC=7.8, double to_gain_mVfC=14.0,
				    double from_shaping=1.0*units::us, double to_shaping=2.0*units::us);


	    /// Set the RC+RC time constant in the system of units for the
	    /// digitization sample time ("tick").  Default is for microboone.
	    void set_rcrc_constant(const std::vector<int>& channels, double rcrc=2000.0);

	
	    /// Set a constant scaling to a band covering the given
	    /// frequency bins (inclusively) for the given channels.
	    /// Frequency bin "i" is from i*f to (i+1)*f where f is
	    /// 1.0/(nsamples*tick).  The largest meaningful frequency bin
	    /// is nsamples/2.  The frequency band is *inclusive* of both
	    /// min and max frequency bins.  Note, it's up to caller to
	    /// appropriately segment multiple masks across multiple
	    /// channels.  For any given channel, last call to this method
	    /// wins.
	    typedef std::tuple<double, int, int> mask_t;
	    typedef std::vector<mask_t> multimask_t;
	    void set_filter(const std::vector<int>& channels, const multimask_t& mask);

	    /// Set the channel groups
	    void set_channel_groups(const std::vector< channel_group_t >& channel_groups) {
		m_channel_groups = channel_groups;
	    }
	    
	    /// Set "bad" channels.
	    void set_bad_channels(const channel_group_t& bc) {
		m_bad_channels = bc;
	    }


	private:
	    double m_tick;
	    int m_nsamples;

	    double m_default_baseline, m_default_gain;
	    std::vector<double> m_baseline, m_gain;


	    typedef std::shared_ptr<filter_t> shared_filter_t;
	    typedef std::vector<shared_filter_t> filter_vector_t;
	    filter_vector_t m_rcrc, m_config, m_masks;
	    shared_filter_t m_default_filter;

	    mutable std::unordered_map<int,int> m_ch2ind;
	    int chind(int ch) const;

	    const IChannelNoiseDatabase::filter_t& get_filter(int channel, const filter_vector_t& fv) const;

	    std::vector< channel_group_t > m_channel_groups;
	    channel_group_t m_bad_channels;
	};
    }

}

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
