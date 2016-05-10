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

namespace WireCellSigProc {

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
	/// units of the digitization sampling time ("tick").
	/// Defaults are microboone (but you need to give channels).
	void set_gains_shapings(const std::vector<int>& channels,
				double from_gain_mVfC=7.8, double to_gain_mVfC=14.0,
				double from_shaping=2.0, double to_shaping=4.0);


	/// Set the RC+RC time constant in units (multiples) of the
	/// digitization sample time ("tick").  Default is for microboone.
	void set_rcrc_constant(const std::vector<int>& channels, double rcrc=2000.0);

	
	/// Add a constant scaling to the filter for the given
	/// channels between the frequencies.  Frequencies are in
	/// units (multiples) of frequency bin which is
	/// 1.0/(nsamples*tick).  The frequency band is *inclusive* of
	/// both min and max frequency bins.  Note, it's up to caller
	/// to appropriately segment multiple masks across multiple
	/// channels.  For any given channel, last setting wins.
	typedef std::tuple<double, int, int> mask_t;
	typedef std::vector<mask_t> multimask_t;
	void set_filter(const std::vector<int>& channels, const multimask_t& mask);

    private:
	int m_nsamples;
	double m_tick;

	double m_default_baseline, m_default_gain;
	std::vector<double> m_baseline, m_gain;


	typedef std::shared_ptr<filter_t> shared_filter_t;
	typedef std::vector<shared_filter_t> filter_vector_t;
	filter_vector_t m_rcrc, m_config, m_masks;
	shared_filter_t m_default_filter;

	mutable std::unordered_map<int,int> m_ch2ind;
	int chind(int ch) const;

	const IChannelNoiseDatabase::filter_t& get_filter(int channel, const filter_vector_t& fv) const;
    };
}

#endif
