#ifndef WIRECELLSIGPROC_OMNICHANNELNOISEDB
#define WIRECELLSIGPROC_OMNICHANNELNOISEDB

#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/WirePlaneId.h"

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Units.h"

#include <vector>
#include <tuple>
#include <unordered_map>
#include <memory>

namespace WireCell {
    namespace SigProc {

	class OmniChannelNoiseDB : public WireCell::IChannelNoiseDatabase , public WireCell::IConfigurable {
	public:

	    /// Create a configurable channel noise DB for digitized
	    /// waveforms with the given size and number of samples.
	    /// Default is for microphone.
	    OmniChannelNoiseDB();
	    virtual ~OmniChannelNoiseDB();

	    // IConfigurable
	    virtual void configure(const WireCell::Configuration& config);
	    virtual WireCell::Configuration default_configuration() const;

	    // IChannelNoiseDatabase
	    virtual int number_samples() const;
	    virtual double sample_time() const;

	    virtual double nominal_baseline(int channel) const;
	    virtual double gain_correction(int channel) const;
            virtual double response_offset(int channel) const;

	    virtual double min_rms_cut(int channel) const;
	    virtual double max_rms_cut(int channel) const;

	    virtual int pad_window_front(int channel) const;
	    virtual int pad_window_back(int channel) const;

	    virtual const filter_t& rcrc(int channel) const;
	    virtual const filter_t& config(int channel) const;
	    virtual const filter_t& noise(int channel) const;
            virtual const filter_t& response(int channel) const;

            // todo:

	    virtual std::vector<channel_group_t> coherent_channels() const {
		return m_channel_groups;
	    }
	    virtual channel_group_t bad_channels() const {
		return m_bad_channels;
	    }


	private:
            double m_tick;
            int m_nsamples;
            IAnodePlane::pointer m_anode;

	    typedef std::shared_ptr<filter_t> shared_filter_t;
	    typedef std::vector<shared_filter_t> filter_vector_t;

	    shared_filter_t m_default_filter;
	    shared_filter_t m_default_response;



            // Embody the "database" entry for one channel. 
            struct ChannelInfo {
                int chid;
    
                // direct scalar values
                double nominal_baseline, gain_correction, response_offset, min_rms_cut, max_rms_cut;
                int pad_window_front, pad_window_back;
    
                // parameters
    
                // frequency space filters
                shared_filter_t rcrc, config, noise, response;
    
                ChannelInfo();
            };

            std::vector<ChannelInfo> m_db;
            //std::unordered_map<int, ChannelInfo*> m_db;

            const ChannelInfo& dbget(int ch, int defch=-1) const {
                // auto it = m_db.find(ch);
                // if (it == m_db.end()) {
                //     it = m_db.find(defch);
                //     return *(it->second);
                // }
                // return *(it->second);
                return m_db.at(ch);
            }

	    std::vector< channel_group_t > m_channel_groups;
	    channel_group_t m_bad_channels;

            // JSON parsing.  Exhausting.
            std::vector<int> parse_channels(const Json::Value& jchannels);
            shared_filter_t make_filter(std::complex<float> defval = std::complex<float>(1,0));
            shared_filter_t default_filter();
            shared_filter_t parse_freqmasks(Json::Value jfm);
            shared_filter_t parse_rcrc(Json::Value jrcrc);
            shared_filter_t parse_reconfig(Json::Value jreconfig);
            shared_filter_t parse_response(Json::Value jreconfig);
            //ChannelInfo* make_ci(int chid, Json::Value jci);
            void update_channels(Json::Value cfg);
            ChannelInfo& get_ci(int chid);


            // Reuse the same filter spectra for matching input parameters.

            // lookup by truncated rcrc value
            std::unordered_map<int, shared_filter_t> m_rcrc_cache;
            // lookup by OR of the four truncated values
            std::unordered_map<int, shared_filter_t> m_reconfig_cache;
            // lookup by explicit waveform id
            std::unordered_map<int, shared_filter_t> m_waveform_cache;
            // lookup by WirePlaneId::ident()
            std::unordered_map<int, shared_filter_t> m_response_cache;

	};
    }

}

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End: