#include "WireCellSigProc/SimpleChannelNoiseDB.h"
#include "WireCellSigProc/Response.h"

using namespace WireCell;
using namespace WireCellSigProc;

SimpleChannelNoiseDB::SimpleChannelNoiseDB(double tick, int nsamples)
    : m_nsamples(0)
    , m_tick(0.0)
    , m_default_baseline(0.0)
    , m_default_gain(1.0)
{
    set_sampling(tick, nsamples);
}
SimpleChannelNoiseDB::~SimpleChannelNoiseDB()
{
}


void SimpleChannelNoiseDB::configure(const WireCell::Configuration& config)
{
}
WireCell::Configuration SimpleChannelNoiseDB::default_configuration() const
{
}

double SimpleChannelNoiseDB::nominal_baseline(int channel) const
{
    const int ind = chind(channel);
    if (0 <= ind && ind < m_baseline.size()) {
	return m_baseline[ind];
    }
    return m_default_baseline;
}

double SimpleChannelNoiseDB::gain_correction(int channel) const
{
    const int ind = chind(channel);
    if (0 <= ind && ind < m_gain.size()) {
	return m_gain[ind];
    }
    return m_default_gain;
}



const IChannelNoiseDatabase::filter_t& SimpleChannelNoiseDB::get_filter(int channel, const filter_vector_t& fv) const
{
    const int ind = chind(channel);
    if (0 <= ind && ind < fv.size()) {
	const shared_filter_t sf = fv[ind];
	const filter_t* filtp = sf.get();
	return *filtp;
    }
    return *(m_default_filter.get());
}

const IChannelNoiseDatabase::filter_t& SimpleChannelNoiseDB::rcrc(int channel) const
{
    return get_filter(channel, m_rcrc);
}

const IChannelNoiseDatabase::filter_t& SimpleChannelNoiseDB::config(int channel) const
{
    return get_filter(channel, m_config);
}

const IChannelNoiseDatabase::filter_t& SimpleChannelNoiseDB::noise(int channel) const
{
    const int ind = chind(channel);
    return get_filter(channel, m_config);
}
	

void SimpleChannelNoiseDB::set_sampling(double tick, int nsamples)
{
    if (m_nsamples == nsamples && tick == m_tick) {
	return;
    }

    m_rcrc.clear();
    m_config.clear();
    m_default_filter = std::make_shared<filter_t>(nsamples);
}
	
// set one thing in a vector at the index, resizing if needed, use def
// to back fill in case of resizing
template<typename T>
void set_one(int ind, T val, std::vector<T>& vec, T def)
{
    if (ind >= vec.size()) {
	vec.resize(ind+1, def);
    }
    vec[ind] = val;
}

void SimpleChannelNoiseDB::set_nominal_baseline(const std::vector<int>& channels, double baseline)
{
    for (auto ch : channels) {
	set_one(chind(ch), baseline, m_baseline, m_default_baseline);
    }
}
void SimpleChannelNoiseDB::set_rcrc_constant(const std::vector<int>& channels, double rcrc)
{
    Response::SimpleRC rcres(rcrc);
    auto signal = rcres.generate(WireCell::Waveform::Domain(0, m_nsamples*m_tick), m_nsamples);
    auto filt = std::make_shared<filter_t>(Waveform::dft(signal));
    for (auto ch : channels) {
	set_one(chind(ch), filt, m_rcrc, m_default_filter);
    }
}

void SimpleChannelNoiseDB::set_gains_shapings(const std::vector<int>& channels,
					      double from_gain, double to_gain,
					      double from_shaping, double to_shaping)
{
    const double gain_ratio = to_gain/from_gain;
    Response::ColdElec from_ce(from_gain, from_shaping);
    Response::ColdElec to_ce(from_gain, from_shaping);
    auto to_sig   =   to_ce.generate(WireCell::Waveform::Domain(0, m_nsamples*m_tick), m_nsamples);
    auto from_sig = from_ce.generate(WireCell::Waveform::Domain(0, m_nsamples*m_tick), m_nsamples);
    auto to_filt   = Waveform::dft(to_sig);
    auto from_filt = Waveform::dft(from_sig);
    Waveform::shrink(to_filt, from_filt); // divide
    auto filt = std::make_shared<filter_t>(to_filt);
    for (auto ch : channels) {
	int ind = chind(ch);
	set_one(ind, filt, m_config, m_default_filter);
	set_one(ind, gain_ratio, m_gain, m_default_gain);
    }
}

void SimpleChannelNoiseDB::set_filter(const std::vector<int>& channels, const multimask_t& masks)
{
    auto filt = std::make_shared<filter_t>(m_nsamples);
    for (auto m : masks) {
	for (int ind=get<1>(m); ind <= get<2>(m); ++ind) {
	    filt->assign(ind, get<0>(m));
	}
    }
    for (auto ch : channels) {
	set_one(chind(ch), filt, m_masks, m_default_filter);
    }
}



int SimpleChannelNoiseDB::chind(int ch) const
{
    auto it = m_ch2ind.find(ch);
    if (it == m_ch2ind.end()) {
	int ind = m_ch2ind.size();
	m_ch2ind[ch] = ind;
	return ind;
    }
    return it->second;
}
