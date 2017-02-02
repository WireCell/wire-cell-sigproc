#include "WireCellSigProc/SimpleChannelNoiseDB.h"
#include "WireCellUtil/Response.h"
//#include <iostream> //debug

using namespace WireCell;
using namespace WireCell::SigProc;

SimpleChannelNoiseDB::SimpleChannelNoiseDB(double tick, int nsamples)
    : m_tick(-1)
    , m_nsamples(-1)
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
    //std::cerr << "ch=" << channel << " ind=" << ind << " " << fv.size() << std::endl;
    if (0 <= ind && ind < fv.size()) {
	const shared_filter_t sf = fv[ind];
	if (sf == nullptr) {
	    return *(m_default_filter.get());
	}

	const filter_t* filtp = sf.get();
	return *filtp;
    }
    const filter_t* filtp = m_default_filter.get();
    //std::cerr << "Filter: "<< (void*)filtp << std::endl;
    return *filtp;
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
    return get_filter(channel, m_masks);
}
	
const IChannelNoiseDatabase::filter_t& SimpleChannelNoiseDB::response(int channel) const
{
    return get_filter(channel, m_response);
}


void SimpleChannelNoiseDB::set_sampling(double tick, int nsamples)
{

    if (m_nsamples == nsamples && tick == m_tick) {
	//std::cerr << "Sampling unchanged: " << nsamples << " @ " << m_tick/units::ms << " ms" << std::endl;
	return;
    }
    m_nsamples = nsamples;
    m_tick = tick;

    m_rcrc.clear();
    m_config.clear();
    m_response.clear();

    Waveform::compseq_t spectrum;
    spectrum.resize(nsamples,std::complex<float>(1,0));
    m_default_filter = std::make_shared<filter_t>(spectrum);
    Waveform::compseq_t empty;
    m_default_response = std::make_shared<filter_t>(empty);
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
    Response::SimpleRC rcres(rcrc,m_tick);
    auto signal = rcres.generate(WireCell::Waveform::Domain(0, m_nsamples*m_tick), m_nsamples);
    
    
    Waveform::compseq_t spectrum = Waveform::dft(signal);
    
    //std::cout << rcrc << " " << m_tick << " " << m_nsamples << " " << signal.front() << " " << signal.at(1) << " " << signal.at(2) << std::endl;

    // get the square of it because there are two RC filters
    Waveform::compseq_t spectrum2 = spectrum;
    Waveform::scale(spectrum2,spectrum);
    // for (auto it : spectrum){
    //   float real_part = it.real();
    //   float imag_part = it.imag();
    //   std::complex<float> temp(real_part*real_part-imag_part*imag_part,2*real_part*imag_part);
    //   spectrum2.push_back(temp);
    // }

    auto filt = std::make_shared<filter_t>(spectrum2);
    
    for (auto ch : channels) {
	set_one(chind(ch), filt, m_rcrc, m_default_filter);
    }
}

void SimpleChannelNoiseDB::set_response(const std::vector<int>& channels, const filter_t& spectrum)
{
    auto filt = std::make_shared<filter_t>(spectrum);
    for (auto ch : channels) {
        set_one(chind(ch), filt, m_response, m_default_response);
    }

}

void SimpleChannelNoiseDB::set_gains_shapings(const std::vector<int>& channels,
					      double from_gain, double to_gain,
					      double from_shaping, double to_shaping)
{
    const double gain_ratio = to_gain/from_gain;
    Response::ColdElec from_ce(from_gain, from_shaping);
    Response::ColdElec to_ce(to_gain, to_shaping);
    auto to_sig   =   to_ce.generate(WireCell::Waveform::Domain(0, m_nsamples*m_tick), m_nsamples);
    auto from_sig = from_ce.generate(WireCell::Waveform::Domain(0, m_nsamples*m_tick), m_nsamples);
    
    //std::cout << to_gain << " " << from_gain << " " << to_shaping << " " << from_shaping << " " << to_sig.at(1) << " " << from_sig.at(1) << std::endl;
    
    auto to_filt   = Waveform::dft(to_sig);
    auto from_filt = Waveform::dft(from_sig);
    Waveform::shrink(to_filt, from_filt); // divide
    auto filt = std::make_shared<filter_t>(to_filt);
    
    //std::cout << to_filt.at(1) << " " << to_filt.at(2) << std::endl;
    for (auto ch : channels) {
	int ind = chind(ch);
	//	std::cout << ch << " " << ind << std::endl;
	set_one(ind, filt, m_config, m_default_filter);
	set_one(ind, gain_ratio, m_gain, m_default_gain);
    }
}

void SimpleChannelNoiseDB::set_filter(const std::vector<int>& channels, const multimask_t& masks)
{
    Waveform::compseq_t spectrum;
    spectrum.resize(m_nsamples,std::complex<float>(1,0));

    for (auto m : masks) {
	for (int ind=get<1>(m); ind <= get<2>(m); ++ind) {
	    //filt->assign(ind, get<0>(m));
	    spectrum.at(ind) = get<0>(m);
	}
	//std::cout << "Xin: " << get<1>(m) << " " << get<2>(m) << " " << get<0>(m) << " " << spectrum.at(0) << " " << spectrum.at(169) << " " << spectrum.at(170)<< std::endl;
    }
    auto filt = std::make_shared<filter_t>(spectrum);
  

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

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
