#include "WireCellUtil/Waveform.h"

#include "WireCellSigProc/Response.h"

#include "TCanvas.h"
#include "TH1F.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace WireCell;
using namespace WireCellSigProc;

void draw_time_freq(TCanvas& canvas, Waveform::timeseq_t& res, Waveform::freqseq_t& spec,
		    const std::string& title,
		    double tick=0.5, double tick0=0.0);
void draw_time_freq(TCanvas& canvas, Waveform::timeseq_t& res, Waveform::freqseq_t& spec,
		    const std::string& title,
		    double tick, double tick0)
{
    int nticks = res.size();
    TH1F h_wave("response",title.c_str(), nticks, tick0, nticks*tick + tick0);

    h_wave.SetXTitle("Time (microsecond)");
    h_wave.SetYTitle("Gain (mV/fC)");

    for (int ind=0; ind<nticks; ++ind) {
	h_wave.SetBinContent(ind+1, res(ind));
    }

    cerr << nticks << " " << spec.size() << endl;

    const double nyquist = 1.0/(2*tick);
    TH1F h_mag("mag","Magnitude of Fourier transform of response", nticks, 0, nyquist);
    h_mag.SetYTitle("Power");
    h_mag.SetXTitle("MHz");

    TH1F h_phi("phi","Phase of Fourier transform of response", nticks, 0, nyquist);
    h_phi.SetYTitle("Power");
    h_phi.SetXTitle("MHz");

    for (int ind=0; ind<nticks; ++ind) {
	auto c = spec(ind);
	h_mag.SetBinContent(ind+1, std::abs(c));
	h_phi.SetBinContent(ind+1, std::arg(c));
    }

    canvas.Clear();
    canvas.Divide(2,1);
    auto pad = canvas.cd(1);
    pad->SetGridx();
    pad->SetGridy();
    h_wave.Draw();

    pad = canvas.cd(2);
    pad->Divide(1,2);

    auto spad = pad->cd(1);
    spad->SetGridx();
    spad->SetGridy();
    h_mag.Draw();

    spad = pad->cd(2);
    spad->SetGridx();
    spad->SetGridy();
    h_phi.Draw();
    canvas.Print("test_fft.pdf",".pdf");

}

int main()
{
    const std::vector<double> gains = {7.8, 14.0}; // mV/fC
    const std::vector<double> shaping_us = {1.0, 2.0}; // microsecond
    const double tick = 0.5;
    const double maxt = 10.0;
    const double tconst_us = 1000.0; // 1ms

    TCanvas canvas("test_fft", "Response Functions", 500, 500);
    canvas.Print("test_fft.pdf[",".pdf");

    for (int ind=0; ind<gains.size(); ++ind) {
	Response::ColdElec ce(gains[ind], shaping_us[ind]);
	Waveform::timeseq_t res = ce.generate(tick, 0, maxt);
	Waveform::freqseq_t spec = Waveform::fft(res);

	draw_time_freq(canvas, res, spec,
		       Form("Cold Electronics Response at %.0fus shaping", shaping_us[ind]), tick);
    }


    {
	Response::SimpleRC rc(tconst_us);
	Waveform::timeseq_t res = rc.generate(tick, 0.0, maxt);
	Waveform::freqseq_t spec = Waveform::fft(res);
	
	draw_time_freq(canvas, res, spec,
		       "RC Response at 1ms time constant", tick);
    }
    {
	Response::SimpleRC rc(tconst_us); 
	Waveform::timeseq_t res = rc.generate(tick, tick, 1000*maxt+tick); // miss the delta
	Waveform::freqseq_t spec = Waveform::fft(res);
	
	draw_time_freq(canvas, res, spec,
		       "RC Response at 1ms time constant (suppress delta)", tick);
    }

    canvas.Print("test_fft.pdf]",".pdf");
}
