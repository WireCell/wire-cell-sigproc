#include "WireCellUtil/Waveform.h"

#include "WireCellSigProc/Response.h"

#include "TCanvas.h"
#include "TH1F.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace WireCell;
using namespace WireCellSigProc;

int main()
{
    const double gain_par = 7.8; // mV/fC
    const double shaping_us = 1.0; // microsecond
    Response::ColdElec ce(gain_par, shaping_us);
    const double tick = 0.5;
    const double maxt = 10.0;
    Waveform::signal_t res = ce.generate(tick, 0, maxt);
    Waveform::fourier_t spec = Waveform::fft(res);

    int nticks = res.size();
    TH1F h_wave("response","Cold Electronics Response at 1us shaping", nticks, 0, 10.0);
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

    TCanvas canvas("test_response","Response Functions", 500, 500);
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

    canvas.Print("test_fft.pdf");
}
