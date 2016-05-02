#include "WireCellUtil/Waveform.h"

#include "WireCellSigProc/Response.h"

#include "TCanvas.h"
#include "TH1F.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCellSigProc;

int main()
{
    const double gain_par1 = 7.8; // mV/fC
    const double gain_par2 = 14.0;// mV/fC
    const double shaping_us1 = 1.0; // microsecond
    const double shaping_us2 = 2.0; // microsecond

    Response::ColdElec ce1(gain_par1, shaping_us1);
    Response::ColdElec ce2(gain_par2, shaping_us2);

    const double tick_us=0.05, begin_us=0.0, end_us=10.0;
    const int nticks = (end_us-begin_us)/tick_us;

    // exercise the generator
    Waveform::timeseq_t res1 = ce1.generate(tick_us, begin_us, end_us);
    Waveform::timeseq_t res2 = ce2.generate(tick_us, begin_us, end_us);

    TH1F resp1("resp1","Cold Electronics Response at 1us shaping", nticks, begin_us, end_us);
    TH1F resp2("resp2","Cold Electronics Response at 2us shaping", nticks, begin_us, end_us);
    resp1.SetLineColor(2);
    resp2.SetLineColor(4);
    for (int ind=0; ind<res1.size(); ++ind) {
	resp1.SetBinContent(ind+1, res1[ind]);
	resp2.SetBinContent(ind+1, res2[ind]);
    }

    TCanvas canvas("test_response","Response Functions", 500, 500);
    canvas.SetGridx();
    canvas.SetGridy();
    TH1F* frame = canvas.DrawFrame(0,0,10,15,"Cold Electronics Response Functions (1us,7.8mV/fC and 2us,14.0mV/fC)");
    frame->SetXTitle("Time (microsecond)");
    frame->SetYTitle("Gain (mV/fC)");
    resp1.Draw("same");
    resp2.Draw("same");
    canvas.Print("test_response.pdf");
}
