#include "TFile.h"
#include "TH1F.h"

#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 3) {
	cerr << "usage: " << argv[0] << " file-with-th1.root name ..." << endl;
	return 1;
    }

    const string fname = argv[1];
    TFile* tfile = TFile::Open(fname.c_str());
    if (!tfile) {
	cerr << "No such file: " << fname << endl;
	return 1;
    }

    for (int argind=2; argind<argc; ++argind) {
	string hname = argv[argind];
	TH1F* hist = (TH1F*)tfile->Get(hname.c_str());
	if (!hist) {
	    cerr << "failed to get histogram " << hname << " from " << fname << endl;
	    return 2;
	}
	const int nbins = hist->GetNbinsX();

	stringstream ss;
	ss << "std::vector<float> " << hname << "{";
	string comma = "";
	for (int hind=0; hind<nbins; ++hind) {
	    ss << comma << hist->GetBinContent(hind+1);
	    comma = ", ";
	}
	ss << "};";
	cout << ss.str();
    }

    return 0;
}
