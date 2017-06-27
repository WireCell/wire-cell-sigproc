#include "WireCellSigProc/FieldResponse.h"
#include "WireCellUtil/Response.h"
#include "WireCellUtil/Testing.h"

#include <iostream>

using namespace WireCell;
using WireCell::SigProc::FieldResponse;
using namespace std;

int main(int argc, char* argv[])
{
    std::string frfname = "ub-10.json.bz2";
    if (argc>1) {
        frfname = argv[1];
        cerr << "Using command line field response file: " << frfname << endl;
    }
    else {
        cerr << "Using default field response file: " << frfname << endl;
    }
    
    // Normally this is done by WCT for us.  And we would do a
    // Factory::find() to get an IFieldResponse.  For brevity for this
    // test we just make it and configure it manually.
    FieldResponse frsource;
    auto cfg = frsource.default_configuration();
    cfg["filename"] = frfname;
    frsource.configure(cfg);

    Response::Schema::FieldResponse fr = frsource.field_response();
    cerr << "FR with " << fr.planes[0].paths.size() << " responses per plane\n";
    Assert(fr.planes[0].paths.size() == 21*6);
    for (int ind=0; ind<3; ++ind) {
        auto arr = Response::as_array(fr.planes[0]);
        cerr << "FR   : plane " << ind << ": " << arr.rows() << " X " << arr.cols() << endl;
    }

    Response::Schema::FieldResponse fravg = Response::wire_region_average(fr);
    cerr << "FR with " << fravg.planes[0].paths.size() << " responses per plane\n";
    for (int ind=0; ind<3; ++ind) {
        auto arr = Response::as_array(fravg.planes[0]);
        cerr << "FRavg: plane " << ind << ": " << arr.rows() << " X " << arr.cols() << endl;
    }


}
