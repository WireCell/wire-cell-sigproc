#include "WireCellUtil/Testing.h"


/// needed to pretend like we are doing WCT internals
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IConfigurable.h"


#include <iostream>

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
    std::string frfname = "ub-10-wnormed.json.bz2";
    if (argc>1) {
        frfname = argv[1];
        cerr << "Using command line field response file: " << frfname << endl;
    }
    else {
        cerr << "Using default field response file: " << frfname << endl;
    }
    

    /// WCT internals, normally user code does not need this
    {
        PluginManager& pm = PluginManager::instance();
        pm.add("WireCellSigProc");
        auto ifrcfg = Factory::lookup<IConfigurable>("FieldResponse");
        auto cfg = ifrcfg->default_configuration();
        cfg["filename"] = frfname;
        ifrcfg->configure(cfg);
    }

    auto ifr = Factory::find<IFieldResponse>("FieldResponse");

    // Get full, "fine-grained" field responses defined at impact
    // positions.
    Response::Schema::FieldResponse fr = ifr->field_response();

    cerr << "FR with " << fr.planes[0].paths.size() << " responses per plane\n";
    Assert(fr.planes[0].paths.size() == 21*6);

    // Make a new data set which is the average FR
    Response::Schema::FieldResponse fravg = Response::wire_region_average(fr);

    cerr << "FR with " << fravg.planes[0].paths.size() << " responses per plane\n";
    /// fixme: why is this is producing 22 responses per plane and not 21?

    // Convert each average FR to a 2D array
    for (int ind=0; ind<3; ++ind) {
        auto arr = Response::as_array(fravg.planes[0]);
        cerr << "FRavg: plane " << ind << ": " << arr.rows() << " X " << arr.cols() << endl;
    }


}
