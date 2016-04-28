/** Things operating in frequency space
 */
#ifndef WIRECELLSIGPROC_FREQUENCY
#define WIRECELLSIGPROC_FREQUENCY

namespace WireCellSigProc {

    namespace Frequency {

	



	

	class MaskFilter {
	    int beg, end;
	public:
	    MaskFilter(double minfreq, double maxfreq, double sampfreq = 2e6);
	    
	};
    }
    

}

#endif

