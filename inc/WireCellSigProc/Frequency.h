/** Things operating in frequency space
 */
#ifndef WIRECELLSIGPROC_FREQUENCY
#define WIRECELLSIGPROC_FREQUENCY

namespace WireCellSigProc {

    namespace Frequency {

	/// A filter applies a frequency sequence 

	class MaskFilter {
	    int _beg, _end;
	public:
	    MaskFilter(double minfreq, double maxfreq, double sampfreq = 2e6);
	    
	};
    }
    

}

#endif

