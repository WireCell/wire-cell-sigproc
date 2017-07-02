#ifndef WIRECELLSIGPROC_PEAKFINDING
#define WIRECELLSIGPROC_PEAKFINDING

namespace WireCell{
  namespace SigProc{

    class PeakFinding {
    public:
      PeakFinding(int fMaxPeaks = 200,
		  double sigma = 1, double threshold = 0.05,
		  bool backgroundRemove = false,int deconIterations =3 ,
		  bool markov = true, int averWindow = 3);
      ~PeakFinding();
      
    private:
      int fMaxPeaks;
      double sigma;
      double threshold;
      bool backgroundRemove;
      int deconIterations;
      bool markov;
      int averWindow;
      
      int SearchHighRes(double *source,double *destVector, int ssize,
			double *fPositionX);
    };
  } 
}
#endif
