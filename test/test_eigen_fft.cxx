#include <Eigen/Dense>

#include <unsupported/Eigen/FFT>

#include <iostream>
using namespace std;

int main()
{
#ifdef MISSING_FFTW_SINGLE_PRECISION
    typedef double TYPE;
#else
    typedef float TYPE;
#endif

    Eigen::FFT<TYPE> fft;
    std::vector<TYPE> timevec = {1,2,3,2,1};
    std::vector<std::complex<TYPE> > freqvec;

    fft.fwd( freqvec,timevec);
    // manipulate freqvec
    fft.inv( timevec,freqvec);

    for (auto x: timevec) {
	cerr << x << " ";
    }
    cerr << endl;
    for (auto x: freqvec) {
	cerr << x << " ";
    }
    cerr << endl;

    // stored "plans" get destroyed with fft destructor

    return 0;
}
