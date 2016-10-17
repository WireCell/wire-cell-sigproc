#!/usr/bin/env python
'''
Functions related to responses.
'''
import units
import numpy

class ResponseFunction(object):
    '''
    A response function object holds the response wave function and metadata.

    Note: time is assumed to be in Wire Cell system of units (ns).  This is NOT seconds.
    '''
    def __init__(self, plane, region, domainls, response, impact=None):
        plane = plane.lower()
        assert plane in 'uvw'
        self.plane = plane
        self.region = region
        self.domainls = domainls
        self.response = response
        self.times = numpy.linspace(*self.domainls)
        self.impact = impact

    def __call__(self, time):
        return numpy.interp(time, self.times, self.response)

    def rebin(self, nbins):
        newls = (self.times[0], self.times[-1], nbins)
        newtimes = numpy.linspace(*newls)
        newresps = numpy.interp(newtimes, self.times, self.response)
        return ResponseFunction(self.plane, self.region, newls, newresps, self.impact)

    @property
    def nbins(self):
        return self.domainls[2]

    @property
    def asdict(self):
        '''
        Object as a dictionary.
        '''
        return dict(plane=self.plane, region=self.region, domainls=self.domainls,
                    response=self.response.tolist(), impact=self.impact)

    def shaped(self, gain=10.12, shaping=2.0*units.us, nbins=5000):
        '''
        Convolve electronics shaping/peaking response, returning a new ResponseFunction.
        '''
        from scipy.signal import fftconvolve
        import electronics
        newfr = self.rebin(nbins)
        elecr = electronics.response(numpy.linspace(*newfr.domainls), gain, shaping)
        newfr.response = fftconvolve(newfr.response, elecr, "same")
        return newfr

    def __str__(self):
        blah = "<ResponseFunction plane=%s region=%d domainls=%s" % \
               (self.plane, self.region, self.domainls)
        if self.impact:
            blah += " impact=%f" % self.impact
        blah += ">"
        return blah

def group_by(rflist, field):
    '''
    Return a list of lists grouping by like values of the field.  
    '''
    ret = list()
    for thing in sorted(set(getattr(d, field) for d in rflist)):
        bything = [d for d in rflist if getattr(d, field) == thing]
        ret.append(bything)
    return ret



def average(fine):
    '''
    Average fine-grained response functions over multiple impact
    positions in the same plane and wire region.

    Return list of new response.ResponseFunction objects ordered by
    plane, region.
    '''
    ret = list()
    for byplane in group_by(fine, 'plane'):
        for byregion in group_by(byplane, 'region'):
            byregion.sort(key=lambda rf: rf.impact)
            first = byregion[0]
            last = byregion[-1]
            tot = numpy.zeros_like(first.response)
            for rf in byregion:
                tot += rf.response
            tot *= 2.0        # flip responses onto other side of wire
            tot -= first.response # but, don't double
            tot -= last.response # count edge paths
            dat = ResponseFunction(first.plane, first.region, first.domainls, tot)
            ret.append(dat)
        continue
    return ret



def convolve(field, elect):
    '''
    Return the convolution of the two array like responses.
    '''
    from scipy.signal import fftconvolve
    return fftconvolve(field, elect, "same")


def write(rflist, outputfile = "wire-cell-garfield-response.json.bz2"):
    '''
    Write a list of response functions to file.
    '''
    import json
    text = json.dumps([rf.asdict for rf in dat])
    if outputfile.endswith(".json"):
        open(outputfile,'w').write(text)
        return
    if outputfile.endswith(".json.bz2"):
        import bz2
        bz2.BZ2File(outputfile, 'w').write(text)
        return
    if outputfile.endswith(".json.gz"):
        import gzip
        gzip.open(outputfile, "wb").write(text)
        return
    raise ValueError("unknown file format: %s" % outputfile)
# fixme: implement read()
