#!/usr/bin/env python
'''
Functions related to responses.
'''
import units
import numpy

def electronics_no_gain_scale(time, gain, shaping=2.0*units.us):
    '''
    This version takes gain parameter already scaled such that the
    gain actually desired is obtained.  
    '''
    domain=(0, 10*units.us)
    if time <= domain[0] or time >= domain[1]:
        return 0.0

    time = time/units.us
    st = shaping/units.us
        
    from math import sin, cos, exp
    ret = 4.31054*exp(-2.94809*time/st)*gain-2.6202*exp(-2.82833*time/st)*cos(1.19361*time/st)*gain \
          -2.6202*exp(-2.82833*time/st)*cos(1.19361*time/st)*cos(2.38722*time/st)*gain \
          +0.464924*exp(-2.40318*time/st)*cos(2.5928*time/st)*gain \
          +0.464924*exp(-2.40318*time/st)*cos(2.5928*time/st)*cos(5.18561*time/st)*gain \
          +0.762456*exp(-2.82833*time/st)*sin(1.19361*time/st)*gain \
          -0.762456*exp(-2.82833*time/st)*cos(2.38722*time/st)*sin(1.19361*time/st)*gain \
          +0.762456*exp(-2.82833*time/st)*cos(1.19361*time/st)*sin(2.38722*time/st)*gain \
          -2.6202*exp(-2.82833*time/st)*sin(1.19361*time/st)*sin(2.38722*time/st)*gain  \
          -0.327684*exp(-2.40318*time/st)*sin(2.5928*time/st)*gain +  \
          +0.327684*exp(-2.40318*time/st)*cos(5.18561*time/st)*sin(2.5928*time/st)*gain \
          -0.327684*exp(-2.40318*time/st)*cos(2.5928*time/st)*sin(5.18561*time/st)*gain \
          +0.464924*exp(-2.40318*time/st)*sin(2.5928*time/st)*sin(5.18561*time/st)*gain
    
    return ret

def electronics(time, gain_mVfC=14, shaping=2.0*units.us):
    '''
    Electronics response function.

        - gain :: the gain such that the peak is the desired mV/fC

        - shaping :: the shaping time in Wire Cell system of units

        - domain :: outside this pair, the response is identically zero
    '''
    # see wirecell.sigproc.plots.electronics() for these magic numbers.
    if shaping <= 0.5*units.us:
        gain = gain_mVfC*10.146826
    elif shaping <= 1.0*units.us:
        gain = gain_mVfC*10.146828
    elif shaping <= 2.0*units.us:
        gain = gain_mVfC*10.122374
    else:
        gain = gain_mVfC*10.120179
    return electronics_no_gain_scale(time, gain, shaping)

electronics = numpy.vectorize(electronics)


class ResponseFunction(object):
    '''
    A response function object holds the response wave function and metadata.

    Note: time is assumed to be in Wire Cell system of units (ns).  This is NOT seconds.
    '''
    def __init__(self, plane, region, pos, domainls, response, impact=None):
        plane = plane.lower()
        assert plane in 'uvw'
        self.plane = plane
        self.region = region
        self.pos = tuple(pos)
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
        return self.dup(domainls=newls, response=newresps)

    def dup(self, **kwds):
        '''
        Return a new ResponseFunction which is a copy of this one and
        with any values in kwds overriding.
        '''
        return ResponseFunction(**dict(self.asdict, **kwds))

    @property
    def nbins(self):
        return self.domainls[2]

    @property
    def asdict(self):
        '''
        Object as a dictionary.
        '''
        return dict(plane=self.plane, region=self.region, pos=self.pos,
                    domainls=self.domainls, response=self.response.tolist(), impact=self.impact)

    def shaped(self, gain_mVfC=14, shaping=2.0*units.us, nbins=5000):
        '''
        Convolve electronics shaping/peaking response, returning a new ResponseFunction.
        '''
        from scipy.signal import fftconvolve
        import electronics
        newfr = self.rebin(nbins)
        elecr = electronics(numpy.linspace(*newfr.domainls), gain_mVfC, shaping)
        newfr.response = fftconvolve(newfr.response, elecr, "same")
        return newfr

    def __str__(self):
        blah = "<ResponseFunction plane=%s region=%d domainls=%s pos=%s" % \
               (self.plane, self.region, self.domainls, self.pos)
        if self.impact is not None:
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

def by_region(rflist, region=0):
    ret = [rf for rf in rflist if rf.region == region]
    ret.sort(key=lambda x: x.plane)
    return ret


def average(fine):
    '''
    Average fine-grained response functions over multiple impact
    positions in the same plane and wire region.

    Return list of new response.ResponseFunction objects ordered by
    plane, region.
    '''
    ret = list()
    for inplane in group_by(fine, 'plane'):
        byregion = group_by(inplane, 'region')
        noigeryb = list(byregion)
        noigeryb.reverse()

        center = max([rflist[0].region for rflist in byregion])

        for regp,regm in zip(byregion[center:], noigeryb[center:]):
            regp.sort(key=lambda x: x.impact)
            regm.sort(key=lambda x: x.impact)
            tot = numpy.zeros_like(regp[0].response)
            for one in regp + regm:
                tot += one.response
            tot *= 2.0
            tot -= regp[0].response + regm[0].response # don't double count impact=0
            tot -= regp[-1].response + regm[-1].response # share region boundary path
            # warning: this makes detailed assumptions about where the impact points are!!11!
            dat = regp[0].dup(response=tot, impact=None)
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
    text = json.dumps([rf.asdict for rf in rflist])
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
