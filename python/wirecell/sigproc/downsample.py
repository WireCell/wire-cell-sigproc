#!/usr/bin/env python
'''
downsampling
'''
import response
import units

import numpy
import matplotlib.pyplot as plt


def checkit(rf):


    gain_mVfC=14
    shaping=2*units.us
    new_nbins = 5000
    target_nbins = 200

    orig = rf.rebin(target_nbins).response
    orig /= numpy.sum(numpy.abs(orig))


    # 1) downsample( ifft( fft(er)*fft(fr) ))
    shaped = rf.shaped(gain_mVfC, shaping, new_nbins)
    onerf = shaped.rebin(target_nbins)
    times = numpy.linspace(*onerf.domainls)
    one = onerf.response
    one /= numpy.sum(numpy.abs(one))
    onepeak = numpy.argmax(one)
    


    # 2) ifft( truncate(fft(er)) * truncate(fft(fr) ))
    rfrebin = rf.rebin(new_nbins)
    elecr = response.electronics(numpy.linspace(*rfrebin.domainls), gain_mVfC, shaping)

    ffter = numpy.fft.fft(elecr)
    fftrf = numpy.fft.fft(rfrebin.response)

    nh = target_nbins/2
    ffter_trunc = numpy.hstack((ffter[:nh], ffter[-nh:]))
    fftrf_trunc = numpy.hstack((fftrf[:nh], fftrf[-nh:]))
    two = numpy.fft.ifft(ffter_trunc * fftrf_trunc)
    two /= numpy.sum(numpy.abs(two))
    twopeak = numpy.argmax(two)
    oneshift = twopeak-onepeak
    one = numpy.hstack((one[-oneshift:], one[:oneshift]))
    
    fig, axes = plt.subplots(3,1, sharex=True)
    axes[0].plot(times, one)
    axes[0].plot(times, two)
    axes[0].plot(times, orig)

    axes[1].plot(times, orig-one)
    axes[1].plot(times, orig-two)
    axes[2].plot(times, one-two)
