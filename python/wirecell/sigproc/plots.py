#!/usr/bin/env python
import units
import response

import numpy
import matplotlib.pyplot as plt


    

def response_by_wire_region(rflist):
    '''
    Plot response functions as 1D graphs.
    '''
    one = rflist[0]
    byplane = response.group_by(rflist, 'plane')

    nwires = map(len, byplane)
    print "%d planes, nwires: %s" % (len(nwires), str(nwires))

    nwireshalf = nwires[0]//2
    wire0 = nwires[0] - nwireshalf - 1

    central_collection = byplane[2][wire0]
    central_sum_field = sum(central_collection.response)
    central_shaped = central_collection.shaped()
    central_sum_shape = sum(central_shaped.response)
    
    print '#regions: %d' % (nwireshalf+1,)
    fig, axes = plt.subplots(nwireshalf+1, 2, sharex=True)
    print axes

    for wire_offset in range(nwireshalf+1):
        axf = axes[wire_offset][0]
        axf.set_title('Wire region %d (field)' % wire_offset);
        axs = axes[wire_offset][1]
        axs.set_title('Wire region %d (shaped)' % wire_offset);

        iwirep = wire0 + wire_offset
        iwirem = wire0 + wire_offset

        for iplane in range(3):
            plane = byplane[iplane]
            fieldp = plane[iwirep]
            fieldm = plane[iwirem]
            shapep = fieldp.shaped() # fixme, just taking
            shapem = fieldm.shaped() # default shaping here....

            field = fieldp.response
            shape = shapep.response
            if wire_offset:     # not central
                field = 0.5*(fieldp.response + fieldm.response)
                shape = 0.5*(shapep.response + shapem.response)
            field /= central_sum_field
            shape /= central_sum_shape
            
            ftime = 1.0e6*numpy.linspace(*fieldp.domainls)
            stime = 1.0e6*numpy.linspace(*shapep.domainls)

            axf.plot(ftime, field)
            axs.plot(stime, shape)


def response_averages_colz(avgtriple, time):
    '''
    Plot averages as 2D colz type plot
    '''
    use_imshow = False
    mintbin=700
    maxtbin=850
    nwires = avgtriple[0].shape[0]
    maxwires = nwires//2    
    minwires = -maxwires
    mintime = time[mintbin]
    maxtime = time[maxtbin-1]
    ntime = maxtbin-mintbin
    deltatime = (maxtime-mintime)/ntime

    x,y = numpy.meshgrid(numpy.linspace(mintime, maxtime, ntime),
                          numpy.linspace(minwires, maxwires, nwires))
    x *= 1.0e6                  # put into us

    print x.shape, mintbin, maxtbin, mintime, maxtime, nwires, minwires, maxwires

    fig = plt.figure()
    cmap = 'seismic'

    toplot=list()
    for iplane in range(3):
        avg = avgtriple[iplane]
        main = avg[:,mintbin:maxtbin]
        edge = avg[:,maxtbin:]
        ped = numpy.sum(edge) / (edge.shape[0] * edge.shape[1])
        toplot.append(main - ped)

    maxpix = max(abs(numpy.min(avgtriple)), numpy.max(avgtriple))
    clim = (-maxpix/2.0, maxpix/2.0)

    ims = list()
    axes = list()

    for iplane in range(3):
        ax = fig.add_subplot(3,1,iplane+1) # two rows, one column, first plot
        if use_imshow:
            im = plt.imshow(toplot[iplane], cmap=cmap, clim=clim,
                            extent=[mintime, maxtime, minwires, maxwires], aspect='auto')
        else:
            im = plt.pcolormesh(x,y, toplot[iplane], cmap=cmap, vmin=clim[0], vmax=clim[1])
        ims.append(im)
        axes.append(ax)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(ims[0], ax=axes[0], cmap=cmap, cax=cbar_ax)

def all_elect_reponses():
    '''
    Plot all combos of gain and peaking times.
    '''
    from electronics import response
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    times = numpy.linspace(0,10*units.us,5000)    
    for gain in [4.7, 7.8, 14.0, 25.0]:
        for peaking in [0.5, 1.0, 2.0, 3.0]:
            peaking = peaking*units.us
            resps  = response(times/units.ms, gain, peaking)
            ax.plot(times, resps)
    
    
def shaping_study(fields, gain=10.12, shaping=2.0*units.us, nbins=5000):
    import electronics
    from scipy.signal import fftconvolve

    nfields = len(fields)

    fig, axes = plt.subplots(nfields, 4, sharex=True)

    for row, field in enumerate(fields):

        newfrf = field.rebin(nbins)
        elecrf = electronics.response(numpy.linspace(*newfrf.domainls), gain, shaping)
        shaped = fftconvolve(newfrf.response, elecrf, "same")

        axes[row][0].plot(field.times/units.us,  field.response)
        axes[row][1].plot(newfrf.times/units.us, newfrf.response)
        axes[row][2].plot(newfrf.times/units.us, elecrf)
        axes[row][3].plot(newfrf.times/units.us, shaped)
        
    
    
