#!/usr/bin/env python
'''
Process Garfield field response output files to produce Wire Cell
field response input files.
'''
import os.path as osp

import numpy
import matplotlib.pyplot as plt
import tarfile

def fromtarfile(filename):
    '''
    Iterate on tarfile, returning (name,text) pair of each file.
    '''

    tf = tarfile.open(filename, 'r')
    for name,member in sorted([(m.name,m) for m in tf.getmembers()]):
        if member.isdir():
            continue
        yield (member.name, tf.extractfile(member).read())

def split_text_records(text):
    for maybe in text.split("\n% "):
        if maybe.startswith("Created"):
            yield maybe

def parse_text_record(text):
    '''
    Iterate on garfield text, returning one record.
    '''
    lines = text.split('\n')

    ret = dict()

    created = lines[0].split()
    ret['created'] = '%s %s' %(created[1], created[3])
    ret['signal'] = None
    if 'Direct signal' in lines[0]:
        ret['signal'] = 'collection'
    if 'Cross-talk' in lines[0]:
        ret['signal'] = 'induction'

    ret['group'] = int(lines[2].split()[1])

    wire = lines[3].split()
    ret['wire'] = int(wire[1])
    ret['label'] = wire[4]

    pos = map(float, wire[6].split('=')[1][1:-1].split(','))
    ret['wirepos'] = tuple([10.0*p for p in pos]) # save as mm
    ret['wirepos_unit'] = 'mm'
    ret['voltage'] = float(wire[9])

    ret['nbins'] = nbins = int(lines[4].split()[4])

    xunit, yunit = lines[5].split(":")[1].split(",")
    xunit = [x.strip() for x in xunit.split("in")]
    yunit = [y.strip() for y in yunit.split("in")]

    xscale = 1.0 # float(lines[7].split("=")[1]);
    if "micro second" in xunit[1]:
        xscale *= 1.0e-6
    ret['x_unit'] = 's'

    yscale = 1.0 # float(lines[8].split("=")[1]);
    if "micro Ampere" in yunit[1]:
        yscale = 1.0e-6
    ret['y_unit'] = 'V'

    ret['xlabel'] = xunit[0]
    ret['ylabel'] = yunit[0]

    xdata = list()
    ydata = list()
    for line in lines[9:9+nbins]:
        xy = line[4:].split()
        xdata.append(float(xy[0]))
        ydata.append(float(xy[1]))
    if nbins != len(xdata) or nbins != len(ydata):
        raise ValueError('parse error for "%s"' % wire)
    ret['x'] = numpy.asarray(xdata)*xscale
    ret['y'] = numpy.asarray(ydata)*yscale
    return ret

def asgenerator(source):
    '''
    If string, assume file, open proper generator, o.w. just return
    '''
    if type(source) not in [type("") or type(u"")]:
        return source
    if osp.splitext(source)[1] in [".tar", ".gz", ".tgz"]:
        return fromtarfile(source)
    raise ValueError('unknown garfield data source: "%s"' % source)


def parse_filename(filename):
    '''
    Try to parse whatever data is encoded into the file name.
    '''
    fname = osp.split(filename)[-1]
    dist, plane = osp.splitext(fname)[0].split('_')
    plane = plane.lower()
    if plane == 'y':
        plane = 'w'
    return dict(impact=float(dist), plane=plane, filename=filename)

def load(source):
    '''
    Load Garfield data source (eg, tarball).

    Return list of per-drift-path dictionaries holding response.
    '''
    source = asgenerator(source)

    ret = list()
    from collections import defaultdict
    uniq = defaultdict(dict)
    for filename, text in source:

        fnamedat = parse_filename(filename)

        plane_letter = None
        for get,want in zip('uvy','uvw'):
            if get+'.dat' in filename.lower():
                plane_letter = want

        gen = split_text_records(text)
        for rec in gen:
            dat = parse_text_record(rec)

            key = tuple([filename] + [dat[k] for k in ['group', 'wire', 'label']])

            old = uniq.get(key, None)
            if old:             # sum up all signal types
                old['y'] += dat['y']
                continue

            dat.pop('signal')                
            dat.update(fnamedat)
            uniq[key] = dat

    return [u[1] for u in sorted(uniq.items())]

def average(dat):
    '''
    Return field response averaged over wire regions.

    Return as tuple (u,v,w) where each is a 2D array shape: (#regions, #responses).

    No time array is returned.  Take it from, eg dat[0]['x'].
    '''
    ret = list()
    for plane in 'uvw':
        byplane = [d for d in dat if d['plane'] == plane]
        responses = list()

        this_plane = list()
        for wire in sorted(set([d['wire'] for d in byplane])):
            bywire = [d for d in byplane if d['wire'] == wire]
            res = list()
            for impact, one in sorted([(d['impact'],d) for d in bywire]):
                res.append(one['y'])
            tot = numpy.zeros_like(res[0])
            for one in res:
                tot += one
            tot *= 2.0          # flip responses onto other side of wire
            tot -= 0.5*(res[0] + res[1]) # don't double count edge paths
            this_plane.append(tot)
        ret.append(numpy.vstack(this_plane))
    return tuple(ret)

def plot_average(avgtriple, time):
    '''
    Plot averages
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


