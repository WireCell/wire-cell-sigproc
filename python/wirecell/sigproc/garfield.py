#!/usr/bin/env python
'''
Process Garfield field response output files to produce Wire Cell
field response input files.

Garfield input is provided as a tar file.  Internal structure does not
matter much but the files are assumed to be spelled in the form:

<impact>_<plane>.dat

where <impact> spells the impact position in mm and plane is from the
set {"U","V","Y"}.

Each .dat file may hold many records.  See parse_text_record() for
details of assumptions.
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
    '''
    Return a generator that splits text by record separators.
    '''
    for maybe in text.split("\n% "):
        if maybe.startswith("Created"):
            yield maybe

def parse_text_record(text):
    '''
    Iterate on garfield text, returning one record.
    '''
    lines = text.split('\n')

    ret = dict()

    # Created 31/07/16 At 19.52.20 < none > SIGNAL   "Direct signal, group   1     "
    created = lines[0].split()
    ret['created'] = '%s %s' %(created[1], created[3])
    ret['signal'] = None
    if 'Direct signal' in lines[0]:
        ret['signal'] = 'collection'
    if 'Cross-talk' in lines[0]:
        ret['signal'] = 'induction'

    #   Group 1 consists of:
    ret['group'] = int(lines[2].split()[1])

    #      Wire 243 with label X at (x,y)=(-3,0.6) and at -110 V
    wire = lines[3].split()
    ret['wire_region'] = int(wire[1])
    ret['label'] = wire[4]

    pos = map(float, wire[6].split('=')[1][1:-1].split(','))
    ret['wire_region_pos'] = tuple([10.0*p for p in pos]) # save as mm
    ret['wire_region_pos_unit'] = 'mm'
    ret['bias_voltage'] = float(wire[9])

    #  Number of signal records:  1000
    ret['nbins'] = nbins = int(lines[4].split()[4])

    #  Units used: time in micro second, current in micro Ampere.
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
    #  + (  0.00000000E+00   0.00000000E+00
    #  +     0.10000000E+00   0.00000000E+00
    # ...
    #  +     0.99800003E+02   0.00000000E+00
    #  +     0.99900002E+02   0.00000000E+00 )
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

    Return list dictionaries of field response functions.

    Each dictionary has keys:

        - plane :: a letter from {"u","v","w"} giving the plane from
          which the wire-of-interest exists.

        - region :: the wire region number in which the path
          corresponding drift path exists with 0 indicating region
          around the wire-of-interest.

        - impact :: the distance in mm to the wire at the center of
          the given region.

        - response :: a Numpy array shaped (2, Nsamples).  First row
          holds the time values (seconds), second row hold the
          response function values (Ampere) for each sample.
    '''
    source = asgenerator(source)

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

            key = tuple([filename] + [dat[k] for k in ['group', 'wire_region', 'label']])

            old = uniq.get(key, None)
            if old:             # sum up all signal types
                old['y'] += dat['y']
                continue

            dat.pop('signal')                
            dat.update(fnamedat)
            uniq[key] = dat

    ret = list()
    for plane in 'uvw':
        byplane = [one for one in uniq.values() if one['plane'] == plane]
        zeros = [one for one in byplane if one['wire_region_pos'][0] == 0.0 and one['impact'] == 0.0]
        if len(zeros) != 1:
            raise ValueError("got too many zeros: %d" % len(zeros))
        zero_wire_region = zeros[0]['wire_region']
        for one in byplane:
            d = dict(plane=plane,
                     region = one['wire_region']-zero_wire_region,
                     impact = one['impact'],
                     response = numpy.asarray((one['x'], one['y'])))
            ret.append(d)
    return ret


def average(fine):
    '''
    Average fine-grained responses over multiple impact positions.  

    Return list dictionaries of field response functions with keys:

        - plane :: a letter from {"u","v","w"} giving the plane from
          which the wire-of-interest exists.

        - region :: the wire region number in which the path
          corresponding drift path exists with 0 indicating region
          around the wire-of-interest.

        - response :: a Numpy array shaped (2, Nsamples).  First row
          holds the time values (seconds), second row hold the
          response function values (Ampere) for each sample.
    '''
    ret = list()
    for plane in 'uvw':
        byplane = [d for d in fine if d['plane'] == plane]
        for region in set([d['region'] for d in byplane]):
            resp = list()
            inregion = [d for d in byplane if d['region'] == region]
            for d in inregion:
                resp.append(d['response'])
            tot = numpy.zeros_like(resp[0][1])
            for one in resp:
                tot += one[1];
            tot *= 2.0        # flip responses onto other side of wire
            tot -= 0.5*(resp[0][1] + resp[-1][1]) # don't double count edge paths
            dat = dict(plane=plane,
                       region = region,
                       response = numpy.asarray((resp[0][0], tot)))
            ret.append(dat)
        continue
    return ret

def dumps(dat):
    '''
    Return JSON string for data.
    '''
    import json
    tmp = list()
    for d in dat:
        d = dict(d)
        d['response'] = d['response'].tolist()
        tmp.append(d)
    return json.dumps(tmp)


def toarrays(dat):
    '''
    Return field response current waveforms as 3 2D arrays.

    Return as tuple (u,v,w) where each is a 2D array shape: (#regions, #responses).

    '''
    ret = list()
    for plane in 'uvw':
        byplane = [d for d in dat if d['plane'] == plane]
        this_plane = list()
        for region in sorted(set([d['region'] for d in byplane])):
            byregion = [d for d in byplane if d['region'] == region]
            if len(byregion) != 1:
                raise ValueError("unexpected number of regions: %d" % len(byregion))
            wave = byregion[0]['response'][1]
            this_plane.append(wave)
        ret.append(numpy.vstack(this_plane))
    return tuple(ret)


def plot_by_region(avgtriple, time):
    '''
    Plot response functions as 1D graphs.
    '''
    nwires = avgtriple[0].shape[0]
    nwireshalf = nwires//2
    wire0 = nwires - nwireshalf - 1

    central_collection = avgtriple[2][wire0]
    central_sum = sum(central_collection)
    time *= 1.0e6
    
    fig, axes = plt.subplots(nwireshalf+1, sharex=True)

    for wire_offset in range(nwireshalf+1):
        ax = axes[wire_offset]
        ax.set_title('Wire region %d' % wire_offset);

        iwirep = wire0 + wire_offset
        iwirem = wire0 + wire_offset

        for iplane in range(3):
            plane = avgtriple[iplane]
            wirep = plane[iwirep]
            wirem = plane[iwirem]

            wave = wirep
            if wire_offset:
                wave = 0.5*(wirep + wirem)
            wave /= central_sum
            ax.plot(time, wave)
        


def plot_average_colz(avgtriple, time):
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


def convert(inputfile, outputfile = "wire-cell-garfield-response.json.bz2", average=True):
    '''
    Convert an input Garfield file pack into an output wire cell field response file.
    '''
    dat = load(inputfile)
    if average:
        dat = average(dat)
    text = dumps(dat)
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

    
