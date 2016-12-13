#!/usr/bin/env python
'''
Export wirecell.sigproc functionality to a main Click program.
'''

import click

@click.group("sigproc")
@click.pass_context
def cli(ctx):
    '''
    Wire Cell Signal Processing Features
    '''


@cli.command("convert-garfield")
@click.argument("garfield-fileset")
@click.argument("wirecell-field-response-file")
@click.pass_context
def convert_garfield(ctx, garfield_fileset, wirecell_field_response_file):
    '''
    Convert an archive of a Garfield fileset (zip, tar, tgz) into a
    Wire Cell field response file (.json with optional .gz or .bz2
    compression).
    '''
    import wirecell.sigproc.garfield as gar
    import wirecell.sigproc.response as res
    import wirecell.sigproc.response.persist as per

    rflist = gar.load(garfield_fileset)
    fr = res.rf1dtoschema(rflist)
    per.dump(wirecell_field_response_file, fr)


@cli.command("plot-field-response")
@click.argument("wcfrfile")
@click.argument("pdffile")
@click.pass_context
def plot_field_response(ctx, wcfrfile, pdffile):
    import wirecell.sigproc.response.persist as per
    import wirecell.sigproc.response.plots as plt

    fr = per.load(wcfrfile)
    # ...


@cli.command("plot-track-response")
@click.option("-o", "--output", default=None,
              help="Set output data file")
@click.option("-g", "--gain", default=14.7,
              help="Set gain.")
@click.option("-s", "--shaping", default=2.0,
              help="Set shaping time in us.")
@click.option("-t", "--tick", default=0.5,
              help="Set tick time in us (0.1 is good for no shaping).")
@click.argument("garfield-fileset")
@click.argument("pdffile")
@click.pass_context
def plot_track_response(ctx, output, gain, shaping, tick,
                            garfield_fileset, pdffile):
    import wirecell.sigproc.garfield as gar
    import wirecell.sigproc.response as res
    import wirecell.sigproc.plots as plots
    from wirecell.sigproc import units

    shaping *= units.us
    tick *= units.us

    rflist = gar.load(garfield_fileset)
    uvw = res.line(rflist)
    fig, data = plots.plot_digitized_line(uvw, gain, shaping, tick)
    fig.savefig(pdffile)

    if output:
        with open(output, 'w') as fp:
            for t, u, v, w in data:
                fp.write('%f %e %e %e\n' % (t, u, v, w))


def main():
    cli(obj=dict())
