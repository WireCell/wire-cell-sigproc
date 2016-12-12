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


def main():
    cli(obj=dict())
