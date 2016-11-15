from wirecell.sigproc import units, garfield, response, plots

garfield_tarball = "/home/bviren/projects/wire-cell/garfield-data/ub_10.tar.gz"

def figure_adc():
    dat = garfield.load(garfield_tarball)    
    uvw = response.line(dat)
    fig = plots.plot_digitized_line(uvw)
    fig.savefig('paper-noise-figure-adc.pdf')

if '__main__' == __name__:
    figure_adc()
    
