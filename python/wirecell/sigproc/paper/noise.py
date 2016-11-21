from wirecell.sigproc import units, garfield, response, plots
import numpy

garfield_tarball = "/home/bviren/projects/wire-cell/garfield-data/ub_10.tar.gz"

def figure_adc():
    dat = garfield.load(garfield_tarball)    
    uvw = response.line(dat)
    fig,data = plots.plot_digitized_line(uvw)
    fig.savefig('paper-noise-figure-adc.pdf')

    with open('paper-noise-figure-adc.txt','w') as fp:
        for t,u,v,w in data:
            fp.write('%f %e %e %e\n' % (t,u,v,w))
            



if '__main__' == __name__:
    figure_adc()
    
