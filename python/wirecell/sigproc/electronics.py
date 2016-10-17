#!/usr/bin/env python
import units
import numpy
from math import sin, cos, exp

# fixme: handle units correctly

def response(time, gain=10.12, shaping=2.0*units.us):
    if time <= 0.0 or time >= 10*units.us:
        return 0.0

    time = time/units.us
    st = shaping/units.us
        
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

response = numpy.vectorize(response)

