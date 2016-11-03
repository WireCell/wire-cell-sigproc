import units
import schema

import json
import numpy


def dumps(obj):
    '''
    Dump object to JSON text.
    '''
    return json.dumps(schema.todict(obj))

def loads(text):
    '''
    Load object from JSON text.
    '''
    return schema.fromdict(json.loads(text))

def dump(filename, obj):
    '''
    Save a response object (typically response.schema.FieldResponse)
    to a file of the given name.
    '''
    text = dumps(obj)
    if filename.endswith(".json"):
        open(filename,'w').write(text)
        return
    if filename.endswith(".json.bz2"):
        import bz2
        bz2.BZ2File(filename, 'w').write(text)
        return
    if filename.endswith(".json.gz"):
        import gzip
        gzip.open(filename, "wb").write(text)
        return
    raise ValueError("unknown file format: %s" % filename)

def load(filename):
    '''
    Return response.schema object representation of the data in the
    file of the given name.
    '''
    if filename.endswith(".json"):
        return loads(open(filename,'r').read())

    if filename.endswith(".json.bz2"):
        import bz2
        return loads(bz2.BZ2File(filename, 'r').read())

    if filename.endswith(".json.gz"):
        import gzip
        return loads(gzip.open(filename, "rb").read())
    
    raise ValueError("unknown file format: %s" % filename)



class Persist(object):
    '''
    Persist (read/write) a set of responses.

    This produces data objects and serializes them following the Wire
    Cell Toolkit field response schema.
    '''

    

    @property
    def export_dict(self):
        '''
        As dictionary for export.

            - response :: the sampled values of the induced current.
              Current in system of units

            - period :: the sampling period of the induced current.
              Time should be in microseconds.

            - start :: the absolute starting location of the drift
              path expressed as the distances in the (pitch,
              antidrift, [wire]) directions.  Note, this is NOT
              (x,y,z).  Distance should be in millimeters.  The 3rd
              element (wire direction) is optional (eg, in the case of
              2D fields).

            - toffset :: the time at which the drift path started.
              Time is in microseconds.

            - plane :: the wire plane identifier as a letter {u, v, w}
              that is holding wire-of-interest.

            - wire :: the absolute location of the wire nearest to the
              start expressed as the distances in the (pitch,
              antidrift) directions from the wire of interest.
              Distance is in millimeters.

            - region :: the wire region number in which the path
              starts counting from the wire-of-interest (which is in
              wire region 0).

            - impact :: the distance from nearest wire to the
              projection of the start point to the wire plane.
              Distance is in millimeters.

        Externally specified:

            - map from plane letter to wire angle

            - nominal drift velocity

            - locations of wires besides the wire-of-interest and the nearest-wire.

        fixme: move this definition into some documentation.
        '''
        return dict(response = self.response.tolist(),
                    period = (self.times[1] - self.times[0])/units.us,
                    start = ((self.pos[0] + self.impact)/units.mm, 10*units.cm/units.mm),
                    toffset = self.times[0]/units.us,
                    plane = self.plane.lower(),
                    wire = (self.pos[0]/units.mm, self.pos[1]/units.mm),
                    region = self.region,
                    impact = self.impact/units.mm)

