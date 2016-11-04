import units
import schema

import json
import numpy


def dumps(obj):
    '''
    Dump object to JSON text.
    '''
    return json.dumps(schema.todict(obj), indent=2)

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

