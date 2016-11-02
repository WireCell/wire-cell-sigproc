#!/usr/bin/env python
'''
This module defines an object schema which strives to be generic
enough to describe various sets of field responses including:

    - 0D :: responses which do not extend to inter-pitch regions and
      have no intra-pitch variation.  Responses come from averaged 2D
      field calculations, eg with Garfield.  This is the type
      originally used for LArSoft simulation and deconvoultion for
      some time.

    - 1D :: responses defined on drift paths starting from
      fine-grained points on a line perpendicular to drift and wire
      directions and which spans multiple wire regions.  Responses
      come from 2D field calculations, eg with Garfield.  This is the
      type used in the Wire Cell simulation as developed by Xiaoyue Li
      and Wire Cell deconvolution as developed by Xin Qian.

    - 2D :: responses defined on drift paths starting from
      fine-grained points on a plane perpendicular to nominal drift
      direction and spanning multiple wire regions.  Responses come
      from 3D field calculations, eg with LARF.  Simulation and
      deconvolution using these type of responses are not yet
      developed.

The schema is defined through a number of `namedtuple` collections.

Units Warning: time in microseconds, distance in millimeters, current
in Amperes.  This differs from the Wire Cell system of units!
'''

from collections import namedtuple

class FieldResponse(namedtuple("FieldResponse","planes axis origin tstart period")):
    __slots__ = ()
    '''
    :param list planes: List of PlaneResponse objects.
    :param list axis: A normalized 3-vector giving direction of axis (anti)parallel to nominal drift direction.
    :param float origin: location in millimeters on the axis where drift paths begin.
    :param float tstart: time in microseconds at which drift paths begin.
    :param float period: the sampling period in microseconds.
    '''

class PlaneResponse(namedtuple("PlaneResponse","paths planeid pitchdir wiredir")):
    __slots__ = ()
    '''
    :param list paths: List of PathResponse objects.
    :param int planeid: A numerical identifier for the plane.
    :param list pitchdir: A normalized 3-vector giving direction of the wire pitch.
    :param list wiredir: A normalized 3-vector giving direction of the wire run.

    Along with FieldResponse.axis, the following should hold: axis X wiredir = pitchdir
    '''
    

class PathResponse(namedtuple("PathResponse", "current regionid impact pitchpos wirepos")):
    __slots__ = ()
    '''
    :param:
    '''
