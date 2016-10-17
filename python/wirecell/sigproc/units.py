#!/usr/bin/env python
'''
Provide Python access to Wire Cell system of units.

'''
#fixme: move this someplace more generic.

# distance
millimeter = mm = 1.0
micrometer = um = mm/1000.0
nanometer = nm = um/1000.0
centimeter = cm = mm*10.0
meter = m = cm*100.0
kilometer = km = m*1000.0

# time
nanosecond = ns = 1.0
microsecond = us = ns*1000.0
millisecond = ms = us*1000.0
second = s = ms*1000.0

# current
amp = 1.0
milliamp = amp/1000.0
microamp = milliamp/1000.0
nanoamp = microamp/1000.0

# ... add more as needed.  See and keep consistent with
# util/inc/WireCellUtil/Units.h.
