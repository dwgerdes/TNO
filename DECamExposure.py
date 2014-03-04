#!/usr/bin/env python

import ephem
from DECamField import DECamField

class DECamExposure(object):
#
    def __init__(self, expnum=0, UTobs='2013-01-01 00:00:00', exptime=0, band='r', ra=ephem.degrees(0), dec=ephem.degrees(0), tag='None', obj='None'):
        self.expnum = expnum
        self.UTobs = ephem.date(UTobs)
        self.exptime = exptime
        self.band = band
        self.ra = ephem.degrees(ra)
        self.dec = ephem.degrees(dec)
        self.tag = tag
        self.obj = obj
    
    def contains(self, ra1, dec1):
        # returns True if the point (ra1, dec1) lies inside the field
        return DECamField(self.ra, self.dec).contains(ra1, dec1)
    
    def ellipse(self):
        return DECamField(self.ra, self.dec).ellipse()
    
def main():
    pass

if __name__=="__main__":
    main()