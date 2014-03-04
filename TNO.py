#!/usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pylab
import sys
import ephem
import pywcs
import csv
from DECamField import DECamField
from DECamExposure import DECamExposure
from ccdBounds import ccdBounds

plots = False

ctio = ephem.Observer()
ctio.lon = ephem.degrees(-70.806525)
ctio.lat = ephem.degrees(-30.169661)
ctio.elevation = 2207


def define_footprint(polydef):
    ra = []
    dec = []
    
    with open(polydef) as f:
        for line in f:
            if line[0] != '#':
                s = line.split(' ')
                ra_i = -999
                for i in range(len(s)):
                    if s[i] != '' and ra_i==-999:
                        ra_i = float(s[i])
                        ra.append(ra_i*np.pi/180)
                    elif s[i] != '' and ra_i != -999:
                        dec.append(float(s[i])*np.pi/180)
    
    return zip(ra, dec)

def compute_TNOs(date):
# date is a pyephem date object
    comment='#'
    alltno=[]
    with open('Soft03Distant.txt') as cat:
        for line in cat:
            if line[0]!=comment: alltno.append(ephem.readdb(line))
    for transNeptunianFlyingRock in alltno:
        transNeptunianFlyingRock.compute(date)
    tno = [t for t in alltno if t.dec<ephem.degrees('05:00:00')]  # dec>5 is outside out footprint for sure.
    return tno
    
def point_inside_polygon(x,y,poly):

    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def SNfields():
    C1 = ephem.readdb("C1,f,03:37:05.83,-27:06:41.8,23.0,2000")
    C2 = ephem.readdb("C2,f,03:37:05.83,-29:05:18.2,23.0,2000")
    C3 = ephem.readdb("C3,f,03:30:35.62,-28:06:00.0,24.0,2000")
    X1 = ephem.readdb("X1,f,02:17:54.17,-04:55:46.2,23.0,2000")
    X2 = ephem.readdb("X2,f,02:22:39.48,-06:24:43.6,23.0,2000")
    X3 = ephem.readdb("X3,f,02:25:48.00,-04:36:00.0,24.0,2000")
    S1 = ephem.readdb("S1,f,02:51:16.80, 00:00:00.0,23.0,2000")
    S2 = ephem.readdb("S2,f,02:44:46.66,-00:59:18.2,23.0,2000")
    E1 = ephem.readdb("E1,f,00:31:29.86,-43:00:34.6,23.0,2000")
    E2 = ephem.readdb("E2,f,00:38:00.00,-43:59:52.8,23.0,2000")
    fields = [C1, C2, C3, X1, X2, X3, S1, S2, E1, E2]
    for f in fields:
        f.compute()
    return fields

def readExposures(expfile):
    exps = []
    with open(expfile,'rb') as f:
        reader = csv.reader(f, delimiter=',', quotechar='|')
        for line in reader:
            expnum = line[0]
            UTdate = line[1]
            exptime = line[2]
            band = line[3]
            ra = line[4]
            dec = line[5]
            tag = line[6]
            obj = line[7]
            dates = UTdate.split('T')
            expdate = dates[0]+' '+dates[-1]
            exps.append(DECamExposure(expnum, expdate, exptime, band, ra, dec, tag, obj))
    return exps

def compute_chip(rockra, rockdec, expra, expdec):
    deltara = 180/np.pi*ephem.degrees(rockra-expra).znorm  # compute difference in degrees (normalized between -180, +180)
    deltadec = 180/np.pi*ephem.degrees(rockdec-expdec).znorm  # the 180/pi is because ephem.Angle objects are natively in radians
    ccdname = 'None'
    for k in ccdBounds.keys():
        if deltara > ccdBounds[k][0] and deltara < ccdBounds[k][1] and deltadec > ccdBounds[k][2] and deltadec < ccdBounds[k][3]:
            ccdname = k
    return ccdname
    

# These rocks wandered into an SV exposure:
SV_TNO_names = ['87269 2000 OO67', '145474 2005 SA278', '145480 2005 TB190', '2006 QZ180', '2007 TZ417', '2007 TD418', '2007 VK302', '2007 VL305', '2010 GW64']


poly_13 = define_footprint('round13-poly.txt')    
foot = patches.Polygon(poly_13, facecolor='blue', edgecolor='none', alpha=0.3)
poly_sptw = define_footprint('sptw-poly.txt')
foot_sptw = patches.Polygon(poly_sptw, facecolor='none', edgecolor='k', alpha=0.6, linewidth=2)

# Draw the SN fields:
SNellipse = [DECamField(f.a_ra, f.a_dec).ellipse() for f in SNfields()]

rocks = compute_TNOs('2013/01/01')
SV_TNOs = [t for t in rocks if t.name in SV_TNO_names]
#for s in SV_TNOs: print s.name, s.mag, s.ra, s.dec, s.earth_distance

exposures = readExposures('SV_exposures.csv')
iexp = 0
Nexp = len(exposures)
for rock in SV_TNOs:
    rock.compute('2013/01/01')
    for exp in exposures:
#        print exp.ra, exp.dec
        ctio.date = exp.UTobs
        rock.compute(ctio)
        if exp.contains(ephem.degrees(rock.ra), rock.dec) and rock.mag < 25:
            ccd = compute_chip(ephem.degrees(rock.ra), rock.dec, exp.ra, exp.dec)
            print rock.name, ',', exp.UTobs, ',', ephem.degrees(rock.ra), ',', ephem.degrees(rock.dec), ',', rock.earth_distance, ',', rock.mag, ',', exp.expnum, ',', ccd, ',', exp.ra, ',', exp.dec, ',', exp.exptime, ',',exp.band, ',', exp.obj
            
            
def do_plots():  
    plt.figure(1, facecolor='w', edgecolor='k')
    #plt.suptitle('TNOs ')
    ax = plt.subplot(111, projection='mollweide')
    
    
    Nept = ephem.Neptune()
    d = ephem.date('2013/01/01')
    end_date = ephem.date('2013/01/02')
    while d < end_date:
        TNOs = compute_TNOs(d)
        TNO_ra = [t.ra for t in TNOs]
        TNO_dec = [t.dec for t in TNOs]
        TNO_dist = [t.earth_distance for t in TNOs]
        for i in range(len(TNO_ra)):
            if TNO_ra[i]>ephem.degrees('180:00:00'): TNO_ra[i]-=2*np.pi
          
        ax.scatter(TNO_ra, TNO_dec, color='r', marker='o', s=6, alpha=0.5)
        ax.add_patch(foot)
        ax.add_patch(foot_sptw)
        for e in SNellipse:    
            ax.add_artist(e)
            e.set_facecolor('g')
            e.set_alpha(0.75)
        for exp in exposures:
            if exp.ra > ephem.degrees('180:00:00'): exp.ra -= 2*np.pi
            ell = exp.ellipse()
            ax.add_artist(ell)
            ell.set_facecolor('b')
            ell.set_alpha(0.3)
        ax.grid(True)
    
        name = 'images/TNO'+str(int(d))+'.png'
    
        plt.show()
    #    plt.cla()  # clear axis
        d += 5     # skip ahead N days
        print d

if plots: do_plots() 