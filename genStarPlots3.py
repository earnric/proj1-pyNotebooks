#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 13 May 2015
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and computes the SFR for
#  the star particles formed during that step. The units
#  are Msol / yr / Mpc^3
#
# Method:
#  Sum the mass of star particles created, determine the
#  time span of star particle formation:
#       ( min t_form - max t_form )
#  and divide. Normalize by boxsize (nominally 3Mpc on a
#  side, so 9 Mpc^3.
#
# Revision history
#

#%matplotlib inline
import os
import re # Regular expression matching
import matplotlib as mpl
mpl.use('Agg') # Needed if you don't have an X-server
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import pynbody.plot as pp
import pynbody.plot.sph as sph
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import sys, os, glob, pickle, gc
import math
#import mmap

mpl.rcParams['figure.figsize'] = (12,8)
mpl.rcParams['font.size'] = 18
pynbody.ramses.multiprocess_num = 16
pynbody.config['number_of_threads'] = 128

# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
files = [f for f in os.listdir('.') if re.match(r'output_000[0-9][0-9]', f)]
files =["output_00007",
"output_00008",
"output_00011",
"output_00016",
"output_00020",
"output_00026",
"output_00033",
"output_00037",
"output_00043",
"output_00050",
"output_00058",
"output_00060",
"output_00068"
]
#files.sort()
print "Files to process %s"%files

# Initialize array with sfr, 2D array: (z,sfr)
print "Files to process %d"%len(files)

# Loop and load data
# If we wanted to just loop over filename - "for file in files:"
for indx, file in enumerate(files):
    gc.collect()
    print "Processing file %s" % file
    s = pynbody.load("./"+file,maxlevel=13)
    s['pos'] -= 0.5
    s.physical_units();
    
    z = 1/s.properties['a']-1
    print "Redshift = %.2lf" % z
    boxsizestring = "%.2lf kpc" % s.properties['boxsize'].in_units('kpc')
    print "Boxsize @ this redshift %s"%boxsizestring
    print "aexp %.3lf" % s.properties['a']
    bs = float(s.properties['boxsize'])

    print "Normalize Z data..."
    s.g['metal'][s.g['metal']<1e-7]   = 1e-10 # Since we divide by Z below, don't use 0.0
    s.g['pgf'][s.g['pgf']>(1.0-1e-7)] = 1.0
    s.g['pzf'][s.g['pzf']<1e-7]       = 0.0

    s.s['metal'][s.s['metal']<1e-7]   = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-7)] = 1.0
    s.s['pzf'][s.s['pzf']<1e-7]       = 0.0

    z_crit = 2.0e-7
    # Fix birth times
    print s.derivable_keys()[1:5]
    pynbody.analysis.ramses_util.get_tform(s,"/home1/02744/earnric/bin/part2birth");

    if (len(s.s) > 0) :
        x,y,zz =  s.s['pos'][1] # Get the coords of the first SP
        print "Star Particle 1 coords: ", x,y,zz
        
        coords= [0, 0, -zz]
        with pynbody.transformation.translate(s,coords):
            try:
                titleStr = "Z - z = %.1lf" % z + "\ncoords %s"%str(coords)
                sph.image(s.g,qty="metal",width=boxsizestring,cmap="nipy_spectral", denoise=True ,av_z=False,
                          log=True, vmax=1.0, vmin=1e-7, approximate_fast=False,title=titleStr,
                          filename=fileOut); #vmin=0.006, vmax=1.0,
            except:
                print "EXCEPTION while plotting Z file: %s"%file
                pass
            try:
                fileOut="PGF-z=%.1lf%s.png"% (z,str(coords))
                titleStr = "PGF - z = %.1lf" % z + "\ncoords %s"%str(coords)
                sph.image(s.g,qty="pgf",width=boxsizestring,cmap="nipy_spectral", denoise=True ,av_z=False,
                          log=True, vmin=1e-7, approximate_fast=False,title=titleStr,
                          filename=fileOut); #vmin=0.006, vmax=1.0,
            except:
                print "EXCEPTION while plotting PFG file: %s"%file
                pass
            try:
                fileOut="Density-z=%.1lf%s.png"% (z,str(coords))
                titleStr = "Density - z = %.1lf" % z + "\ncoords %s"%str(coords)
                sph.image(s.g,qty="rho",width=boxsizestring,cmap="nipy_spectral", denoise=True ,av_z=False,
                          log=True, approximate_fast=False,title=titleStr,
                          filename=fileOut); #vmin=0.006, vmax=1.0,
            except:
                print "EXCEPTION while plotting density file: %s"%file
                pass
            try:
                fileOut="Temp-z=%.1lf%s.png"% (z,str(coords))
                titleStr = "Temp - z = %.1lf" % z + "\ncoords %s"%str(coords)
                sph.image(s.g,qty="temp",width=boxsizestring,cmap="nipy_spectral", denoise=True ,av_z=False,
                          log=True, approximate_fast=False,title=titleStr,
                          filename=fileOut); #vmin=0.006, vmax=1.0,
            except:
                print "EXCEPTION while plotting temp file: %s"%file
                pass
            try:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title("Star particle locations z=%.1f\n2D projection, xy plane, %s boxsize"%(z,boxsizestring))
                plt.plot(s.s['x'],s.s['y'],'o')
                # x1,x2,y1,y2 = plt.axis()
                axes().set_aspect('equal')
                plt.axis((-bs/2.0,bs/2.0,-bs/2.0,bs/2.0))
                plt.savefig("SP-locsz%.1f.png"%z)
            except:
                print "EXCEPTION while plotting sp locations file: %s"%file
                pass
    # end if len(s.s) > 0
    del s
    gc.collect()
    gc.collect()
