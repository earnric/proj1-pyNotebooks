#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 09 Sep 2015
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and computes the various
#  pynbody as well as star-location plots.
#
# Method:
#  Sample 3 sps from the s.s collection and plot Z, PGF,
#  Temp, Rho and locations for each... 
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
pynbody.ramses.multiprocess_num = 32
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
"output_00068",
"output_00074",
"output_00084",
"output_00097",
"output_00104"
]
files.sort()
print "Files to process %s"%files

# Initialize array with sfr, 2D array: (z,sfr)
print "Files to process %d"%len(files)

np.set_printoptions(precision=3,suppress=True)

# Loop and load data
# If we wanted to just loop over filename - "for file in files:"
for indx, file in enumerate(files):
    gc.collect()
    print "Processing file %s" % file
    s = pynbody.load("./"+file)
    s['pos'] -= 0.5
    s.physical_units();
    
    z = 1/s.properties['a']-1
    print "Redshift = %.2lf" % z
    boxsizestring = "%.2lf kpc" % s.properties['boxsize'].in_units('kpc')
    print "Boxsize @ this redshift %s"%boxsizestring
    print "aexp %.3lf" % s.properties['a']
    bs = float(s.properties['boxsize'])

    z_crit = 2.0e-7

    sbox = 40.0 / ( 1.0 + z) / 0.71 # fix this! Needs to be the same we use for in genStarPlots6.py
    
    step = len(s.s)/3
    if step == 0: step = 1
    for i in range(0,len(s.s),step):
        x,y,zz = s.s['pos'][i] # Get the last one
        print "Star offset [%.2lf %.2lf %.2lf]"%(x,y,zz)
        print "Boxsize ", bs
        # Ensure we center such that we can depict a 40 kpc box
        if (abs(x) > bs/2.0 - sbox):
            x = np.sign(x) * (bs/2.0 - sbox)
        if (abs(y) > bs/2.0 - sbox):
            y = np.sign(y) * (bs/2.0 - sbox)
#        if (abs(zz) > bs/2.0 - sbox):
#            zz = np.sign(zz) * (bs/2.0 - sbox)
                    
        # Create the base name for plots
        coords= [x,y,zz] # NOTE WE ARE TRANSLATING... May not be centered on sp
        print coords,i

        np.savetxt("z%.2lf_SpCoord_%i.txt"%(z,i), coords,comments='') # comments='' gets rid of "#" in output
        
        gc.collect()

    del s
    gc.collect()
    gc.collect()

print "Done"
