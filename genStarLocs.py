#
# Compute histogram of metallicity (weighted by sp mass)
# 31 Aug 2015
# Rick Sarmento
#
# Purpose:
#  Read in RAMSES data and compute a Z histogram weighted
#  by the particles' pristine fraction
#
# Method:
#
# Revision history
#

import os
import re # Regular expression matching
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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
files.sort()
#print "All files %s"%files
# Here are the files that correspond to the outputs we want... 
#files = files[3:len(files)]
files = [
"output_00007",                                                                                                             
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
"output_00066",                                                                                         
"output_00068",
"output_00074",
"output_00084",
"output_00097",
"output_00104"
]

print "Files to process %s"%files
headerLine = "z,boxsize_kpc"
key        = np.empty([len(files),2]) 

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

    zCrit = 2.0e-7
    print "Normalize Z data..."
    s.g['metal'][s.g['metal']<zCrit]    = 1e-10 # Since we divide by Z below, don't use 0.0
    s.g['pgf'][s.g['pgf']>(1.0-1e-10)]  = 1.0
    s.g['pzf'][s.g['pzf']<zCrit]        = 0.0

    s.s['metal'][s.s['metal']<zCrit]   = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-10)] = 1.0
    s.s['pzf'][s.s['pzf']<zCrit]       = 0.0

    # # Fix birth times
    # print s.derivable_keys()[1:5]
    # pynbody.analysis.ramses_util.get_tform(s,"/Users/earnric/bin/part2birth");
    print "Saving data..."
    np.savetxt("spLoc_%.2lf.txt"%z,s.s['pos'],header="x\ty\tz",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spMass_%.2lf.txt"%z,s.s['mass'],header="mass",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spZ_%.2lf.txt"%z,s.s['metal'],header="Z",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spPPF_%.2lf.txt"%z,s.s['ppf'],header="ppf",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spPZF_%.2lf.txt"%z,s.s['pzf'],header="pzf",comments='') # comments='' gets rid of "#" in output
    print "Files saved"
    key[indx] = ["%.2lf"%z,"%.2lf"%s.properties['boxsize'].in_units('kpc')]
    del s

np.savetxt("zKeysForSPfiles.txt",key,header="z",comments='')
