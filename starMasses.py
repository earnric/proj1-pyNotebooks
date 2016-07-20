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
import pynbody.plot.sph
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import sys, os, glob, pickle, gc
import math
#import mmap

mpl.rcParams['figure.figsize'] = (12,8)
mpl.rcParams['font.size'] = 18
pynbody.ramses.multiprocess_num = 8
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
files = ["output_00068"]

print "Files to process %s"%files

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

    # Fix birth times
#    print s.derivable_keys()[1:5]
#    pynbody.analysis.ramses_util.get_tform(s,"/Users/earnric/bin/part2birth");
   
    # Create the base name for plots
#    plotfile = "3Mpc_z%.2lf"%z

#    print "Num stars %i"%len(s.s)
    #print "Star masses %.2lf"%s.s['mass']
    
    np.savetxt("SM_z%.2lf.txt"%z,s.s['mass'],comments='') # comments='' gets rid of "#" in output
    np.savetxt("SM_z%.2lf-Msol.txt"%z,s.s['mass'].in_units('Msol'),comments='') # comments='' gets rid of "#" in output
   
    del s
    gc.collect()
    gc.collect()
