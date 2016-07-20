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
pynbody.ramses.multiprocess_num = 16
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
files = ["output_00007",
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
"output_00068"]

print "Files to process %s"%files

print "Number of files to process %d"%len(files)

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

    print "Normalize Z data..."
    s.g['metal'][s.g['metal']<1e-7]   = 1e-10 # Since we divide by Z below, don't use 0.0
    s.g['pgf'][s.g['pgf']>(1.0-1e-7)] = 1.0
    s.g['pzf'][s.g['pzf']<1e-7]       = 0.0

    s.s['metal'][s.s['metal']<1e-7]   = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-7)] = 1.0
    s.s['pzf'][s.s['pzf']<1e-7]       = 0.0

    # Fix birth times
    print s.derivable_keys()[1:5]
    pynbody.analysis.ramses_util.get_tform(s,"/home1/02744/earnric/bin/part2birth");
   
    # Create the base name for plots
    plotfile = "3Mpc_z%.2lf"%z

    print "Num stars %i"%len(s.s)
    
    print "Plot SP Z hist..."
    titleStr = "SP Z @ z = %.3lf" % z
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(titleStr)
    print "Right before hist plt"
    plt.hist(s.s['metal'],bins=np.logspace(-10, 0, 51))
    plt.gca().set_xscale("log")
#        plt.ylim([0,1])
    plt.savefig(plotfile+'.starParticle.Zhist.png')

    del s
    gc.collect()
    gc.collect()


