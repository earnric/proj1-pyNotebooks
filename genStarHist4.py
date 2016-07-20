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

    zCrit = 2.0e-7
    print "Normalize Z data..."
    s.g['metal'][s.g['metal']<zCrit]    = 1e-10 # Since we divide by Z below, don't use 0.0
    s.g['pgf'][s.g['pgf']>(1.0-1e-10)]  = 1.0
    s.g['pzf'][s.g['pzf']<zCrit]        = 0.0

    s.s['metal'][s.s['metal']<zCrit]   = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-10)] = 1.0
    s.s['pzf'][s.s['pzf']<zCrit]       = 0.0

    # Fix birth times
    print s.derivable_keys()[1:5]
    pynbody.analysis.ramses_util.get_tform(s,"/Users/earnric/bin/part2birth");
   
    # Create the base name for plots
    plotfile = "3Mpc_z%.2lf"%z

    # Loop over the bins and grab the star particles
    # with the metallicity in the bin range.
    # Multipling by the sp mass will give us the mass
    # of pristing stars in that bin.
    bins = np.logspace(-10, 0, 51)
    psm  = np.zeros(len(bins)-1) # Pristine Stellar mass in that bin
    tsm  = np.zeros(len(bins)-1) # total mass
    primsm  = np.zeros(len(bins)-1) # Primordial stellar mass
    for indx2,i in enumerate(bins):
        if indx2 < len(bins)-1:
            cond = (s.s['metal'] >= i) & (s.s['metal'] < bins[indx2+1]) # This selects for sp's in the bin
            psm[indx2] = np.sum(s.s['ppf'][cond] * s.s['mass'][cond])
            tsm[indx2] = np.sum(s.s['mass'][cond])
            # For sp's that are in our bin (above):
            #   Compute the polluted fraction (1-ppf)
            #   Compute the fraction of pristine metals: pzf/Z
            #   Compute the mass of polluted stars that are polluted only by pristine metals
            primsm[indx2] = np.sum((1.0-s.s['ppf'][cond]) * (s.s['pzf'][cond] / s.s['metal'][cond]) * s.s['mass'][cond])
            print "Bin left edge %.3e has mass %.4e" % (i,primsm[indx2])
    
    np.savetxt("primordSM_z%.2lf"%z,primsm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("pristSM_z%.2lf"%z,psm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("totSM_z%.2lf"%z,tsm,comments='') # comments='' gets rid of "#" in output

    print "Plot pristine sp mass in Z bins..."
    titleStr = "Pristine sp mass in sp Z bins @ z = %.3lf" % z
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(titleStr)
    ax.set_xlabel("Z")
    ax.set_ylabel("$M_{sun}$")
    plt.gca().set_xscale("log")
    #plt.gca().set_yscale("log")
    xpsm = np.ma.log10(psm)
    plt.fill_between(bins[0:50],0,psm)
    plt.plot(bins[0:50], psm,'bo')
#        plt.ylim([0,1])
    plt.savefig(plotfile+'.spPristMass-Zhist.png')

    print "Plot Primordial sp mass in Z bins..."
    titleStr = "Primordial sp mass in sp Z bins @ z = %.3lf" % z
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(titleStr)
    ax.set_xlabel("Z")
    ax.set_ylabel("$M_{sun}$")
    plt.gca().set_xscale("log")
    #plt.gca().set_yscale("log")
    xprimsm = np.ma.log10(primsm) # Use numpy math on the array. Handles zero values
    plt.fill_between(bins[0:50],0,primsm)
    plt.plot(bins[0:50], primsm,'ro')
    plt.ylim([0,1.0e7])
    plt.savefig(plotfile+'.spPrimordMass-Zhist.png')

    print "Plot total sp mass in Z bins..."
    titleStr = "Total sp mass in sp Z bins @ z = %.3lf" % z
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    ax.set_title(titleStr)
    ax.set_xlabel("Z")
    ax.set_ylabel("$M_{sun}$")
    plt.gca().set_xscale("log")
#    plt.gca().set_yscale("log")
    xtsm = np.ma.log10(tsm)
    plt.fill_between(bins[0:50],0,tsm)
    plt.plot(bins[0:50], tsm,'go')
#    plt.ylim([0,3.0e7])
    plt.savefig(plotfile+'.spTotalMass-Zhist.png')

    del s
    gc.collect()
    gc.collect()