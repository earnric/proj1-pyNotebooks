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
mpl.use('Agg') # Needed if you don't have an X-server
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
pynbody.ramses.multiprocess_num = 16
pynbody.config['number_of_threads'] = 64

# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
files = [f for f in os.listdir('.') if re.match(r'output_000[0-9][0-9]', f)]
files.sort()
print "All files %s"%files
# Here are the files that correspond to the outputs we want... 
files = files[3:len(files)]
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
"output_00060",
"output_00066",
"output_00068"]

print "Files to process %s"%files

# Initialize array with sfr, 2D array: (z,sfr)
data      = np.empty([len(files),9])
headerStr = 'z\ttStart\ttEnd\ttotalStarMass\ttotalPop3StarMass\ttotalPollStarMass\ttotalPrimordStarMass\ttotalGasMass\ttotalPristGasMass'

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

    print "Normalize Z data..."
    s.g['metal'][s.g['metal']<2e-7]   = 1e-10 # Since we divide by Z below, don't use 0.0
    s.g['pgf'][s.g['pgf']>(1.0-1e-7)] = 1.0
    s.g['pzf'][s.g['pzf']<2e-7]       = 0.0

    s.s['metal'][s.s['metal']<2e-7]   = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-7)] = 1.0
    s.s['pzf'][s.s['pzf']<2e-7]       = 0.0

    # Fix birth times
    print s.derivable_keys()[1:5]
    pynbody.analysis.ramses_util.get_tform(s,"/home1/02744/earnric/bin/part2birth");

    # Get total stellar mass, then total gas mass
    totalStarMass = np.sum(s.s['mass'].in_units('Msol'))
    print "Total stellar mass is %.2e Msol @ z~%.2lf" % (totalStarMass,z)
    totalPop3StarMass = np.sum(s.s['mass'].in_units('Msol') * s.s['ppf'])
    print "Total POPIII stellar mass is %.2e Msol @ z~%.2lf" % (totalPop3StarMass,z)
    totalPolStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0-s.s['ppf']))
    print "Total polluted stellar mass is %.2e Msol @ z~%.2lf" % (totalPolStarMass,z)

    # Compute the mass of stars with only primordial metals
    # 1 - ppf is the polluted fraction of stars, by mass
    # (1-ppf) * pzf is then the fraction polluted by primordial metals only
    # (1-ppf) * pzf / Z normalizes the fraction relative to the metallicity of the star particle
    totalPriStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0 - s.s['ppf']) * s.s['pzf']/ s.s['metal'])
    print "Total primordial stellar mass is %.2e Msol @ z~%.2lf" % (totalPriStarMass,z)

    totalGasMass  = np.sum(s.g['mass'].in_units('Msol'))
    print "Total gas mass is %.2e Msol @ z~%.2lf" % (totalGasMass,z)
    totalPristGasMass = np.sum(s.g['mass'].in_units('Msol') * s.g['pgf'])
    print "Total pristing gas mass is %.2e Msol @ z~%.2lf" % (totalPristGasMass,z)

    data[indx] = [z, min(s.s['tform'].in_units('yr')), max(s.s['tform'].in_units('yr')),
                  totalStarMass,totalPop3StarMass,totalPolStarMass,totalPriStarMass,totalGasMass,totalPristGasMass]
    
    # Create the base name for plots
    plotfile = "3Mpc_z%.2lf"%z
    
    try:
        print "SFH plot..."
        titleStr = "Star Form History @ z = %.2lf" % z
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(titleStr)
        pp.sfh(s,nbins=100) # Star Formation History
        plt.savefig(plotfile+'.sfh.png')

        print "SFH for Primordial Stars..."
        titleStr = "Primordial Star Form History @ z = %.2lf" % z
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(titleStr)
        plt.title("Primordial Star Formation History @ z=%.2lf"%z)
        pp.sfh_Z_p(s,nbins=100) # Star Formation History - primordial Z 
        plt.savefig(plotfile+'.sfhPrimordOnly.png')
    
        print "Plot Z..."
        titleStr = "Z @ z = %.3lf" % z + "\nVol ave"
        pp.sph.image(s.g,qty="metal",width=boxsizestring,cmap="nipy_spectral", denoise=True ,av_z=True, log=True, vmin=1e-7,
                approximate_fast=False,title=titleStr,qtytitle="$log_{10}\; Z$",filename=plotfile+'.Z.png');
        
        print "Plot Primord Z..."
        titleStr = "Primord_Z @ z = %.3lf" % z + "\nVol ave"
        pp.sph.image(s.g,qty="pzf",width=boxsizestring,cmap="nipy_spectral", denoise=True ,av_z=True, log=True, vmin=1e-7,
                approximate_fast=False,title=titleStr,qtytitle="$log_{10}\; Z_{pri}$",filename=plotfile+'.primordZ.png');
        
        print "Plot PGF..."
        titleStr = "Polluted gas fraction @ z = %.3lf" % z + "\nVol ave"
        s.g['ompgf'] = 1.0 - s.g['pgf']
        pp.sph.image(s.g,qty="ompgf",width=boxsizestring,cmap="nipy_spectral", denoise=True ,av_z=True, log=True, vmin=1e-7,
               approximate_fast=False,title=titleStr,qtytitle="$log_{10}\; 1-PGF$",filename=plotfile+'.1-pgf.png');
        
        print "Plot Rho..."
        titleStr = "rho @ z = %.3lf" % z + "\nVol ave"
        pp.sph.image(s.g,qty="rho",units="g cm^-3",width=boxsizestring,cmap="PuBuGn", denoise=True ,av_z=True, log=True, 
                approximate_fast=False,title=titleStr,filename=plotfile+'.rho.png');
        
        print "Plot DM Rho..."
        titleStr = "DM rho @ z = %.3lf" % z + "\nVol ave"
        pp.sph.image(s.dm,qty="rho",units="g cm^-3",width=boxsizestring,cmap="PuBuGn", denoise=True ,av_z=True, log=True, 
                approximate_fast=False,title=titleStr,filename=plotfile+'.DM-rho.png');
        
        print "Plot Temp..."
        titleStr = "temp @ z = %.3lf" % z + "\nVol ave"
        pp.sph.image(s.g,qty="temp",width=boxsizestring,cmap="YlGnBu", denoise=True ,av_z=True, log=True, 
                approximate_fast=False,title=titleStr,filename=plotfile+'.temp.png');

        print "Plot gas phase..."
        titleStr = "GAS PHASE @ z = %.2lf" % z
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(titleStr)
        phasePlot=pp.gas.rho_T(s.g,rho_units='g cm**-3', # rho_range=(rhoMin,rhoMax), t_range=(tempMin,tempMax),
                colorbar=True,  title=titleStr) # scalemin=0.1, scalemax=8.0e5,
        plt.savefig(plotfile+'.phase.png')
    
    except:
        print "EXCEPTION while processing file: %s"%file
        pass

    np.savetxt(file+'_data3Mpc.txt',data,comments='',header=headerStr) # comments='' gets rid of "#" in output
    del s
    gc.collect()
    gc.collect()

# Write out the array to a file
print data
np.savetxt("data3Mpc.txt",data,comments='',header=headerStr) # comments='' gets rid of "#" in output
del data


