#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 12 Nov 2015
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and outputs masses for sp.
#
# Method:
#


# ##########################################################
# Output the data needed to generate the Mathematica histograms
# ##########################################################
def outputHistogramData(z,i,halo):
    bins    = np.logspace(-10, 0, 51) # Log bins for histogram data 
    psm     = np.zeros(len(bins)-1) # Pristine Stellar mass in that bin
    tsm     = np.zeros(len(bins)-1) # total mass
    primsm  = np.zeros(len(bins)-1) # Primordial stellar mass

    # By using zsolar, we are binning in units of solar Z.
    for indx2,j in enumerate(bins):
        if indx2 < len(bins)-1:
            cond = (halo.s['zsolar'] >= j) & (halo.s['zsolar'] < bins[indx2+1]) # This selects for sp's in the bin
            psm[indx2] = np.sum(halo.s['ppf'][cond] * halo.s['mass'][cond]) # Pop III
            tsm[indx2] = np.sum(halo.s['mass'][cond]) # Total star particle mass in this bin
            # For sp's that are in our bin (above):
            #   Compute the polluted fraction (1-ppf)
            #   Compute the fraction of pristine metals: pzf/Z
            #   Compute the mass of polluted stars that are polluted only by pristine metals
            # Note that 'metal' maybe 1e-10, but that means pzf = 0.0 and we get the correct
            # mass
            primsm[indx2] = np.sum((1.0-halo.s['ppf'][cond]) * 
                                   (halo['pzf'][cond] / halo.s['metal'][cond]) * halo.s['mass'][cond])
    # Save the total masses
    np.savetxt("primordSM_z%.1lf-%i.txt"%(z,i),primsm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("pristSM_z%.1lf-%i.txt"%(z,i),psm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("totSM_z%.1lf-%i.txt"%(z,i),tsm,comments='') # comments='' gets rid of "#" in output

    return

# ##########################################################
# ##########################################################
##
## Main program
##
# ##########################################################
# ##########################################################
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

mpl.rcParams['figure.figsize'] = (12,12)
mpl.rcParams['font.size'] = 22
pynbody.ramses.multiprocess_num = 32
pynbody.config['number_of_threads'] = 128

# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
files = [f for f in os.listdir('.') if re.match(r'output_000[0-9][0-9]', f)]
files =[
    ## "output_00007",
    ## "output_00008",
    ## "output_00011",
    "output_00016",
    ## "output_00020",
    ## "output_00026",
    ## "output_00033",
    ## "output_00037",
    ## "output_00043",
    ## "output_00050",
    ## "output_00058",
    ## "output_00066",
    ## "output_00068",
    ## "output_00073",
    ## "output_00084",
    ## "output_00097",
    ## "output_00107",
    "output_00121",
    ## "output_00136",
    ## "output_00152",
    ## "output_00169",
    "output_00191",
    ## "output_00214",
    "output_00241"
    ]
files.sort()
print "Files to process %d"%len(files)

np.set_printoptions(precision=3,suppress=True)

key = np.empty([len(files),2]) # a map of z to boxsize, 

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

    z_crit = 2.0e-7 # Mass fraction, 1e05 solar
    print "Normalize Z data..."

    s.s['metal'][s.s['metal'] < z_crit]    = 1e-10 # Since we divide by metal when computing pz mass don't use 0.0
    #s.s['ppf'][s.s['ppf'] > (1.0-z_crit)]  = 1.0
    if (len(s.s['ppf'][s.s['ppf'] < 0.0])):
        print("**** ERROR -- Prist particle fraction < 0.0")
    s.s['pzf'][s.s['pzf'] < z_crit]        = 0.0
    s.s['zsolar']  = s.s['metal'] * 50.0        # Solar units
    s.s['pzsolar'] = s.s['pzf'] * 50.0          # Solar units

    # Fix birth times
    print s.derivable_keys()[1:5]
    pynbody.analysis.ramses_util.get_tform(s,"~/bin/part2birth");
    
    sbox = 40.0 / (1.0 + z) * 0.71 # 40 kpc comoving box
    print "40 kpc comoving Small box @ z=%.2lf is %.2lf"%(z,sbox)

    step = len(s.s)/3 # Note integer division
    if step == 0: step = 1
    for i in range(0,len(s.s),step):
        x,y,zz = - s.s['pos'][i] # get a star, this is an offset (-1 * pos)
        print "Star offset [%.2lf %.2lf %.2lf]"%(x,y,zz)
        print "Boxsize @ this z ", bs
        # Ensure we center such that we can depict a 40 kpc box
        if (abs(x) > bs/2.0 - sbox/2.0):
            x = np.sign(x) * (bs/2.0 - sbox/2.0) # closest x coord to sp
        if (abs(y) > bs/2.0 - sbox/2.0):
            y = np.sign(y) * (bs/2.0 - sbox/2.0) # closest y coord to sp
#        if (abs(zz) > bs/2.0 - sbox/2.0): # Commented out so we get the exact z-coord of the sp
#            zz = np.sign(zz) * (bs/2.0 - sbox/2.0) # we could go off the end... ok
                    
        # Create the size of our box for sph image plots
        # May not be centered on sp - we have adjusted x,y such that we
        # don't go 'off the edge' for a 'smallbox' sized image
        coords= [x,y,zz] # NOTE WE ARE TRANSLATING... coords is -1 * location of sp of interest
        if sbox > 1.0: smallbox= str(sbox) + ' kpc'
        else: smallbox= str(sbox*1000.0) + ' pc'
        print "Plotsize: ",smallbox

        # ##########################################################
        # Save the index and coordinates of the center of our image...
        np.savetxt("z%05.2lf_SpCoord_%i.txt"%(z,i), -1.0*np.array(coords),fmt='%.2lf',comments='') 
        # ##########################################################

        rx,ry,rz = -1 * np.array([x,y,zz]) # For filters, we need the actual coords, not translation coords.
        print "Cuboid center [%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        # Note that the thickness in z is only 1/2 that of the other dimensions
        halo = s[pynbody.filt.Cuboid(str((rx-sbox/2.0)) + " kpc", str((ry-sbox/2.0)) + " kpc",str((rz-sbox/4.0)) + " kpc",
                                        str((rx+sbox/2.0)) + " kpc", str((ry+sbox/2.0)) + " kpc",str((rz+sbox/4.0)) + " kpc")]
        
        # Compute stellar masses and save for the halo histogram plot
        outputHistogramData(z,i,halo)

    del s
    gc.collect()
    gc.collect()

print "Done"
