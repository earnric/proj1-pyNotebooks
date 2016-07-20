#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 11 Jan 2016
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and generates gas images.
#
#
# Revision history
#

# ##########################################################
# Generate and save images... 
# ##########################################################
def outputImages(z,i,snBlast):
    # Generate sph images for out box...

    with pynbody.transformation.translate(snBlast,coords):
        try:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            fileOut="img_log_Z-z=%.1lf-%i.png"% (z,i)
            titleStr = "$<Z>/Z_{\odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="zsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False, subplot=ax,
                      log=True, vmax=1.0, vmin=1e-7, qtytitle="$log_{10}\;<Z>/Z_{\odot}$", approximate_fast=False)
#            ,title=titleStr
#            filename=fileOut); #vmin=0.006, vmax=1.0,
            plt.savefig(fileOut)
            plt.close(fig)

            fig = plt.figure(rect)
            ax = fig.add_subplot(111)
            fileOut="img_log_PGF-z=%.1lf-%i.png"% (z,i)
            titleStr = "PGF - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="pgf",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,subplot=ax,
                      log=True, vmax = 1.0, vmin=1e-7, qtytitle="$log_{10}\; P$",approximate_fast=False)
#                      title=titleStr,
#                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            plt.savefig(fileOut)
            plt.close(fig)

            fig = plt.figure(rect)
            ax = fig.add_subplot(111)
            fileOut="img_log_PZ-z=%.1lf-%i.png"% (z,i)
            titleStr = "$Z_{P, \odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="pzsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,subplot=ax,
                      log=True, vmax=1.0, vmin=1e-7, qtytitle="$log_{10}\; <Z_{P}>/Z_{\odot}$",approximate_fast=False)
#                      title=titleStr,
#                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            plt.savefig(fileOut)
            plt.close(fig)

            fig = plt.figure(rect)
            ax = fig.add_subplot(111)
            fileOut="img_log_Density-z=%.1lf-%i.png"% (z,i)
            titleStr = "Density - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="rho",width=smallbox,cmap="terrain", denoise=True ,av_z=False,subplot=ax,
                      log=True, approximate_fast=False)
#                      title=titleStr,
#                      filename=fileOut); #vmin=0.006, vmax=1.0,
            
            plt.savefig(fileOut)
            plt.close(fig)

            fig = plt.figure(rect)
            ax = fig.add_subplot(111)
            fileOut="img_log_Temp-z=%.1lf-%i.png"% (z,i)
            titleStr = "Temp - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="temp",width=smallbox,cmap="RdYlBu", denoise=True ,av_z=False,subplot=ax,
                      log=True, approximate_fast=False)
            #,title=titleStr,
            #          filename=fileOut); #vmin=0.006, vmax=1.0,
            
            plt.savefig(fileOut)
            plt.close(fig)

            fig = plt.figure(rect)
            ax = fig.add_subplot(111)
            fileOut="img_log_vt-z=%.1lf-%i.png"% (z,i)
            titleStr = "$v_{t}$ - z = %.1lf" % z # + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
            print (titleStr)
            sph.image(snBlast.g,qty="tv",width=smallbox,cmap="RdYlBu", denoise=True ,av_z=False,subplot=ax,
                      log=True, qtytitle="$log_{10}\; v_{t} [km/s]$",approximate_fast=False)
            #,title=titleStr,
            #          filename=fileOut); #vmin=0.006, vmax=1.0,
            
            plt.savefig(fileOut)
            plt.close(fig)

            gc.collect()
        except:
            print ("EXCEPTION while processing file: %s"%file)
            print ("Plot title: %s"%titleStr)
            gc.collect()
            pass
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

mpl.rcParams['figure.figsize'] = (12,8)
mpl.rcParams['font.size'] = 18
pynbody.ramses.multiprocess_num = 12
pynbody.config['number_of_threads'] = 64

#rect = [0.15,0.15,0.85,0.9]
# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
files = [f for f in os.listdir('.') if re.match(r'output_000[0-9][0-9]', f)]
files =[
    # "output_00007",
    # "output_00008",
    # "output_00011",
    "output_00016",
    # "output_00020",
    # "output_00026",
    # "output_00033",
    # "output_00037",
    # "output_00043",
    # "output_00050",
    # "output_00058",
    # "output_00066",
    # "output_00068",
    # "output_00073",
    # "output_00084",
    # "output_00097",
    # "output_00107",
    "output_00121",
    #"output_00136",
    #"output_00152",
    #"output_00169",
    #"output_00191",
    #"output_00214",
    #"output_00237"
    ]
files.sort()
print ("Files to process %d"%len(files))

np.set_printoptions(precision=3,suppress=True)


# Loop and load data
# If we wanted to just loop over filename - "for file in files:"
for indx, file in enumerate(files):
    gc.collect()
    print ("Processing file %s" % file)
    s = pynbody.load("./"+file)
    s['pos'] -= 0.5
    s.physical_units();
    
    z = 1/s.properties['a']-1
    print ("Redshift = %.2lf" % z)
    boxsizestring = "%.2lf kpc" % s.properties['boxsize'].in_units('kpc')
    print ("Boxsize @ this redshift %s"%boxsizestring)
    print ("aexp %.3lf" % s.properties['a'])
    bs = float(s.properties['boxsize'])

    z_crit = 2.0e-7 # Mass fraction
    print ("Normalize Z data...")
    s.g['metal'][s.g['metal']<2e-7]    = 1e-10 # Since we divide by Z below, don't use 0.0
    s.g['pgf'][s.g['pgf']>(1.0-1e-10)]  = 1.0
    s.g['pzf'][s.g['pzf']<2e-7]        = 1e-10
    s.g['zsolar'] = s.g['metal'] * 50.0         # Solar units
    s.g['pzsolar'] = s.g['pzf'] * 50.0          # Solar units

    s.s['metal'][s.s['metal']<2e-7]    = 1e-10
    s.s['ppf'][s.s['ppf']>(1.0-1e-10)]  = 1.0
    s.s['pzf'][s.s['pzf']<2e-7]        = 1e-10
    s.s['zsolar'] = s.s['metal'] * 50.0         # Solar units
    s.s['pzsolar'] = s.s['pzf'] * 50.0          # Solar units

    # Fix birth times
    print (s.derivable_keys()[1:5])
    pynbody.analysis.ramses_util.get_tform(s); # Don't need path...

    sbox = 40.0 / (1.0 + z) * 0.71 # 40 kpc comoving box
    print ("40 kpc comoving Small box @ z=%.2lf is %.2lf"%(z,sbox))

    step = int(len(s.s)/3) # Note integer division
    if step == 0: step = 1
    for i in range(0,len(s.s),step):
        x,y,zz = - s.s['pos'][i] # get a star
        print ("Star offset [%.2lf %.2lf %.2lf]"%(x,y,zz))
        print ("Boxsize @ this z ", bs)
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
        print ("Plotsize: ",smallbox)

        # ##########################################################
        # Save the index and coordinates of the center of our image...
#        np.savetxt("z%05.2lf_SpCoord_%i.txt"%(z,i), -1.0*np.array(coords),comments='') 
        # ##########################################################

        rx,ry,rz = -1 * np.array([x,y,zz]) # For filters, we need the actual coords, not translation coords.
        print ("Cuboid center [%.2lf %.2lf %.2lf]"%(rx,ry,rz))
        # Note that the thickness in z is only 1/2 that of the other dimensions
        snBlast = s[pynbody.filt.Cuboid(str((rx-sbox/2.0)) + " kpc", str((ry-sbox/2.0)) + " kpc",str((rz-sbox/4.0)) + " kpc",
                                        str((rx+sbox/2.0)) + " kpc", str((ry+sbox/2.0)) + " kpc",str((rz+sbox/4.0)) + " kpc")]
        
        outputImages(z,i,snBlast)

    del s
    gc.collect()
    gc.collect()

print ("Done")
