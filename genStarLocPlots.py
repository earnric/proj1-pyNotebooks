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
import sys, os, glob, pickle, gc
import math
#import mmap

mpl.rcParams['figure.figsize'] = (20,20)
mpl.rcParams['font.size'] = 18

## ## ## Setup colormap ## ## ##
# define the colormap
cmap = plt.cm.jet
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
cmaplist[0] = (.5,.5,.5,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

# define the bins and normalize
bounds = np.linspace(-10,0,11)
norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)
## ## ##

# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
locFiles  = [f for f in os.listdir('.') if re.match(r'spLoc_[0-9]*.*.txt', f)]
massFiles = [f for f in os.listdir('.') if re.match(r'spMass_[0-9]*.*.txt', f)]
ZFiles    = [f for f in os.listdir('.') if re.match(r'spZ_[0-9]*.*.txt', f)]
ppfFiles  = [f for f in os.listdir('.') if re.match(r'spPPF_[0-9]*.*.txt', f)]
ppzFiles  = [f for f in os.listdir('.') if re.match(r'spPZF_[0-9]*.*.txt', f)]

# Sort the files
locFiles.sort();massFiles.sort();ZFiles.sort();ppfFiles.sort();ppzFiles.sort()

zs  = np.loadtxt("zKeysForSPfiles.txt",skiprows=1)

if len(locFiles) != len(zs):
    print "Diff # of files to process than I have redshifts"
    os._exit(1)

dotNorm = 3.0
i = 0
for locF,massF,ZF,ppfF,pzfF in zip(locFiles,massFiles,ZFiles,ppfFiles,ppzFiles):
    z    = zs[i][0]
    bs   = zs[i][1]
    sbox = 5. / (1.0 + z) / 0.71 # Create a box that's 20kpc comoving
    print "z=%.3lf"%z
    spPosExp = r'z%.2lf_SpCoord_[0-9]*.txt' %z # These are the x,y,z locations of interest from the genStarPlots6.py
    poisFiles  = [f for f in os.listdir('.') if re.match(spPosExp, f)]
    poisFiles.sort()

    for indx, poi in enumerate(poisFiles):
        print "Star locs: ", locF
        print "Coord file: ", poi
        locs = np.loadtxt(locF,skiprows=1)
        mass = np.loadtxt(massF,skiprows=1)
        Z    = np.loadtxt(ZF,skiprows=1)
        ppf  = np.loadtxt(ppfF,skiprows=1)
        pzf  = np.loadtxt(pzfF,skiprows=1)
        x,y,zz = np.loadtxt(poi) # star particle location of interest
        xo,yo,zo = x,y,zz
        print "Star offset [%.2lf %.2lf %.2lf]"%(x,y,zz)
        print "Boxsize ", bs
        # Ensure we center such that we can depict a 40 kpc box
        if (abs(x) > bs/2.0 - sbox/2.0):
            x = np.sign(x) * (bs/2.0 - sbox/2.0)
        if (abs(y) > bs/2.0 - sbox/2.0):
            y = np.sign(y) * (bs/2.0 - sbox/2.0)
        if (abs(zz) > bs/2.0 - sbox/2.0):
            zz = np.sign(zz) * (bs/2.0 - sbox/2.0)
        print "Star offset adjusted for boxsize [%.2lf %.2lf %.2lf]"%(x,y,zz)
            
        xmin = x - sbox/2.0
        xmax = x + sbox/2.0
        ymin = y - sbox/2.0
        ymax = y + sbox/2.0
        zmin = zz - sbox/2.0 
        zmax = zz + sbox/2.0
        print "Box range x(%.2lf,%.2lf), y(%.2lf,%.2lf), z(%.2lf,%.2lf)"%(xmin,xmax,ymin,ymax,zmin,zmax)
        xcond = ((locs[:,0] >= xmin) & (locs[:,0] <= xmax))
        ycond = ((locs[:,1] >= ymin) & (locs[:,1] <= ymax))
        zcond = ((locs[:,2] >= zmin) & (locs[:,2] <= zmax))

        locsFilt = locs[xcond & ycond & zcond] # filter: just get stars in plot region
        if len(locsFilt) == 0: continue
        massFilt = mass[xcond & ycond & zcond] # filter: just get stars in plot region
        Zfilt    = Z[xcond & ycond & zcond] # filter: just get stars in plot region
        ppfFilt  = ppf[xcond & ycond & zcond] # filter: just get stars in plot region
        pzfFilt  = pzf[xcond & ycond & zcond] # filter: just get stars in plot region

        # Normalize locations
        locsFilt = locsFilt - np.array([x,y,zz])

        print "Z values all ", len(Zfilt)
        rng1 = (Zfilt < 2.e-7)
        rng2 = ((Zfilt >= 2.e-7) & (Zfilt < 1.e-4))
        rng3 = ((Zfilt >= 1.e-4) & (Zfilt < 1.e-2))
        rng4 = (Zfilt >= 1.e-2)
        z1=np.log10(Zfilt[rng1])
        z2=np.log10(Zfilt[rng2])
        z3=np.log10(Zfilt[rng3])
        z4=np.log10(Zfilt[rng4])

        # Create the starlocation Z plot
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10,10))
        # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
        plt.setp([a.get_xticklabels() for a in [ax1,ax2]], visible=False)
        plt.setp([a.get_yticklabels() for a in [ax2,ax4]], visible=False)
        plt.setp([a.set_xlim([-sbox/2.0,sbox/2.0]) for a in [ax1,ax2,ax3,ax4]])
        plt.setp([a.set_ylim([-sbox/2.0,sbox/2.0]) for a in [ax2,ax4,ax3,ax4]])

        xcoord = locsFilt[:,0]; ycoord = locsFilt[:,1]
        ax1.scatter(xcoord[rng1], ycoord[rng1], s=massFilt[rng1]/dotNorm, c=z1, cmap=cmap,vmin=-10, vmax=0)
        ax2.scatter(xcoord[rng2], ycoord[rng2], s=massFilt[rng2]/dotNorm, c=z2, cmap=cmap,vmin=-10, vmax=0)
        ax3.scatter(xcoord[rng3], ycoord[rng3], s=massFilt[rng3]/dotNorm, c=z3, cmap=cmap,vmin=-10, vmax=0)
        ax4.scatter(xcoord[rng4], ycoord[rng4], s=massFilt[rng4]/dotNorm, c=z4, cmap=cmap,vmin=-10, vmax=0)
        # create a second axes for the colorbar
        ax5 = fig.add_axes([0.87, 0.1, 0.025, 0.75])
        cb = mpl.colorbar.ColorbarBase(ax5, cmap=cmap, norm=norm, spacing='proportional',
                            ticks=bounds, boundaries=bounds, format='%1i')
        fig.suptitle('Star Particles z=%.1lf\n[%.2lf,%.2lf,%.2lf]'%(z,xo,yo,zo),size=24)
        ax1.set_title('Z < 2e-7',size=20)
        ax2.set_title('2e-7 < Z < 1e-4',size=20)
        ax3.set_title('1e-4 < Z < 1e-2',size=20)
        ax4.set_title('1.e-2 < Z',size=20)
        
        if sbox > 1.0:
            ax3.set_xlabel('x kpc')
            ax4.set_xlabel('x kpc')
            ax1.set_ylabel('y kpc')
            ax3.set_ylabel('y kpc')
        else :
            ax3.set_xlabel('x pc')
            ax4.set_xlabel('x pc')
            ax1.set_ylabel('y pc')
            ax3.set_ylabel('y pc')
            
        ax5.set_ylabel('Log Z', size=18)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.85, top=0.85, wspace=.15, hspace=.15)
        plt.savefig("SP_Z_locs_z=%.1lf-%i.png"%(z,indx))
#        plt.show()
        print "Files saved: SP_Z_locs_z=%.1lf-%i.png"%(z,indx)
        print "\n"

    i=i+1