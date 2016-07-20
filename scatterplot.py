# coding: utf-8
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_ylim([-250,250])
ax1.set_xlim([-250,250])
ax1.scatter(locs[:,0],locs[:,1],s=mass/40.0,c=tag, cmap=cmap,norm=norm)
ax1.scatter(locs[:,0],locs[:,1],s=mass/40.0,c=Z, cmap=cmap,norm=norm)
# create a second axes for the colorbar
ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%.2lf')
ax.set_title('Well defined discrete colors')
ax1.set_title('Well defined discrete colors')
ax2.set_ylabel('Very custom cbar [-]', size=12)
plt.show()
