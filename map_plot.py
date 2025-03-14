#_._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._..
#
#* File Name : map_plot.py
#
#* Purpose : plotting geospatial data (either globally or given specific extent)
#
#* Creation Date : Sat 29 Jul 2023
#
#* Last Modified : Mon 20 Jan 2025 10:49:43 (CET)
#
#* Created By : Aqeel Piracha 	Github. apiracha1
#								gmail. piracha.aqeel1
#
#
#_._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._..*/
#Imports
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cf
import sys
from warnings import filterwarnings
filterwarnings('ignore')

"""
Extent is;
minimum latitude, maximum longitude, minimum longitude, maximum latitude
"""

#setting axis projection
def map_plot(ds, vmax='0', vmin='0', cmap='turbo', cbar=0, title='', savefig=0, imname='figure',
             isglobal=1, central_lon=0, extent=[0,0,0,0], isrect=1, rect=[0,0,0,0]):
#checking for correct dimensions
    if len(ds.shape) > 2 and ds.shape[0] > 1:
        sys.exit('MAKE SURE 3RD DIMENSION IS SHAPE 1')
#setting axis projection
    proj = ccrs.Robinson(central_longitude=central_lon)
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(projection=proj)
#regional plot if chosen
    if isglobal == 0:
        ax.set_extent(extent, ccrs.PlateCarree())
    ax.coastlines(color='.5')
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
#making sure current figure and axis is focused
    plt.sca(ax)
#plotting
    if cbar == 1:
        ds.plot(vmax=vmax,vmin=vmin, transform=ccrs.PlateCarree(), cmap=cmap, ax=ax, 
                cbar_kwargs={'fraction': 0.046, 
                             'pad': 0.06})
    if cbar == 0:
        ds.plot(vmax=vmax,vmin=vmin, transform=ccrs.PlateCarree(), cmap=cmap, ax=ax, 
                             add_colorbar=False)
        
    if isrect == 1:
        ax.add_patch(mpatches.Rectangle(xy=[rect[1], rect[0]], width=rect[2], height=rect[3],
                                    facecolor='none',
                                    edgecolor='k',
                                    linewidth=3,
                                    transform=ccrs.PlateCarree())
                 )
    plt.title(title, fontdict={'fontsize': 15,'fontweight': 'bold'})
#saving if chosen
    if savefig == 1:
        plt.savefig(imname + '.jpg', dpi=600)
        plt.close()

