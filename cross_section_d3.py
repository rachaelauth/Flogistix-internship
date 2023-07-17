##This section is to load in the modules that we need 

import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import scipy.ndimage
import matplotlib.colors as colors #This gives us access to the matplot color library.
import matplotlib 
import glob
import os 
import cartopy.crs as crs #This allows us to use mapping features so we can have the US map on our plot.
import metpy.calc as mpcalc
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,cartopy_xlim,cartopy_ylim,
                 interpline, CoordPair)
from math import acos, degrees, pi, cos, sin
import numpy as np
import xarray as xr

################### Function to truncate color map ###################

def truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100):
    '''truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100)'''    
    cmapIn = plt.get_cmap(cmapIn)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    return new_cmap

# Read in the file(s)
i = 0
filelist = glob.glob(os.path.join('/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023041800', 'wrfout_d03_2023-*:00:00'))
#filelist = ['/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023060500/wrfout_d03_2023-06-05_02:00:00']
print(filelist)
print("retrieved all wrfout files")
for filename in sorted(filelist): 
    if i == 0:
        print(filename)
        # Open the NetCDF file
        ncfile = Dataset(filename)
        fileid = Dataset(filename, mode = 'r', format='cdf')

        # Get the WRF variables
        ht = getvar(ncfile, "z")
        ht = ht[0:20,:,:]
        ter = getvar(ncfile, "ter")
        
        w = getvar(ncfile, "wa")
        w = w[0:20,:,:]
        
        u = getvar(ncfile, "ua")
        u = u[0:20,:,:]
        
        v = getvar(ncfile, "va")
        v = v[0:20,:,:]
        
        spddir = getvar(ncfile,"wspd_wdir")
        wspd = spddir[0,:]
        wspd = wspd[0:20,:,:]
        wdir = spddir[1,:]
        wdir = wdir[0,:,:]
        
        LonRaw = getvar(ncfile, "XLONG")
        XLon = fileid.variables["XLONG"][0,:,:]
        
        LatRaw = getvar(ncfile, "XLAT")
        XLat = fileid.variables["XLAT"][0,:,:]
        
        CH4_Ant_Raw = getvar(ncfile, "CH4_ANT")
        CH4_ANT = fileid.variables['CH4_ANT'][0,:,:]
        
        CH4_Bio_Raw = getvar(ncfile, "CH4_BIO")
        CH4_BIO = fileid.variables['CH4_BIO'][0,:,:]
        
        CH4_Bck_Raw = getvar(ncfile, "CH4_BCK")
        CH4_BCK = fileid.variables['CH4_BCK'][0,:,:]
        
        CH4_Tst_Raw = getvar(ncfile, "CH4_TST")
        CH4_TST = fileid.variables['CH4_TST'][0,:,:]
        
        CO2_Tst_Raw = getvar(ncfile, "CO2_TST")
        CO2_TST = fileid.variables['CO2_TST'][0,:,:]
          
        # Doing the math so that we can get the data to be accurate.     
        CH4_2 = (CH4_Bio_Raw + CH4_Ant_Raw - CH4_Bck_Raw + (CH4_Tst_Raw - CH4_Bck_Raw) + (CO2_Tst_Raw - CH4_Bck_Raw))
        CH4_2 = CH4_2[0:20,:,:]
        
        # The cross section is calculated with a pivot point and an angle. We can either use the mean wind direction of the domain, 
        #or choose our own angle. 
        
        # Define the pivot point for the cross section. Right now, its set for the El Reno site
        pivot_point = CoordPair(lat=35.5348573, lon=-98.0984316)
        
        # Define the angle to take the cross section with. 0 is NS and 90 is WE. OR, uncomment the line marked below
        angle = 0 
        
        # Calculate average wind direction
        V_east = np.mean(wspd * np.sin(wdir * np.pi/180))
        V_north = np.mean(wspd * np.cos(wdir * np.pi/180))
        unit_WD = np.arctan2(V_east, V_north) * 180/np.pi
        #angle = unit_WD #uncomment if using mean wind direction

        
        # Take the CH4 cross section
        cross_section = vertcross(CH4_2, ht,  wrfin=ncfile, pivot_point= pivot_point, 
                                  angle=angle, latlon=True, meta=True) 
        
        # Take the terrain cross section
        ter_line = interpline(ter, wrfin=ncfile, pivot_point= pivot_point, 
                                  angle=angle, latlon=True, meta=True) 
        
        # Make a copy to work with - we use this to fill in the gap between the terrain and the start of the data
        cross_section_filled = np.ma.copy(to_np(cross_section))
        
        # For each cross section column, find the first index with non-missing
        # values and copy these to the missing elements below.
        for l in range(cross_section_filled.shape[-1]):
            column_vals = cross_section_filled[:,l]
        # Let's find the lowest index that isn't filled. The nonzero function
        # finds all unmasked values greater than 0. Since 0 is a valid value
        # for CH4, let's change that threshold to be -200 instead.
            first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
            cross_section_filled[0:first_idx, l] = cross_section_filled[first_idx, l]
            
        # Take the wind cross sections for the wind vectors
        total_cross = vertcross(wspd, ht, wrfin=ncfile,
                            pivot_point= pivot_point, 
                            angle=angle, latlon=True, meta=True)  
        
        w_cross = vertcross(w, ht, wrfin=ncfile,
                            pivot_point= pivot_point, 
                            angle=angle, latlon=True, meta=True)   
        
        
        ##### PLOTING THE CROSS SECTION CURTAIN PLOT #####
        
        # Get the cartopy projection object
        cart_proj = get_cartopy(w)
        
        cmap_mod = truncate_colormap(minval=.2, maxval=1.0)  # calls function to truncate colormap
        
        # Create the interval for color bar
        cmin = CH4_2.min()
        cmax = CH4_2.max()
        cint = 0.00225
        clevs = np.round(np.arange(cmin,cmax,cint),3)
        
        # Create the figure
        fig = pyplot.figure(figsize=(8,6))
        ax_cross = pyplot.axes()
        
        # Make the cross section plot for CH4
        xs = np.arange(0, cross_section.shape[-1], 1)
        ys = to_np(cross_section.coords["vertical"])
        
        # Plot the cross section
        w_contours = ax_cross.contourf(xs, ys, to_np(cross_section_filled),
                                  cmap=cmap_mod,levels= clevs, extend="both")
        # Add the color bar
        cbar = fig.colorbar(w_contours, ax=ax_cross)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label('CH4 (ppmv)', rotation=-270, fontsize=12)
        
        ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),
                                facecolor="black")
        
        # Quiver for the wind direction - rotate wind direction because vertcross only handles magnitude 
        Q = ax_cross.quiver(xs[::7], ys[::7],
                          (to_np(total_cross[::7, ::7])), to_np(w_cross[::7, ::7]))
        qk = plt.quiverkey(Q, .95, 1.025, 4, r'$4 \frac{m}{s}$', labelpos='W', labelsep = 0.075, coordinates = 'axes')
        
        # Set the x-ticks to use latitude and longitude labels
        coord_pairs = to_np(cross_section.coords["xy_loc"])
        x_ticks = np.arange(coord_pairs.shape[0])
        x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]
        
        # Set the desired number of x ticks below
        num_ticks = 5
        thin = int((len(x_ticks) / num_ticks) + .5)
        ax_cross.set_xticks(x_ticks[::thin])
        ax_cross.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)
        
        # Set the x-axis and  y-axis labels
        ax_cross.set_xlabel("Latitude, Longitude", fontsize=12)
        ax_cross.set_ylabel("Altitude (m)", fontsize=12)
        ax_cross.set_ylim(ymax = 2000)
        
        # Add the title
        str_xhu1 = b''.join(ncfile.variables['Times'][0,:]).decode()
        str_xhu2 = "Vertical Cross-Section @"
        str_xhu = str_xhu2+str_xhu1
        plt.title(str_xhu,fontsize=10)
        
        # Show and save the figure (change file path to desired directory)
        figname='wrfout_d03_CH4_BIO+ANT+Wet+agwaste_python_cross_section'+str(i)+'.png'
        print(figname)
        plt.savefig(os.path.join('/scratch/09283/rauth/wrf',figname),dpi=300,format='png', bbox_inches = 'tight')
        plt.show()
        plt.clf()
        plt.close()
    i = i +1