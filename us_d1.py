##This section is to load in the modules that we need 
from netCDF4 import Dataset #This helps us open up and mess with the data files.
import matplotlib #This is a matplot library so that we can later be able to plot the data.
#matplotlib.use('agg') 
import matplotlib.pyplot as plt #This is where we can plot the data.
import matplotlib.colors as colors #This gives us access to the matplot color library.
import glob #This lets us dig through the directory to find the data we need.
import os  #This lets us do stuff inside the file we access.
import numpy as np #This lets us do mathematical formulas.
from matplotlib.cm import get_cmap #This brings in the option to use a colormap for matplotlib.
import cartopy.crs as crs #This allows us to use mapping features so we can have the US map on our plot.
from cartopy.feature import NaturalEarthFeature #This allows us to bring in country borders and stuff.
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords) #This allows us to use the wrf library to smooth our data.
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator,LongitudeLocator) #This lets us plot the coordinates on the plot.
import numpy.ma as ma # This allows us to mask out values in the plot that we do not want.
from metpy.plots import USCOUNTIES 
import metpy.calc as mpcalc


################### Function to truncate color map ###################
def truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100):
    '''truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100)'''    
    cmapIn = plt.get_cmap(cmapIn)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    return new_cmap

cmap_mod = truncate_colormap(minval=.2, maxval=1.0)  # calls function to truncate colormap 

# Read in the file(s)
i = 0
#filelist = glob.glob(os.path.join('/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023041800/', 'wrfout_d01_2022-*:00:00'))
filelist = ['/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023041800/wrfout_d01_2023-04-18_01:00:00']
print(filelist)
print("retrieved all wrfout files")
for filename in sorted(filelist): 
    if i == 0:
        print(filename)

#We are now starting with the first iteration (0=1 in code).
#We are digging into the certain directory with that date at any time and having it
#equal to filelist.
#We then begin to iterate each hour.

        ncfile = Dataset(filename)
        fileid = Dataset(filename, mode = 'r', format='cdf')
        
        
        # We opened the NetCDF above and now we are digging into the file to find the variables needed.
        # getvar is finding us the variable under its "name." We then define the array with [0,:,:].
        # We then print the array shape of the variable to make for sure that the arrays are similar.
                
        URaw = getvar(ncfile, "U10")
        U_Wind = fileid.variables["U10"][0,:,:]
        #print(np.shape(U_Wind))
        
        VRaw = getvar(ncfile, "V10")
        V_Wind = fileid.variables["V10"][0,:,:]
        
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
          
        CH4 = (CH4_BIO[0,:,:] + CH4_ANT[0,:,:] - CH4_BCK[0,:,:] + (CH4_TST[0,:,:] - CH4_BCK[0,:,:]) + (CO2_TST[0,:,:] - CH4_BCK[0,:,:]))
        
    
        # Get the latitude and longitude points#
        lats, lons = latlon_coords(CH4_Ant_Raw)
                      
        # Get the cartopy mapping object
        cart_proj =  get_cartopy(CH4_Ant_Raw)

        # Create a figure
        fig = plt.figure(figsize=(5,4))
        
        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
        
        # This is the highest and lowest values of CH4 we want to see. Anything above or below the values
        # will be masked out.
        
        MAXCHVAL = 2.10
        MINCHVAL = 1.85
        
        #Running iterations to plot out the data. We also have the array of numbers for our colorbar.
        
        ispeed = 1
        
        if ispeed: 
            ret = ax.projection.transform_points(crs.PlateCarree(), np.array(lons),
                np.array(lats)) # This method only accept ndarray!
            xx = ret[..., 0]
            yy = ret[..., 1]
        
        # If you want full color plot comment out the lines with cropped in it.            
         
            #cropped = (slice(65, 120, None), slice(150, 250, None))
            
            CH4[np.where(CH4 >= MAXCHVAL)] = ma.masked
            CH4[np.where(CH4 <= MINCHVAL)] = ma.masked
            
            #m = plt.contourf(xx[cropped], yy[cropped], to_np(CH4[cropped]),levels = 50,
                             #cmap=cmap_mod, extend ='both')            
            m = plt.contourf(xx, yy, to_np(CH4),levels = 50,
                              cmap=cmap_mod, extend ='both')
             
        #This is making a blue star on the map where El Reno and Pampa is so that
        #we know where they are.
        
 		# Add a star for the city of El Reno                    
        name =['El Reno']
        lat1 = [35.54122]
        lon1 = [-97.95494]
        plt.plot(lon1, lat1, 'm*', markersize = 10,label = name, transform=crs.PlateCarree())    

		# Add a title
        str_xhu1 = b''.join(fileid.variables['Times'][0,:]).decode()
        str_xhu2 = "CH4 @"
        str_xhu = str_xhu2+str_xhu1
        plt.title(str_xhu,fontsize=8)
        
        # This is adding a color bar and determining its shape and location on the plot.
        
        cbaxes = fig.add_axes([0.08, 0.15, 0.86, 0.04]) 
        
        cb=plt.colorbar(cax=cbaxes, orientation='horizontal',drawedges=True, shrink = .75)
        cb.ax.tick_params(labelsize=6,direction="in")
        
        cb.set_label("CH4 Concentration (ppm)",y=-100,fontsize=6, rotation=0)
                       
        #Download and add the states and coastlines
        states = NaturalEarthFeature(category="cultural", scale="50m",
                                      facecolor="none",
                                      name="admin_1_states_provinces_lines")
        ax.add_feature(states, linewidth=.5, edgecolor="black")
        country_borders = NaturalEarthFeature(category="cultural", scale="50m",
                                      facecolor="none",
                                      name="admin_0_boundary_lines_land")
        ax.add_feature(country_borders, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)
        
        # Set the map bounds
        ax.set_xlim(cartopy_xlim(CH4_Ant_Raw))
        ax.set_ylim(cartopy_ylim(CH4_Ant_Raw))    
        
        
		# Add wind vectors and a key
        skip = (slice(None, None, 15), slice(None, None, 15))
        
        Q=ax.quiver(XLon[skip], XLat[skip], U_Wind[skip], V_Wind[skip], transform = crs.PlateCarree())
        
        qk = plt.quiverkey(Q, .95, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', labelsep = 0.05, coordinates = 'axes')
        
		# Retrieve lat/lon points, add labels and formating 
        gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True, dms=True, x_inline=False, y_inline=False)
        gl.top_labels = False
        gl.left_labels = True
        gl.bottom_labels = True
        gl.right_labels = True
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = LongitudeLocator()
        gl.ylocator = LatitudeLocator()
        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter()
        
        #This will print out the figure name and add on what time it is. 
        #The plt.save line will also be saving the figure in the
        #set directory. It will then close it out and stop the iteration of
        #making plots once out of times.
        
        figname='wrfout_d01_US_CH4_BIO+ANT+Wet+agwaste_python'+str(i)+'.png'
        print(figname)
        plt.show()
        #plt.savefig((figname),dpi=300,format='png', bbox_inches = 'tight')
        plt.clf()
        plt.close()

