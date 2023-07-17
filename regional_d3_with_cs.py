# This section is to load in the modules that we need 
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
                 cartopy_ylim, latlon_coords,vertcross, CoordPair, xy_to_ll) #This allows us to use the wrf library to smooth our data.
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
filelist = glob.glob(os.path.join('/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023041800/', 'wrfout_d03_2023-*:00:00'))
#filelist = ['/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023061100/wrfout_d03_2023-06-11_01:00:00']
#filelist = ['/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023060500/wrfout_d03_2023-06-05_02:00:00']
print(filelist)
print("retrieved all wrfout files")
for filename in sorted(filelist):
    if i >-1:
        print(filename)

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

        ht = getvar(ncfile, "z")
        ht = ht[0:20,:,:]          
        
	    # Doing the math so that we can get the data to be accurate.     
          
        CH4 = (CH4_BIO[0,:,:] + CH4_ANT[0,:,:] - CH4_BCK[0,:,:] + (CH4_TST[0,:,:] - CH4_BCK[0,:,:]) + (CO2_TST[0,:,:] - CH4_BCK[0,:,:]))
        
        # Get the latitude and longitude points#
        lats, lons = latlon_coords(CH4_Ant_Raw)
                      
        # Get the cartopy mapping object
        cart_proj =  get_cartopy(CH4_Ant_Raw)
         
        cart_proj = crs.PlateCarree()
        
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
            
            m = plt.contourf(xx, yy, to_np(CH4),levels = 50,
                              cmap=cmap_mod, extend ='both')
        
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
        
        cb=plt.colorbar(cax=cbaxes, orientation='horizontal',drawedges=True)
        cb.ax.tick_params(labelsize=6,direction="in")
        
        cb.set_label("CH4 Concentration (ppm)",y=-100,fontsize=6, rotation=0)
                       
        
        #This section here will add county borders so that our map
        #looks complete.
        
        ax.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='black')
                
        # Set the map bounds
        ax.set_xlim(cartopy_xlim(CH4_Ant_Raw))
        ax.set_ylim(cartopy_ylim(CH4_Ant_Raw))    
        
        #This line here is making the plot smaller so we can look at the
        #map more closely.
        
        ax.set_extent([-98.3, -97.9, 35.46, 35.66], crs=crs.PlateCarree()) 
        
        skip = (slice(None, None, 12), slice(None, None, 12))
        
		# Add wind vectors and key
        Q=ax.quiver(XLon[skip], XLat[skip], U_Wind[skip], V_Wind[skip],scale = 200, width= .002, transform = crs.PlateCarree())
        qk = plt.quiverkey(Q, .95, 1.05, 4, r'$4 \frac{m}{s}$', labelpos='W', labelsep = 0.05, coordinates = 'axes')
        
		### CALCULATIONS TO PLOT CROSS SECTION LINE ###

        # The cross section is calculated with a pivot point and an angle. We can either use the mean wind direction of the domain, 
        #or choose our own angle. 
        
        # Define the pivot point for the cross section. Right now, its set for the El Reno site
        pivot_point = CoordPair(lat=35.5348573, lon=-98.0984316)
        # Define the angle to take the cross section with. 0 is NS and 90 is WE. OR, uncomment the line marked below
        angle = 0
        
        # Calculate average wind direction
        V_east = np.mean(wspd * np.sin(wdir * np.pi/180))
        V_north = np.mean(wspd * np.cos(wdir * np.pi/180))
        mean_WD = np.arctan2(V_east, V_north) * 180/np.pi
        unit_WD = mean_WD
		#angle = unit_WD #uncomment if using mean wind direction
        
        # Calculation below retains meta data
        CH4_2 = (CH4_Bio_Raw + CH4_Ant_Raw - CH4_Bck_Raw + (CH4_Tst_Raw - CH4_Bck_Raw) + (CO2_Tst_Raw - CH4_Bck_Raw))
        CH4_2 = CH4_2[0:20,:,:]
        
        # Calculate the cross section to over plot 
        cross_section = vertcross(CH4_2, ht,  wrfin=ncfile, pivot_point= CoordPair(lat=35.534573, lon=-98.0984316), 
                                  angle=70, latlon=True, meta=True) 
        
        # Plot the cross section line  
        coord_pairs = to_np(cross_section.coords["xy_loc"])
        all_x = [pair.x for pair in to_np(cross_section.coords["xy_loc"])]
        all_y = [pair.y for pair in to_np(cross_section.coords["xy_loc"])] 
        lat_lon = xy_to_ll(ncfile, all_x, all_y)
        x = lat_lon[1,:]
        y = lat_lon[0,:]
        ax.plot(x, y, color = 'red')
        
        # This is a manual way for you to be able to add the lat and lon points onto the plot.
        
        ax.set_xticks([-97.9, -98.0, -98.1, -98.2, -98.3])
        
        ax.set_yticks([35.46, 35.50, 35.54, 35.58, 35.62, 35.66])
        
        # This line allows us to put our points on the right side of our plot as well.
        
        secax = ax.secondary_yaxis(1.0)
        secax.set_yticks([35.46, 35.50, 35.54, 35.58, 35.62, 35.66])        
        
        # These next two lines allow us to change the font size of all the points.
        
        ax.tick_params(axis = 'both', labelsize = 8)
        
        secax.tick_params(labelsize = 8)
        
        # All three of these sections allow us to be able to put words onto our lat and lon data so that we can 
        # know that it is degrees north and degrees west.
        
        labels = ax.get_xticks().tolist()
        new_labels = [str(l) + " \u00b0W" for l in labels]
        ax.set_xticklabels(new_labels)
        
        labels = ax.get_yticks().tolist()
        new_labels = [str(l) + " \u00b0N " for l in labels]
        ax.set_yticklabels(new_labels)
        
        labels = secax.get_yticks().tolist()
        new_labels = [str(l) + " \u00b0N " for l in labels]
        secax.set_yticklabels(new_labels)
        
        #This will print out the figure name and add on what time it is. 
        #The plt.save line will also be saving the figure in the
        #set directory. It will then close it out and stop the iteration of
        #making plots once out of times.
        
        figname='wrfout_d03_CH4_BIO+ANT+Wet+agwaste_python'+str(i)+'.png'
        print(figname)
        plt.savefig(os.path.join('/scratch/09283/rauth/wrf',figname),dpi=300,format='png', bbox_inches = 'tight',pad_inches = 0)
        plt.show()
        plt.clf()
        plt.close()
    i = i +1
