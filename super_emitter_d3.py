from netCDF4 import Dataset
from wrf import to_np, getvar,  latlon_coords, CoordPair, ll_to_xy, xy_to_ll
import numpy as np

# I created a copy so I didn't break the original created with bare cp
my_filename = '/scratch/09283/rauth/wrf/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023041800/wrfchemi_00z_d03'
#this is being used for domain 3 coordinates
my_fdomain = "/scratch/09283/rauth/wrf/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2023041800/wrfout_d03_2023-04-18_01:00:00"
my_file = Dataset(my_filename, mode="r+")
fdomain = Dataset(my_fdomain, mode ="r")

E_CH4 = my_file.variables["E_CH4"][:,:,:,:]
loc = ll_to_xy(fdomain,35.5348573,-98.0984316)
loc.astype(int)
i = loc[0]
j = loc[1]
print(i,j)
v = 0
my_file.variables["E_CH4"][:,v,j,i] = 2000
my_file.close()

# Verify
my_file = Dataset(my_filename, mode="r")
E_CH4 = my_file.variables["E_CH4"][:,:,:,:]
print(E_CH4[:,0,loc[1],loc[0]])
my_file.close()