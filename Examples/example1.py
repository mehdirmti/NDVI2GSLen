# importing librarries and packages
import numpy as np
import netCDF4
import pandas as pd

import sys
sys.path.append('../Codes/')

# importing the personal module of LFD_NDVI_Method
from algorithm import LFD_NDVI

# define a function to apply LFD method over different pixels
def LFD_Over_Pixels(y, date):
    # Define an internal function to apply LFD over different years
    def LFD_Over_Years(YR):

        # indexation to seperate data for each specific year
        ind = Year == YR
        x_in = x[ind]
        y_in = y[ind]

        # call LFD method to do calculation for each individual year
        OG, OD, OG_ndviC, OD_ndviC, peak_Timing, PearsonCorrelation = LFD_NDVI(x_in, y_in, YR, 'off')
        return OG, OD, OG_ndviC, OD_ndviC, peak_Timing, PearsonCorrelation
    # end of function

    ###### Do the calculations ####################################################################
    # get corresponding day of year for dates
    x = np.array(date.day_of_year)

    # get corresponding years of dates
    Year = np.array(date.year)

    # find unique years during the examined period
    Unique_year = np.unique(np.array(date.year))
    Unique_year = Unique_year.reshape((1, len(Unique_year)))

    # Call internal Function to apply LFD method over different years
    OG, OD, OG_ndviC, OD_ndviC, peak_Timing, PearsonCorrelation = np.apply_along_axis(LFD_Over_Years, 0, Unique_year)
    return OG, OD, OG_ndviC, OD_ndviC, peak_Timing, PearsonCorrelation
# end of function

####################### calculations: MODIS data #############################
nc = netCDF4.Dataset('../Data/ndvi-Euro-MODIS.nc')
ndvi = nc.variables['ndvi']
lat = np.array(nc.variables['lat'])
lon = np.array(nc.variables['lon'])

# Define the date vector in Monthly resolution
date = pd.period_range(start="2001-01-01",end="2020-12-31",freq="M")

# Call LFD method and apply it for axis = 0 (time axis)
OG, OD, OG_ndviC, OD_ndviC, peak_Timing, PearsonCorrelation = np.apply_along_axis(LFD_Over_Pixels, 0, ndvi, date)

# Store the outputs for further use
np.savez('LFD_Results', 
        OG = OG, OD = OD, 
        OG_ndviC = OG_ndviC, 
        OD_ndviC = OD_ndviC, 
        peak_Timing=peak_Timing, 
        PearsonCorrelation = PearsonCorrelation, 
        lat=lat, lon=lon, year=np.unique(date.year))