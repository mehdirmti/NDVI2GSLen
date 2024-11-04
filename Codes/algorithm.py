"""
Created on Wed Sep 28 08:00:00 2022
Revised on Tuesday May 23 10:15:00 2023

@author: Mehdi Rahmati
Email: mehdirmti@gmail.com, m.rahmati@fz-juelich.de

Description:
This script provides functions to implement the logistic function derivative (LFD) NDVI method.
The method uses the seasonality of the NDVI curve to determine the onset of greening and vegetation dormancy. 

More detailed information can be found in the article published at the following link:
https://www.XXX.XXX

Inputs and Outputs:
The main function is called: LFD_NDVI(DoY, ndvi, year, plotting).

    @Inputs:
    --------
    DoY: is a numpy.ndarray with a shape of (n,) indicating the days for which NDVI data is provided
    ndvi: is a numpy.ndarray with the shape of (n,) containing NDVi values for each day represented in DoY
    year: is a scaler int specifying the year for which NDVI data is provided
    ploting: is a string. If "on" or "On" or "ON" is specified, the results will be displayed in the graph

    @outputs:
    ---------
    OG: is a numpy.float64 indicating the day of the year when greening starts.
    OD: is a numpy.float64 indicating the day of the year when vegetation dormancy begins.
    OG_ndviC: is a numpy.float64 indicating the critical NDVI values on the days of revegetation.
    OD_ndviC: is a numpy.float64 that provides the critical NDVI values on the days of vegetation dormancy.
    Peaking_time: is a numpy.float64 that gives the day of the year when revegetation peaks.
    R: is a numpy.float64 indicating the Pearson correlation between the provided NDVI data and the fitted logical function.

To be cited as: 
Rahmati et al., 2023.The continuous increase in evaporative demand shortened the growing season of European ecosystems 
in the last decade. Comunications Earth & Environement, XX, XX-XX. 
"""
import numpy as np
import matplotlib.pylab as plt
from scipy import optimize
from scipy.signal import find_peaks
from scipy.optimize import root_scalar
import warnings
np.warnings = warnings
from scipy.interpolate import interp1d
import calendar
from matplotlib import ticker
import pandas as pd

#warnings.simplefilter('ignore', np.RankWarning)
warnings.filterwarnings('ignore')

# determine the maximum curvature point(s) of cumulative curve
def curvature_calc(y):
    # y is a numpy.ndarray with the shape of (n,) containing estimated NDVi values for all individual days of year
    m = len(y)
    dt = 1 # delta t is assumed to be 1 as y provides daily data

    # Allocation of containers for the first and second derivatives of y
    yp = np.full((m), fill_value=np.nan)
    ypp = np.full((m), fill_value=np.nan)

    # calculate the first derivative of the y
    yp[2:m-1] = (y[3:m] - y[1:m-2])/(2*dt)  
    
    # Second derivative of the cy
    ypp[2:m-1] = (y[3:m] + y[1:m-2] - 2*y[2:m-1]) / (dt**2) 

    # time dependent curvature of y
    k = np.abs(ypp)/(1 + yp**2)**1.5 

    # find maximum curvature points (we will have two curvature points)
    pks = find_peaks(k)[0]
    vls = find_peaks(-1.*k)[0]
    

    return k, pks, vls

""" 
Define a function that applies Sequential Linear Approximation Method to determine 
the crtitical time below/above which a linear relationshop exist between t and y 
"""
def Sequential_Linear_Approximator(t, y, portion):

    """
    Inputs:
    t: is a numpy.ndarray with a shape of (n,) indicating the days for which NDVI data is provided
    y: is a numpy.ndarray with the shape of (n,) containing NDVi values for each day represented in t
    portion: determines in which part of year we look for linearity: "start_of_year" or "end_of_year"
    """

    # if we the linearity is overlooked at the end of the year, the data should be fliped before calculation
    if portion == "end_of_year":
        ind = np.flip(np.argsort(t))
        t = t[ind]
        y = y[ind]
    
    # Apply Sequential Linear Approximation Method to determine the crtitical time
    # we start with minimum data points of 3.

    # Allocation of container for the correlation coefiicient between y and y_hat
    RR = np.full((len(t)), fill_value=np.nan)
    for i in range(3,len(t)):
        
        # fit a linear line over data at intial steps
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                p = np.polyfit(t[0:i], y[0:i], 1)
            except np.RankWarning:
                print("not enough data for polyfit function [line 109]")
        
        # do predictions applying the fitted linear line 
        y_hat = []
        y_hat = np.polyval(p, t[:i])

        # Calculate correlation coefficient between y and y_prd
        rr = (np.corrcoef(y_hat, y[:i]))**2
        
        # Store R^2 for further comparisons
        RR[i] = rr[0,1]

        # Check if folloiwng criteria occurs to break the loop to save the time
        if RR[i] < 0.98 and RR[i-1] >= 0.98:
            break
    
    if RR[i] < 0.98:
        try:
            t_critical = np.interp(0.98, RR[~np.isnan(RR)], t[~np.isnan(RR)])
        except (RuntimeError, TypeError, NameError):
            t_critical = np.nan
    else:
        t_critical = t[i]

    return t_critical

# Define a function to determine the bisector line between two lines 
def Finding_Bisector_between_Two_Line(line1, line2, intersection_point):
    
    # compute angle between lines
    alpha1 = np.degrees(np.arctan(line1[0]))
    alpha2=np.degrees(np.arctan(line2[0]))
    alpha = abs(alpha1-alpha2)

    # compute bisector line of acute angle 
    slope_acute = np.tan((alpha1+alpha/2)*np.pi/180)
    intercept_acute = intersection_point[1] - slope_acute*intersection_point[0]

    # compute perpendicular line to bisector line of acute angle: obtuse angle
    slope_obtuse = np.tan((alpha1+alpha/2+90)*np.pi/180)
    intersept_obtuse = intersection_point[1] - slope_obtuse*intersection_point[0]

    return [slope_obtuse, intersept_obtuse]

# define a function that does the main job and determines the onsets of greening and dromancy    
def LFD_NDVI(DoY, ndvi, year, plotting):

    # find number of days in year of attention
    n_day = 366 if calendar.isleap(year) else 365

    # Set Xtick and Xticklabel for plotting
    date = pd.date_range(str(year)+'/01/01', str(year)+'/12/31', freq='M')
    date = np.insert(np.array(date.days_in_month), 0, 0)
    date = np.cumsum(date)
    date = date[:-1]
    XTick = (date)/(n_day)

    XtickLabe = []
    for i in range(1,13):
        XtickLabe.append(calendar.month_abbr[i])

    # introduce a array of time providing a full range of DoY for the year of attention 
    time = np.arange(1, n_day+1, 1)

    # Allocate space for outputs 
    OG = np.nan
    OD = np.nan
    OG_ndviC = np.nan
    OD_ndviC = np.nan
    Peaking_time = np.nan
    R = np.nan

    # Before progressing, check to see if ndvi does no contain only NaN, if it is the case, then return 
    if len(ndvi[np.isnan(ndvi)]) == len(ndvi):
        print('The provided NDVI data contains only NaN values. Please check and prevent the occurrence of NaN')
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    
    # plot Original NDVI data versus time
    if plotting =='on' or plotting == 'On' or plotting == 'ON':     

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)  

        # set some initial settings on the figure
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams.update({'font.size': 12})
        fig.set_figheight(8)
        fig.set_figwidth(14) 

        ax1.plot(DoY, ndvi)
        ax1.set_ylabel('NDVI [-]',**{'fontname':'Times New Roman', 'size':'12'})
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['left'].set_visible(True) 
        ax1.set_ylim(-0.05, 0.8)
        ax1.set_xlim(0,365)
        ax1.set_xticks(date)
        ax1.set_xticklabels(XtickLabe)
        ax1.text(-40, 0.75, 'a)',**{'fontname':'Times New Roman', 'size':'12'})
        
        # Set the tick labels font
        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontname('Times New Roman')
            label.set_fontsize(12)

    # Scaling of the data of both axes in 0 and 1. 
    # This is necessary to accurately determine the angle bisectors. 
    # For the y-axis, it also helps to remove the negative values.
    y = (ndvi - np.nanmin(ndvi))/(np.nanmax(ndvi)-np.nanmin(ndvi)) 
    x = (DoY - 1)/(n_day - 1)
    t = (time - 1)/(n_day - 1)

    # do cumulative sum for scaled NDVI
    y = np.nancumsum(y)
    
    # Re-scale cumualtive data of y in 0 and 1 once again
    y = (y - np.nanmin(y))/(np.nanmax(y)-np.nanmin(y)) 

    # define logit function
    def logit_func(x, a, b, MU, S):
        return a + b/(1 + np.exp(-(x-MU)/S))
    
    # Set initial values for parameters
    P0=(0, 1, 0.5, 0.1)

    # Screen out NaN data from y and x vectors before fitting
    ind = np.logical_and(~np.isnan(x), ~np.isnan(y))

    # check if y is empty, then return with NaN values for outputs
    if len(y[ind]) == 0:
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # fit logit function over data
    try:
        popt, pcov = optimize.curve_fit(logit_func, x[ind], y[ind], p0=P0, method="lm", maxfev=5000)
    # Sometime the optimum value is not reached after 5000 times iteration, so return NaN for outputs in this case and exit
    except: 
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R
    
    # Check Pearson Correlation between y and y_hat
    R = np.corrcoef(logit_func(x, *popt),y)[0,1]

    # Predict y for all possible days within year which is defined by t vector           
    y_hat = logit_func(t, *popt)

    # plot cumultive NDVI (original and fitted)
    if plotting =='on' or plotting == 'On' or plotting == 'ON':
        ax2.plot(x, y, 'o')
        ax2.plot(t, y_hat)
        ax2.set_ylabel('Re-scaled cumulative NDVI [-]',**{'fontname':'Times New Roman', 'size':'12'})
        ax2.set_ylim(-0.05,1.2)
        ax2.set_xlim(0,1)
        ax2.set_xticks(XTick)
        ax2.set_xticklabels(XtickLabe)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(True)
        ax2.spines['left'].set_visible(True)
        lg = ax2.legend(['Data', 'logistic function [R = {}]'.format('%.3f'%(R))], loc='lower right') 
        lg.get_frame().set_facecolor('none')
        lg.get_frame().set_edgecolor('none')
        ax2.text(-0.1, 1.1, 'b)',**{'fontname':'Times New Roman', 'size':'12'})

        # Set the tick labels font
        for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            label.set_fontname('Times New Roman')
            label.set_fontsize(12)

    # calculate curvature of y_hat curve for each day
    k, pks, vls = curvature_calc(y_hat)

    # plot time-dependent curvarture of data
    if plotting =='on' or plotting == 'On' or plotting == 'ON':
        ax3.plot(t, k)
        #ax3.set_title('Curvature')
        ax3.set_xlabel('DoY',**{'fontname':'Times New Roman', 'size':'12'})
        ax3.set_ylabel('Curvature [-]',**{'fontname':'Times New Roman', 'size':'12'})
        ax3.set_xlim(0,1)
        ax3.set_xticks(XTick)
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True) 
        formatter.set_powerlimits((-1,1)) 
        ax3.yaxis.set_major_formatter(formatter)
        ax3.set_xticklabels(XtickLabe)

        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['bottom'].set_visible(True)
        ax3.spines['left'].set_visible(True) 
        ax3.text(-0.1, 2.2*10**-4, 'c)',**{'fontname':'Times New Roman', 'size':'12'})

        # Set the tick labels font
        for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
            label.set_fontname('Times New Roman')
            label.set_fontsize(12)

    # Check to see if you have got desired numbers of peaks
    # Desired number of peak is 2: one for greening and one for dormancy 
    # If number of peaks are less than 2 or higher than 2, then return NaN for outputs
    if len(pks) != 2:
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # Check to see if you have got desired numbers of valleys, as well.
    # Desired number of valleys is 1 which occurs in peaking day of NDVI 
    # If number of valleys are less than 1 or higher than 1, then return NaN for outputs
    # the valley should occur between first and second peaks
    if np.logical_and(len(vls) == 1, np.logical_and(vls[0] > np.min(pks), vls[0] < np.max(pks))):
        Peaking_time = vls[0]
    else:
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # determine the time and cum scaled NDVI at points when the first and second peaks in curvature occures
    t_1k = t[pks[0]]
    t_2k = t[pks[1]]
    y_1k = y_hat[pks[0]]
    y_2k = y_hat[pks[1]]        

    # determine the criticale time at early days of year where NDVI vs time is linear        
    if len(t[t< t_1k]) > 3:
        t1_critical = Sequential_Linear_Approximator(t[t< t_1k], y_hat[t < t_1k], 'start_of_year')
    else:
        # If the above condtion is not the case, then we simply assume that the linear regression 
        # between first three data points define the linear part, so we consider the DoY=3 as critical time
        t1_critical = t[2] 
    
    # if you have got a NaN value for above mentioned critical time, 
    # then return NaN for all outputs and quite the calculations
    if np.isnan(t1_critical):
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # fit a linear line over data where t < t1_critical
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            p1 = np.polyfit(t[t <= t1_critical], y_hat[t <= t1_critical], 1)
        except np.RankWarning:
            print("not enough data for polyfit function [line 346]")    
    

    # fit a linear line over data falling between peaks of curvature 
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            p2 = np.polyfit((t_1k, (t_1k + t_2k)/2, t_2k), (y_1k, (y_1k + y_2k)/2, y_2k), 1)
        except np.RankWarning:
            print("not enough data for polyfit function [line 355]") 

    # determine the criticale time at late days of year where NDVI vs time is linear        
    if len(t[t> t_2k]) > 3:
        t2_critical = Sequential_Linear_Approximator(t[t> t_2k], y_hat[t > t_2k], 'end_of_year')
    else:
        t2_critical = t[len(t)-3]
    
    # fit a linear line over data at late time steps [t < t1_critical]
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            p3 = np.polyfit(t[t>= t2_critical], y_hat[t>= t2_critical], 1)
        except np.RankWarning:
            print("not enough data for polyfit function [line 369]") 
            
    # if you have got a NaN value for above mentioned critical time, 
    # then return NaN for all outputs and quite the calculations
    if np.isnan(t2_critical):
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # determine the intersection point of first and second lines
    intersection_1 = np.empty((2))
    intersection_1[0] =  (p1[1] - p2[1])/(p2[0] - p1[0])
    intersection_1[1] = np.polyval(p2, intersection_1[0])

    # determine the  of bisector line between lines 1 and 2
    p4 = Finding_Bisector_between_Two_Line(p1, p2, intersection_1)

    # determine a function to be solved to find the intersection of logistic curve and bisector line
    def fun1(t, a1, b1, c1, d1, e1, f1):
        return a1 + b1/(1 + np.exp(-(t-c1)/d1)) - e1*t - f1

    # Solve the above function to find the intersection
    try:
        sol = root_scalar(fun1, args=(*popt, *p4), method='toms748', bracket=[0, 1])
    except: # In the case no slution for the function, return NaN for outputs and quit the calculations
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # Check to see if solution is converged
    if sol.converged == True:
        OG = sol.root
    else: # otherwise return NaN for outputs and quit the calculations
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R
    
    # determine critical cumulative NDVI value at onset day of greening
    y1_critical = logit_func(OG, *popt)
    
    # determine the intersection point of third and second lines
    intersection_2 = np.empty((2))
    intersection_2[0] =  (p3[1] - p2[1])/(p2[0] - p3[0])
    intersection_2[1] = np.polyval(p3, intersection_2[0])

    # determine the bisector line between lines 3 and 2
    p5 = Finding_Bisector_between_Two_Line(p3, p2, intersection_2)
    try: 
        sol = root_scalar(fun1, args=(*popt, *p5), method='toms748', bracket=[0, 365])
    except: # In the case no slution for the function, return NaN for outputs and quit the calculations
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # Check to see if solution is converged   
    if sol.converged == True:
        OD = sol.root
    else: # otherwise return NaN for outputs and quit the calculations
        return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R

    # determine critical cumulative NDVI value at onset day of dromancy
    y2_critical = logit_func(OD, *popt)

    # plot all results together
    if plotting =='on' or plotting == 'On' or plotting == 'ON':
        ax4.plot(t, y_hat) # plot cum scaled ndvi
        ax4.plot(t[pks], y_hat[pks], 'o') # plot peak curvature points
        ax4.plot(t[t< t_1k], np.polyval(p1, t[t< t_1k]), '--');   # plot first linear line
        ax4.plot(np.linspace(t_1k*0.5,t_2k*1.2), np.polyval(p2, np.linspace(t_1k*0.5,t_2k*1.2)), '--')   # plot second linear line
        ax4.plot(t[t> t_2k*0.8], np.polyval(p3, t[t> t_2k*0.8]), '--') # plot third linear line
        ax4.plot(np.linspace(intersection_1[0]-0.04, intersection_1[0]+0.04), np.polyval(p4, np.linspace(intersection_1[0]-0.04, intersection_1[0]+0.04)), '--') # plot first bisector line
        ax4.plot(np.linspace(intersection_2[0]-0.04, intersection_2[0]+0.04), np.polyval(p5, np.linspace(intersection_2[0]-0.04, intersection_2[0]+0.04)), '--') # plot second bisector line
        ax4.plot(OG, y1_critical, '*', markerfacecolor='g', markeredgecolor='g')
        ax4.plot(OD, y2_critical, '*', markerfacecolor='y', markeredgecolor='y')
        ax4.set_ylim(-0.05,1.2)
        ax4.set_xlim(0,1)
        ax4.set_xticks(XTick)
        ax4.set_xticklabels(XtickLabe)
        ax4.set_xlabel('DoY',**{'fontname':'Times New Roman', 'size':'12'})
        ax4.set_ylabel('Cum scaled NDVI [-]',**{'fontname':'Times New Roman', 'size':'12'})
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['bottom'].set_visible(True)
        ax4.spines['left'].set_visible(True)
        lg = ax4.legend(['Data', 'Max curvature points', 'Winter dormant period', 'Active growth period', 'Fall dormant period', "Bisector 1", "Bisector 2", 'Onset of greening', 'Onset of dormancy'], loc='upper left') 
        lg.get_frame().set_facecolor('none')
        lg.get_frame().set_edgecolor('none')
        ax4.text(-0.1, 1.1, 'd)',**{'fontname':'Times New Roman', 'size':'12'})

        # Set the tick labels font
        for label in (ax4.get_xticklabels() + ax4.get_yticklabels()):
            label.set_fontname('Times New Roman')
            label.set_fontsize(12)

        #fig.savefig('E:/ET-SWC paper/#6 Final Check/Extended data Figure 1.jpg', dpi=600, transparent=False)   # save the figure to file

        plt.show()
    
    # Convert results into original range
    OG = np.round(OG * (n_day - 1) + 1)
    OD = np.round(OD * (n_day - 1) + 1)

    ndvi_interp = interp1d(DoY, ndvi)
    if OG >= DoY[0]:
        OG_ndviC = ndvi_interp(OG)
    if OD <= DoY[len(DoY)-1]:
        OD_ndviC = ndvi_interp(OD)
    return OG, OD, OG_ndviC, OD_ndviC, Peaking_time, R