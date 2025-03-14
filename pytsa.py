#_._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._..
#
#* File Name : pytsa.py
#
#* Purpose : time series analysis toolbox (lnear trends, harmonic analysis,...) and statistics
#
#* Creation Date : Sun 04 Jun 2023
#
#* Last Modified : Thu 04 Jan 2024 10:03:16 (CET)
#
#* Created By : Aqeel Piracha 	Github. apiracha1
#								gmail. piracha.aqeel1
#
#
#_._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._..*/
##Loading libraries
import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pymannkendall as mk
from scipy.stats import norm
from warnings import filterwarnings

filterwarnings('ignore')

def linear_trend(ds, ds_name, fx, x, time_name, method='least_squares', factor=12, ext='', 
                 seasonal=True):

    print('computing linear trends for ' + ds[ds_name].attrs['long_name'])
    if method == 'mk':
# obtaining Theil-sen slope and intercept and other statistics related to significance
# via the mann-kendal and p-values assuming 3d dataset with dimension 0 = time, dimension 1 = lat,
# dimension 2 = lon
# [0] trend: This tells the trend-increasing, decreasing, or no trend.
# [1] h: True if the trend is present. False if no trend is present.
# [2] p: The p-value of the test.
# [3] z: The normalized test statistic.
# [4] Tau: Kendall Tau.
# [5] s: Mann-Kendal’s score
# [6] var_s: Variance S
# [7] slope: Theil-Sen estimator/slope
# [8] intercept: Intercept of Kendall-Theil Robust Line
## initialising matrices
        z = np.zeros((fx.shape[1],fx.shape[2]))
        Tau = 0
        var_s = 0
        trend = np.zeros((fx.shape[1],fx.shape[2]))
        s = np.zeros((fx.shape[1],fx.shape[2]))
        intercept = np.zeros((fx.shape[1],fx.shape[2]))
        pValue = np.zeros((fx.shape[1],fx.shape[2]))
        h = np.zeros((fx.shape[1],fx.shape[2]))
        xBar = np.nanmean(x, axis=0)
        xErr = x - xBar
        xErrSqr = np.power(xErr,2)

        if seasonal == True:
            for lat in np.arange(fx.shape[1]):
                print("\rprogress = {}%".format(int(lat/fx.shape[1]*100)), end=' ', flush='true')
                for lon in np.arange(fx.shape[2]):
                    if np.isnan(fx[0,lat,lon]):
                        z[lat,lon], s[lat,lon], h[lat,lon], pValue[lat,lon], trend[lat,lon], \
                                intercept[lat,lon] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                        continue
                    if sum(np.isnan(fx[:,lat,lon])) >= fx.shape[0]-3:
                        s[lat,lon], h[lat,lon], pValue[lat,lon], trend[lat,lon], \
                                intercept[lat,lon] = np.nan, np.nan, np.nan, np.nan, np.nan
                        continue
                    if sum(~np.isnan(fx[:,lat,lon])) > 12:
                        m = mk.seasonal_test(fx[:,lat,lon])
                    else:
                        m = mk.original_test(fx[:,lat,lon])
                    z[lat,lon], s[lat,lon], h[lat,lon], pValue[lat,lon], trend[lat,lon], \
                            intercept[lat,lon] = m[3], m[5], m[1], m[2], m[7], m[8]
        else:
            for lat in np.arange(fx.shape[1]):
                print("\rprogress = {}%".format(int(lat/fx.shape[1]*100)), end=' ', flush='true')
                for lon in np.arange(fx.shape[2]):
                    if np.isnan(fx[0,lat,lon]):
                        z[lat,lon], s[lat,lon], h[lat,lon], pValue[lat,lon], trend[lat,lon], \
                                intercept[lat,lon] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                        continue
                    if sum(np.isnan(fx[:,lat,lon])) >= fx.shape[0]-3:
                        s[lat,lon], h[lat,lon], pValue[lat,lon], trend[lat,lon], \
                                intercept[lat,lon] = np.nan, np.nan, np.nan, np.nan, np.nan
                        continue
                    m = mk.original_test(fx[:,lat,lon])
                    z[lat,lon], s[lat,lon], h[lat,lon], pValue[lat,lon], trend[lat,lon], \
                            intercept[lat,lon] = m[3], m[5], m[1], m[2], m[7], m[8]

        print("\n")
## slope is the actual trend line throughout time
        slope = intercept + trend * x
#calculating residual between the slope and data
        e = fx - slope
#variance of residual about trend line
        se_2 = (1 / (fx.shape[0] - 2)) * np.nansum(np.power(e,2), axis=0)
#calculating standard error of the trend
        sb_2 = se_2 / np.nansum(xErrSqr, axis=0)
#adding extra variables
        signif = abs(z) > 1.7
        ds["mk_s" + ext] = (['lat','lon'], signif)
        ds["mk_s" + ext].attrs['long_name'] = 'mann-kendall s score' + ext
        ds["mk_s" + ext].attrs['units'] = ''
        ds["trend_present" + ext] = (['lat','lon'], pValue)
        ds["trend_present" + ext].attrs['long_name'] = 'trend presence' + ext
        ds["trend_present" + ext].attrs['units'] = \
                'n/a'

    elif method == 'ls':
## simple ordinary least squared regression
## Calculating lineartrends through linear regression
        yBar = np.nanmean(fx, axis=0)
        yErr = fx - yBar
        xBar = np.nanmean(x, axis=0)
        xErr = x - xBar
        xErrSqr = np.power(xErr,2)
        corr = xErr * yErr
#getting trends
        numerator = np.nansum(corr, axis=0)
        denominator = np.nansum(xErrSqr, axis=0)
        trend = numerator / denominator
        intercept = yBar - trend * xBar
        slope = intercept + trend * x

## Analysing significance of trend
#calculating residual between the slope and data
        e = fx - slope
#variance of risidual about trend line
        se_2 = (1 / (fx.shape[0] - 2)) * np.nansum(np.power(e,2), axis=0)
#calculating standard error of the trend
        sb_2 = se_2 / np.nansum(xErrSqr, axis=0)
#t-value for significance
        t = trend/np.sqrt(sb_2)
#pValue according to normal distribution (left tailed)
        pValue = 1-norm.cdf(np.sqrt(np.power(t,2)))

## Adding variables to original xarray dataset 
#trend variable
    ds["trends" + ext] = (['lat','lon'], trend*factor)
    ds["trends" + ext].attrs['long_name'] = 'Linear Trend' + ext
    ds["trends" + ext].attrs['units'] = ds[ds_name].attrs['units'] + '/year'
    ds["p_values" + ext] = (['lat','lon'], pValue)
    ds["p_values" + ext].attrs['long_name'] = 'P value significance test' + ext
    ds["p_values" + ext].attrs['units'] = 'n/a'
    ds["standard_error_trend" + ext] = (['lat','lon'], np.sqrt(sb_2))
    ds["standard_error_trend" + ext].attrs['long_name'] = 'standard error of the linear trend' + ext
    ds["standard_error_trend" + ext].attrs['units'] = 'n/a'
    ds["slope" + ext] = ([time_name,'lat','lon'], slope)
    ds["slope" + ext].attrs['long_name'] = 'Linear Trend Line of best fit' + ext
    ds["slope" + ext].attrs['units'] = ds[ds_name].attrs['units']


def harmonic_analysis(ds, ds_name, freq=(0), axis_lengths=(10,10)):

    print('computing harmonic analysis for ' + ds[ds_name].attrs['long_name'])
## Extracting data and time
    data = ds[ds_name].data
    time = ds.time

## Starting harmonic analysis
    w1t = freq[0] * time
    w2t = freq[-1] * time
#trigonometric functions for harmonics
    sw1 = np.sin(w1t)
    cw1 = np.cos(w1t)
    sw2 = np.sin(w2t)
    cw2 = np.cos(w2t)
    scw1 = sw1 * cw1
    ssw1 = sw1 * sw1
    sw1cw2 = sw1 * cw2
    sw1sw2 = sw1 * sw2
    sw2cw2 = sw2 * cw2
    ccw1 = cw1 * cw1
    cw1cw2 = cw1 * cw2
    sw2cw1 = sw2 * cw1
    ccw2 = cw2 * cw2
    scw2 = sw2 * cw2
    ssw2 = sw2 * sw2
#filling B matrix
    B = np.zeros((4,4))
    B[0,0] = scw1.sum()
    B[0,1] = -ssw1.sum()
    B[0,2] = sw1cw2.sum()
    B[0,3] = -sw1sw2.sum()
    B[1,0] = ccw1.sum()
    B[1,1] = -scw1.sum()
    B[1,2] = cw1cw2.sum()
    B[1,3] = -sw2cw1.sum()
    B[2,0] = sw2cw1.sum()
    B[2,1] = -sw1sw2.sum()
    B[2,2] = scw2.sum()
    B[2,3] = -ssw2.sum()
    B[3,0] = cw1cw2.sum()
    B[3,1] = -sw1cw2.sum()
    B[3,2] = ccw2.sum()
    B[3,3] = -scw2.sum()

#making processing easier by initialising and performing calculation on whole map
    sw1Map = np.zeros(axis_lengths)
    cw1Map = np.zeros(axis_lengths)
    sw2Map = np.zeros(axis_lengths)
    cw2Map = np.zeros(axis_lengths)
    w1tMap = np.zeros(axis_lengths)
    w2tMap = np.zeros(axis_lengths)
    dsTimeMeanMap = np.zeros((axis_lengths[1], axis_lengths[2]))
    for lats in np.arange(axis_lengths[1]):
        print("\rprogress = {}%".format(int(lats/data.shape[1]*100)), end=' ', flush='true')
        for lons in np.arange(axis_lengths[2]):
            if np.isnan(data[0,lats,lons]):
                continue
            #initialising
            sw1Map[:,lats,lons] = sw1[:]
            cw1Map[:,lats,lons] = cw1[:]
            sw2Map[:,lats,lons] = sw2[:]
            cw2Map[:,lats,lons] = cw2[:]
            w2tMap[:,lats,lons] = w2t[:]
            w1tMap[:,lats,lons] = w1t[:]
            dsTimeMeanMap[lats,lons] = np.nanmean(data[:,lats,lons]) 
#matrix initialising
    A = np.zeros((4,axis_lengths[1], axis_lengths[2]))
    X = np.zeros((4,axis_lengths[1], axis_lengths[2]))
    A[0,:,:] = np.nansum(data * sw1Map, axis=0) - (dsTimeMeanMap * np.nansum(sw1, axis=0))
    A[1,:,:] = np.nansum(data * cw1Map, axis=0) - (dsTimeMeanMap * np.nansum(cw1, axis=0))
    A[2,:,:] = np.nansum(data * sw2Map, axis=0) - (dsTimeMeanMap * np.nansum(sw2, axis=0))
    A[3,:,:] = np.nansum(data * cw2Map, axis=0) - (dsTimeMeanMap * np.nansum(cw2, axis=0))
#looping to solve for unknowns at each grid point
    for lats in np.arange(axis_lengths[1]):
        print("\rprogress = {}%".format(int(lats/data.shape[1]*100)), end=' ', flush='true')
        for lons in np.arange(axis_lengths[2]):
            if np.isnan(data[0,lats,lons]):
                X[:,lats,lons] = np.nan
                continue
            #finding unknowns by inversion
            X[:,lats,lons] = np.dot(np.linalg.inv(B), np.transpose(A[:,lats,lons]))

    print("\n")
#unknown variables
    A1 = np.sqrt((X[0,:,:] * X[0,:,:]) + (X[1,:,:] * X[1,:,:]))
    Phi1 = np.arctan2(X[1,:,:],X[0,:,:])
    A2 = np.sqrt((X[2,:,:] * X[2,:,:]) + (X[3,:,:] * X[3,:,:]))
    Phi2 = np.arctan2(X[3,:,:],X[2,:,:])
#calculating harmonic time series as
#A1cos(ω₁t*Φ₁) + A2cos(ω₂t*Φ₂)
    harmonicTS = dsTimeMeanMap + A1 * np.cos(w1tMap + Phi1) + A2 * np.cos(w2tMap + Phi2)
    harmonicTS1 = dsTimeMeanMap + A1 * np.cos(w1tMap + Phi1) 
    harmonicTS2 = dsTimeMeanMap + A2 * np.cos(w2tMap + Phi2)
    R2 = (1-( (ds[ds_name]-harmonicTS).var('time').data  / ds[ds_name].var('time').data )) * 100
    R2_1 = (1-( (ds[ds_name]-harmonicTS1).var('time').data  / ds[ds_name].var('time').data )) * 100
    R2_2 = (1-( (ds[ds_name]-harmonicTS2).var('time').data  / ds[ds_name].var('time').data )) * 100
   # R2_1 = (1-(( np.nanvar(ds[ds_name].mean('time').data-harmonicTS1,axis=0) ) / np.nanvar(ds[ds_name].mean('time').data,axis=0) )) * 100
   # R2_2 = (1-(( np.nanvar(ds[ds_name].mean('time').data-harmonicTS2,axis=0) ) / np.nanvar(ds[ds_name].mean('time').data,axis=0) )) * 100
    R2[np.where(R2<0)] = 0
    R2_1[np.where(R2_1<0)] = 0
    R2_2[np.where(R2_2<0)] = 0
#calculating month of maximum amplitude
    tmax1 = ( (2 * np.pi - Phi1) / (2*np.pi) * 12)
    tmax2 = ( (2 * np.pi - Phi2) / (2*np.pi) * 6)
#correcting tmaxes to fall in correct month ranges
    tmax1[np.where(tmax1>12)] = tmax1[np.where(tmax1>12)]-12
    tmax1[np.where(tmax1<0)] = tmax1[np.where(tmax1<0)]+12
    tmax2[np.where(tmax2>6)] = tmax2[np.where(tmax2>6)]-6
    tmax2[np.where(tmax2<0)] = tmax2[np.where(tmax2<0)]+6

#harmonic analysis variables
    ds["R2"] = (['lat', 'lon'], R2)
    ds.R2.attrs['long_name'] = 'Percentage of variance by both harmonic'
    ds.R2.attrs['units'] = '%'
    ds["R2_1"] = (['lat', 'lon'], R2_1)
    ds.R2_1.attrs['long_name'] = 'Percentage of variance by annual harmonic'
    ds.R2_1.attrs['units'] = '%'
    ds["R2_2"] = (['lat', 'lon'], R2_2)
    ds.R2_2.attrs['long_name'] = 'Percentage of variance by sub-annual harmonic'
    ds.R2_2.attrs['units'] = '%'
    ds["harmonic_time_Series_all"] = (['time', 'lat', 'lon'], harmonicTS)
    ds.harmonic_time_Series_all.attrs['long_name'] = 'Both harmonic (time_series)'
    ds.harmonic_time_Series_all.attrs['units'] = ds[ds_name].attrs['units']
    ds["harmonic_time_Series_1"] = (['time', 'lat', 'lon'], harmonicTS1)
    ds.harmonic_time_Series_1.attrs['long_name'] = 'First harmonic (time_series)'
    ds.harmonic_time_Series_1.attrs['units'] = ds[ds_name].attrs['units']
    ds["harmonic_time_Series_2"] = (['time', 'lat', 'lon'], harmonicTS2)
    ds.harmonic_time_Series_2.attrs['long_name'] = 'Second harmonic (time_series)'
    ds.harmonic_time_Series_2.attrs['units'] = ds[ds_name].attrs['units']
    ds["amplitude_1"] = (['lat', 'lon'], A1)
    ds.amplitude_1.attrs['long_name'] = 'Amplitude of first harmonic'
    ds.amplitude_1.attrs['units'] = ds[ds_name].attrs['units']
    ds["amplitude_2"] = (['lat', 'lon'], A2)
    ds.amplitude_2.attrs['long_name'] = 'Amplitude of second harmonic'
    ds.amplitude_2.attrs['units'] = ds[ds_name].attrs['units']
    ds["max_amplitude_1"] = (['lat', 'lon'], tmax1)
    ds.max_amplitude_1.attrs['long_name'] = 'Month of first harmonic max amplitude'
    ds.max_amplitude_1.attrs['units'] = 'Month'
    ds["max_amplitude_2"] = (['lat', 'lon'], tmax2)
    ds.max_amplitude_2.attrs['long_name'] = 'Month of second harmonic max amplitude'
    ds.max_amplitude_2.attrs['units'] = 'Month'
