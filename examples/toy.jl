## IMPORTS
using Unitful: m, km
using GeoStats
using GeoStatsFunctions
using Plots
using VariogramVerify

## 
maxlag = 1000m
dtol = maxlag/5
angles, ranges, idx_maj, idx_min, g_list, γ_list = rangefinder(rtb, maxlag=maxlag, dtol=dtol, deg_step=10.0)
p = rangerose(angles, ranges, idx_maj, idx_min)