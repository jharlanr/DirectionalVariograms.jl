## ============= ##
## LOAD PACKAGES ##
## ============= ##
import GLMakie as Mke
using Revise
using GeoStats
using Unitful: m, km
using Random
using DataFrames
using GeoIO
using CSV
using NPZ
using Plots
using MultiGrid

## ========= ##
## LOAD DATA ##
## ========= ##
csv_path = "/Users/jrines/stanford_gp/research/mx/MultiGrid/data/radio_lines.csv"
vtb, grid = csv2gtb(csv_path, EPSG_code=32735, gridres=50m)
vtb_ps, grid_ps = gtb_grid2pointset(vtb, grid)
vtb_ps |> viewer

## =============== ##
## GENERATE TRENDS ##
## =============== ##
ttbs_grid, ttbs_obspts = make_trends(input_gtb=vtb_ps, 
                                     n_trends=1, 
                                     geom4trend=grid_ps, 
                                     obspts=vtb_ps.geometry, 
                                     gridsize=size(grid), 
                                     strata_size=300, 
                                     maxlag=500m, 
                                     nlags=100, 
                                     nugget=0.0)

# SELECT TREND FROM TREND OPTIONS
ttbs_grid_copy = copy(ttbs_grid)
ttbs_obspts_copy = copy(ttbs_obspts)
ttb_obspts = ttbs_obspts_copy[1]
ttb_grid = ttbs_grid_copy[1]

ttb_grid |> viewer

## ================= ##
## COMPUTE RESIDUALS ##
## ================= ##
rtb = compute_residuals(ttb_obspts, vtb_ps)


## ============================================= ##
## COMPUTE VARIOGRAM FOR ENTIRE SET OF RESIDUALS ##
## ============================================= ##



