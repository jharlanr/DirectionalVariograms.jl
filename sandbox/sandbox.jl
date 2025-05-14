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

θ = 90*π/180  # angle of anisotropy

# 1) Build two orthogonal direction‐vectors at your target angle θ:
dir1 =  ( cos(θ),  sin(θ) )   # “major” direction
dir2 = ( -sin(θ),  cos(θ) )   # ⟂ to dir1

# 2) Compute directional variograms along those:
g1 = DirectionalVariogram(dir1, rtb, "val",
                          estimator=:cressie,
                          maxlag=1000m, dtol=300m)
γ1 = GeoStatsFunctions.fit(MaternVariogram, g1)

g2 = DirectionalVariogram(dir2, rtb, "val",
                          estimator=:cressie,
                          maxlag=1000m, dtol=300m)
γ2 = GeoStatsFunctions.fit(MaternVariogram, g2)

# 3) Build your anisotropic variogram with those true radii:
ellipsoid = MetricBall((γ1.ball.radii[1], γ2.ball.radii[1]), Angle2d(θ))
γ_ani    = MaternVariogram(ellipsoid,
                           sill   = γ1.sill,
                           nugget = γ1.nugget)

g_ani = EmpiricalVariogram(rtb, "val", estimator=:cressie, maxlag=1500m, nlags=10)
# 4) Plot the variogram:
funplot!(funplot(g_ani), γ_ani, maxlag=1500m, color="red")
