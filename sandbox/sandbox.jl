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






# ## ROSE PLOT
# using GeoStats, GeoStatsFunctions, Plots
# ##
# angles = 0:pi/18:(pi - pi/18)

# ranges = Float64[]

# for θ in angles
#     dir_maj = (cos(θ), sin(θ))
#     g   = DirectionalVariogram(dir_maj, rtb, "val";
#                                estimator=:cressie,
#                                maxlag=1000m, dtol=300m)
#     γ   = GeoStatsFunctions.fit(MaternVariogram, g)

#     # strip units by dividing by 1m → a unitless Float64
#     push!(ranges, γ.ball.radii[1] / (1m))
# end

# for θ_deg in (0, 90, 180, 270)
#     θ     = θ_deg * π/180
#     dir1  = (cos(θ), sin(θ))
#     println("θ = $θ_deg → dir_maj = $dir1")
# end


# ##

# # 2) Now plot using the rotated angles directly:
# using Plots

# p = plot(angles, ranges;
#          proj       = :polar,
#          seriestype = :stem,
#          marker     = :circle,
#          linewidth  = 2,
#          legend     = false,
#          title      = "Anisotropy Range-Rose",
#          # no need for theta_offset/direction now
# )

# display(p)



# # angles sweep in radians from 0 (East) counter‐clockwise
# angles = 0:pi/18:(pi - pi/18)

# ranges = Float64[]
# for θ in angles
#     # 0°→(1,0) on X; 90°→(0,1) on Y
#     dir_maj = (cos(θ), sin(θ))

#     g   = DirectionalVariogram(dir_maj, rtb, "val";
#                                estimator = :cressie,
#                                maxlag     = 1000m,
#                                dtol       = 300m)
#     γ   = GeoStatsFunctions.fit(MaternVariogram, g)

#     push!(ranges, γ.ball.radii[1] / (1m))
# end




##

function rangefinder(gtb::GeoTable; maxlag::Quantity=1km, dtol::Quantity=300m, deg_step::Real=10.0)
    """
    rangefinder(gtb::GeoTable, maxlag::Quantity=1km, dtol::Quantity=300m, deg_step::Real=10.0)
    computes the range of the variogram in all directions (0-180 degrees) and returns a vector of ranges
    for each angle in the range of 0-180 degrees. The angles are in radians and the ranges are unitless.
    """
    angles = 0:pi/(180/deg_step):(pi - pi/(180/deg_step))
    ranges = Float64[]
    g_list  = Any[]
    γ_list  = Any[]

    for θ in angles
        dir_maj = (cos(θ), sin(θ))
        g   = DirectionalVariogram(dir_maj, gtb, "val";
                                estimator=:cressie,
                                maxlag=maxlag, dtol=dtol)
        γ   = GeoStatsFunctions.fit(MaternVariogram, g)

        # strip units by dividing by 1m → a unitless Float64
        push!(ranges, γ.ball.radii[1] / (1m))
        push!(g_list,  g)
        push!(γ_list,  γ)
    end

    # get index of maximum range
    idx_maj = argmax(ranges)

    # get index of minor axis (perp to major axis)
    perp = mod(angles[idx_maj] + pi/2, pi)
    diffs  = abs.(angles .- perp)
    dists  = min.(diffs, pi .- diffs)
    idx_min = argmin(dists)

    return angles, ranges, idx_maj, idx_min, g_list, γ_list
end

function rangerose(angles::AbstractVector{<:Real}, ranges::AbstractVector{<:Real}, idx_maj::Int=1, idx_min::Int=1)
    """
    rangerose(angles, ranges, idx_maxrange, idx_minrange)

    Plots a rose diagram of variogram range vs. direction.
    - Highlights the maximum range direction in red ("Major axis").
    - Highlights the specified minor axis direction in black ("Minor axis").
    - Shows the max- and min-range angles (in degrees) in the title.

    # Arguments
    - `angles`: Vector of angles in radians (0…π).
    - `ranges`: Corresponding variogram ranges (unitless Float64).
    - `idx_maj`: Index of the major axis (max range).
    - `idx_min`: Index of the minor axis (perpendicular direction).

    # Returns
    - A Plots.jl `Plot` object of the anisotropy range-rose.
    """

    # Extract major and minor axes
    θ_max = angles[idx_maj]
    r_max = ranges[idx_maj]
    θ_min = angles[idx_min]
    r_min = ranges[idx_min]

    # Convert to degrees for labels
    deg_max = round(rad2deg(θ_max); digits=1)
    deg_min = round(rad2deg(θ_min); digits=1)

    # Base rose plot of all directions
    p = plot(angles, ranges;
             proj        = :polar,
             seriestype  = :stem,
             marker      = :circle,
             linewidth   = 2,
             legend      = :topright,
             label       = "All directions",
             title       = "Anisotropy Range-Rose (max @ $(deg_max)°, min @ $(deg_min)°)")

    # Major axis in red
    plot!([θ_max], [r_max];
          proj        = :polar,
          seriestype  = :stem,
          linewidth   = 3,
          color       = :red,
          marker      = :circle,
          markercolor = :red,
          label       = "Major axis @ $(deg_max)°")

    # Minor axis in black
    plot!([θ_min], [r_min];
          proj        = :polar,
          seriestype  = :stem,
          linewidth   = 2,
          color       = :black,
          marker      = :circle,
          markercolor = :black,
          label       = "Minor axis @ $(deg_min)°")

    return p
end


## CALL THE FUNCTIONS
maxlag = 1000m
dtol = maxlag/5
angles, ranges, idx_maj, idx_min, g_list, γ_list = rangefinder(rtb, maxlag=maxlag, dtol=dtol, deg_step=10.0)
p = rangerose(angles, ranges, idx_maj, idx_min)





## PLOT variograms



## 
# using VariogramVerify

