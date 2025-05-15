using Unitful: m, km
using GeoStats
using GeoStatsFunctions
using Plots

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
    perp = mod(angles[idx_major] + pi/2, pi)
    diffs  = abs.(angles .- perp)
    dists  = min.(diffs, pi .- diffs)
    idx_min = argmin(dists)

    return angles, ranges, idx_maj, idx_min, g_list, γ_list
end
