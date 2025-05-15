function rangerose(angles::AbstractVector{<:Real}, ranges::AbstractVector{<:Real}, idx_maxrange::Int=1)
    """
    rangerose(angles, ranges, idx_maxrange=1)

    Plots a rose diagram of variogram range vs. direction.
    - Highlights the maximum range direction in red ("Major axis").
    - Highlights the orthogonal direction in black ("Minor axis").
    - Shows the max-range angle (in degrees) in the title.
    """

    θ_max = angles[idx_maxrange]
    r_max = ranges[idx_maxrange]
    θ_orth = mod(θ_max + π/2, π)

    idx_orth = findfirst(θ -> isapprox(θ, θ_orth; atol=1e-3), angles)
    r_orth = isnothing(idx_orth) ? 0.0 : ranges[idx_orth]

    deg_max = round(rad2deg(θ_max); digits=1)

    # Base rose plot
    p = plot(angles, ranges;
             proj        = :polar,
             seriestype  = :stem,
             marker      = :circle,
             linewidth   = 2,
             legend      = :topright,
             label       = "All directions",
             title       = "Anisotropy Range-Rose (max @ $(deg_max)°)")

    # Major axis (red)
    plot!([θ_max], [r_max];
          proj        = :polar,
          seriestype  = :stem,
          linewidth   = 3,
          color       = :red,
          marker      = :circle,
          markercolor = :red,
          label       = "Major axis")

    # Minor axis (black)
    plot!([θ_orth], [r_orth];
          proj        = :polar,
          seriestype  = :stem,
          linewidth   = 2,
          color       = :black,
          marker      = :circle,
          markercolor = :black,
          label       = "Minor axis (perp to Major axis)")

    return p
end