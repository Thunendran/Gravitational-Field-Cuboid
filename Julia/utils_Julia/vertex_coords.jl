# utils_Julia/vertex_coords.jl

const VERTEX_SIGNS = let s = Float64[]
    for ix in 0:1, iy in 0:1, iz in 0:1
        push!(s, (-1.0)^(ix + iy + iz))
    end
    s
end

"""
    vertex_coordinates(x, y, z, L, B, D)

Return (X, Y, Z, S) where X,Y,Z are length-8 Float64 vectors
for the oriented vertex distances, and S are the ±1 signs.
"""
function vertex_coordinates(x::Float64, y::Float64, z::Float64,
                            L::Float64, B::Float64, D::Float64)
    Xb = (L - x, -(L + x))
    Yb = (B - y, -(B + y))
    Zb = (D - z, -(D + z))

    X = Vector{Float64}(undef, 8)
    Y = Vector{Float64}(undef, 8)
    Z = Vector{Float64}(undef, 8)

    idx = 1
    @inbounds for ix in 1:2, iy in 1:2, iz in 1:2
        X[idx] = Xb[ix]
        Y[idx] = Yb[iy]
        Z[idx] = Zb[iz]
        idx += 1
    end

    return X, Y, Z, VERTEX_SIGNS
end
