# utils_Julia/acceleration.jl
#
# using scalar sp_log, sp_arctan and vertex-based sums.

include("safe_math.jl")
using Base.Threads: @threads, nthreads

# ================================================================
# Precomputed vertex signs  S = (-1)^{i_x+i_y+i_z}
# ================================================================
const _VERTEX_SIGNS_ACCEL = let s = Float64[]
    for ix in (0, 1), iy in (0, 1), iz in (0, 1)
        push!(s, (-1.0)^(ix + iy + iz))
    end
    s
end

# ================================================================
# Oriented vertex coordinates (signed face distances)
# ================================================================
function _vertex_coordinates_accel(x, y, z, L, B, D)
    x = Float64(x); y = Float64(y); z = Float64(z)
    L = Float64(L); B = Float64(B); D = Float64(D)

    Xb = (L - x, -(L + x))
    Yb = (B - y, -(B + y))
    Zb = (D - z, -(D + z))

    X = Vector{Float64}(undef, 8)
    Y = Vector{Float64}(undef, 8)
    Z = Vector{Float64}(undef, 8)

    # same pattern as:
    # X = np.repeat(Xb, 4)
    # Y = np.tile(np.repeat(Yb,2),2)
    # Z = np.tile(Zb,4)
    idx = 1
    @inbounds for ix in 1:2
        for iy in 1:2
            for iz in 1:2
                X[idx] = Xb[ix]
                Y[idx] = Yb[iy]
                Z[idx] = Zb[iz]
                idx += 1
            end
        end
    end

    return X, Y, Z, _VERTEX_SIGNS_ACCEL
end


# ================================================================
# Core acceleration kernel  g_core (per-vertex)
# ================================================================
function g_core(x_j, x_k, x_l, X_comp, Y_comp, Z_comp)
    xj = collect(Float64, x_j)
    xk = collect(Float64, x_k)
    xl = collect(Float64, x_l)
    Xc = collect(Float64, X_comp)
    Yc = collect(Float64, Y_comp)
    Zc = collect(Float64, Z_comp)

    out = Vector{Float64}(undef, length(xj))

    @inbounds for i in eachindex(out)
        Xci = Xc[i]
        Yci = Yc[i]
        Zci = Zc[i]

        r = sqrt(Xci*Xci + Yci*Yci + Zci*Zci)

        term1 = xj[i] * sp_log(xk[i] + r)
        term2 = xk[i] * sp_log(xj[i] + r)

        num   = xj[i] * xk[i]
        den   = xl[i] * r
        term3 = -xl[i] * sp_arctan(num, den)

        out[i] = term1 + term2 + term3
    end

    return out
end


# ================================================================
# Gravitational acceleration components
# ================================================================

function gravitational_gx_cuboid(x, y, z, L, B, D, rho_vol, G)
    X, Y, Z, S = _vertex_coordinates_accel(x, y, z, L, B, D)

    # Kernel arguments:
    #   x_l = X  (axis)
    #   x_j = Y
    #   x_k = Z
    kernel = g_core(Y, Z, X, X, Y, Z)

    return -Float64(G) * Float64(rho_vol) * sum(S .* kernel)
end


function gravitational_gy_cuboid(x, y, z, L, B, D, rho_vol, G)
    X, Y, Z, S = _vertex_coordinates_accel(x, y, z, L, B, D)

    # Python: g_core(x_j=X, x_k=Z, x_l=Y, X_comp=X, Y_comp=Y, Z_comp=Z)
    kernel = g_core(X, Z, Y, X, Y, Z)

    return -Float64(G) * Float64(rho_vol) * sum(S .* kernel)
end


function gravitational_gz_cuboid(x, y, z, L, B, D, rho_vol, G)
    X, Y, Z, S = _vertex_coordinates_accel(x, y, z, L, B, D)

    # Python: g_core(x_j=X, x_k=Y, x_l=Z, X_comp=X, Y_comp=Y, Z_comp=Z)
    kernel = g_core(X, Y, Z, X, Y, Z)

    return -Float64(G) * Float64(rho_vol) * sum(S .* kernel)
end


# ================================================================
# Convenience wrapper: acceleration vector at a point
# ================================================================
function gravitational_acceleration_point(x, y, z, L, B, D, rho_vol, G)
    gx = gravitational_gx_cuboid(x, y, z, L, B, D, rho_vol, G)
    gy = gravitational_gy_cuboid(x, y, z, L, B, D, rho_vol, G)
    gz = gravitational_gz_cuboid(x, y, z, L, B, D, rho_vol, G)
    return gx, gy, gz
end


# ================================================================
# Parallel & batch evaluation
# ================================================================

"""
    gravitational_acceleration_batch(points, L, B, D, rho_vol, G;
                                     parallel_threshold=20)

points: N×3 array [[x₁,y₁,z₁], ..., [xN,yN,zN]].

Returns (gx, gy, gz) vectors of length N.

Uses multithreading if N > parallel_threshold.
"""
function gravitational_acceleration_batch(points, L, B, D, rho_vol, G;
                                          parallel_threshold::Int = 20)
    pts = Array{Float64}(points)
    if ndims(pts) != 2 || size(pts, 2) != 3
        error("points must be an array of shape (N, 3)")
    end
    N = size(pts, 1)

    gx = Vector{Float64}(undef, N)
    gy = Vector{Float64}(undef, N)
    gz = Vector{Float64}(undef, N)

    if N > parallel_threshold && nthreads() > 1
        @threads for i in 1:N
            x = pts[i,1]; y = pts[i,2]; z = pts[i,3]
            gx[i], gy[i], gz[i] = gravitational_acceleration_point(x, y, z, L, B, D, rho_vol, G)
        end
    else
        @inbounds for i in 1:N
            x = pts[i,1]; y = pts[i,2]; z = pts[i,3]
            gx[i], gy[i], gz[i] = gravitational_acceleration_point(x, y, z, L, B, D, rho_vol, G)
        end
    end

    return gx, gy, gz
end
