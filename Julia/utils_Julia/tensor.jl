# utils_Julia/tensor.jl
# using scalar sp_log, sp_arctan in per-vertex loops.

include("safe_math.jl")
using Base.Threads: @threads, nthreads

# =========================================================
# VERTEX COORDINATES & SIGNS
# =========================================================

const _VERTEX_SIGNS_TENSOR = let s = Float64[]
    for ix in (0, 1), iy in (0, 1), iz in (0, 1)
        push!(s, (-1.0)^(ix + iy + iz))
    end
    s
end

function _vertex_coordinates_tensor(x, y, z, L, B, D)
    x = Float64(x); y = Float64(y); z = Float64(z)
    L = Float64(L); B = Float64(B); D = Float64(D)

    Xb = (L - x, -(L + x))
    Yb = (B - y, -(B + y))
    Zb = (D - z, -(D + z))

    X = Vector{Float64}(undef, 8)
    Y = Vector{Float64}(undef, 8)
    Z = Vector{Float64}(undef, 8)

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

    return X, Y, Z, _VERTEX_SIGNS_TENSOR
end


# =========================================================
# CORE TENSOR FUNCTIONS (ARRAY INPUT, SCALAR SAFE MATH)
# =========================================================

function tensor_core_ln(x_rem, r_val)
    x = collect(Float64, x_rem)
    r = collect(Float64, r_val)
    out = Vector{Float64}(undef, length(x))

    @inbounds for i in eachindex(out)
        out[i] = sp_log(x[i] + r[i])
    end
    return out
end

function tensor_core_atan(x_i, x_j, x_k, r_val)
    xi = collect(Float64, x_i)
    xj = collect(Float64, x_j)
    xk = collect(Float64, x_k)
    r  = collect(Float64, r_val)

    out = Vector{Float64}(undef, length(xi))

    @inbounds for i in eachindex(out)
        num = xj[i] * xk[i]
        den = xi[i] * r[i]
        out[i] = sp_arctan(num, den)
    end
    return out
end


# =========================================================
# DIAGONAL COMPONENTS (V_xx, V_yy, V_zz)
# =========================================================

function gravitational_Vxx_cuboid(x, y, z, L, B, D, rho, G)
    X, Y, Z, S = _vertex_coordinates_tensor(x, y, z, L, B, D)
    R = similar(X)
    @inbounds @simd for i in eachindex(X)
        R[i] = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i])
    end

    atan_vals = tensor_core_atan(X, Y, Z, R)
    S_sum = sum(S .* atan_vals)
    return -Float64(G) * Float64(rho) * S_sum
end

function gravitational_Vyy_cuboid(x, y, z, L, B, D, rho, G)
    X, Y, Z, S = _vertex_coordinates_tensor(x, y, z, L, B, D)
    R = similar(X)
    @inbounds @simd for i in eachindex(X)
        R[i] = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i])
    end

    atan_vals = tensor_core_atan(Y, X, Z, R)
    S_sum = sum(S .* atan_vals)
    return -Float64(G) * Float64(rho) * S_sum
end

function gravitational_Vzz_cuboid(x, y, z, L, B, D, rho, G)
    X, Y, Z, S = _vertex_coordinates_tensor(x, y, z, L, B, D)
    R = similar(X)
    @inbounds @simd for i in eachindex(X)
        R[i] = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i])
    end

    atan_vals = tensor_core_atan(Z, X, Y, R)
    S_sum = sum(S .* atan_vals)
    return -Float64(G) * Float64(rho) * S_sum
end


# =========================================================
# OFF-DIAGONAL COMPONENTS (V_xy, V_xz, V_yz)
# =========================================================

function gravitational_Vxy_cuboid(x, y, z, L, B, D, rho, G)
    X, Y, Z, S = _vertex_coordinates_tensor(x, y, z, L, B, D)
    R = similar(X)
    @inbounds @simd for i in eachindex(X)
        R[i] = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Y[i] + 0.0) 
    end
    # Actually R must be sqrt(X^2+Y^2+Z^2). 
    @inbounds @simd for i in eachindex(X)
        R[i] = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i])
    end

    ln_vals = tensor_core_ln(Z, R)
    S_sum = sum(S .* (-ln_vals))
    return -Float64(G) * Float64(rho) * S_sum
end

function gravitational_Vxz_cuboid(x, y, z, L, B, D, rho, G)
    X, Y, Z, S = _vertex_coordinates_tensor(x, y, z, L, B, D)
    R = similar(X)
    @inbounds @simd for i in eachindex(X)
        R[i] = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i])
    end

    ln_vals = tensor_core_ln(Y, R)
    S_sum = sum(S .* (-ln_vals))
    return -Float64(G) * Float64(rho) * S_sum
end

function gravitational_Vyz_cuboid(x, y, z, L, B, D, rho, G)
    X, Y, Z, S = _vertex_coordinates_tensor(x, y, z, L, B, D)
    R = similar(X)
    @inbounds @simd for i in eachindex(X)
        R[i] = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i])
    end

    ln_vals = tensor_core_ln(X, R)
    S_sum = sum(S .* (-ln_vals))
    return -Float64(G) * Float64(rho) * S_sum
end


# =========================================================
# FULL TENSOR AT A POINT
# =========================================================

function gravitational_tensor_cuboid(x, y, z, L, B, D, rho, G)
    V_xx = gravitational_Vxx_cuboid(x, y, z, L, B, D, rho, G)
    V_yy = gravitational_Vyy_cuboid(x, y, z, L, B, D, rho, G)
    V_zz = gravitational_Vzz_cuboid(x, y, z, L, B, D, rho, G)

    V_xy = gravitational_Vxy_cuboid(x, y, z, L, B, D, rho, G)
    V_xz = gravitational_Vxz_cuboid(x, y, z, L, B, D, rho, G)
    V_yz = gravitational_Vyz_cuboid(x, y, z, L, B, D, rho, G)

    T = Matrix{Float64}(undef, 3, 3)
    T[1,1] = V_xx; T[1,2] = V_xy; T[1,3] = V_xz
    T[2,1] = V_xy; T[2,2] = V_yy; T[2,3] = V_yz
    T[3,1] = V_xz; T[3,2] = V_yz; T[3,3] = V_zz
    return T
end


# =========================================================
# BATCH + PARALLEL EVALUATION
# =========================================================

"""
    gravitational_tensor_batch(points, L, B, D, rho, G;
                               parallel_threshold=20)

Compute tensor at N points (N×3 array).

Returns Array{Float64,3} of size (N,3,3).
"""
function gravitational_tensor_batch(points, L, B, D, rho, G;
                                    parallel_threshold::Int = 20)
    pts = Array{Float64}(points)
    if ndims(pts) != 2 || size(pts, 2) != 3
        error("points must be an array of shape (N, 3)")
    end
    N = size(pts, 1)

    T = Array{Float64}(undef, N, 3, 3)

    if N > parallel_threshold && nthreads() > 1
        @threads for i in 1:N
            x = pts[i,1]; y = pts[i,2]; z = pts[i,3]
            T[i, :, :] = gravitational_tensor_cuboid(x, y, z, L, B, D, rho, G)
        end
    else
        @inbounds for i in 1:N
            x = pts[i,1]; y = pts[i,2]; z = pts[i,3]
            T[i, :, :] = gravitational_tensor_cuboid(x, y, z, L, B, D, rho, G)
        end
    end

    return T
end
