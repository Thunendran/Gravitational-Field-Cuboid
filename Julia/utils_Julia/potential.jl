# utils_Julia/potential.jl

include("surfaces.jl")
using Base.Threads: @threads, nthreads

# =========================================================
# SINGLE-POINT GRAVITATIONAL POTENTIAL
# =========================================================

function Gravitational_potential_cuboid(x, y, z, L, B, D, rho, G)
    x = Float64(x); y = Float64(y); z = Float64(z)
    L = Float64(L); B = Float64(B); D = Float64(D)
    rho = Float64(rho); G = Float64(G)

    X1 = (L - x)
    X2 = (L + x)
    X3 = (B - y)
    X4 = (B + y)
    X5 = (D - z)
    X6 = (D + z)

    V =
        Front(X1, X2, X3, X4, X5, X6) +
        Back(X1, X2, X3, X4, X5, X6) +
        Right(X1, X2, X3, X4, X5, X6) +
        Left(X1, X2, X3, X4, X5, X6) +
        Top(X1, X2, X3, X4, X5, X6) +
        Bottom(X1, X2, X3, X4, X5, X6)

    return G * rho * V
end


# =========================================================
# PARALLEL + VECTOR BATCH POTENTIAL EVALUATION
# =========================================================

function potential_batch(points, L, B, D, rho, G; parallel_threshold::Int = 20)
    pts = Array{Float64}(points)
    if ndims(pts) != 2 || size(pts, 2) != 3
        error("points must be shape (N, 3)")
    end
    N = size(pts, 1)

    V = Vector{Float64}(undef, N)

    if N > parallel_threshold && nthreads() > 1
        @threads for i in 1:N
            x = pts[i, 1]; y = pts[i, 2]; z = pts[i, 3]
            V[i] = Gravitational_potential_cuboid(x, y, z, L, B, D, rho, G)
        end
    else
        @inbounds for i in 1:N
            x = pts[i, 1]; y = pts[i, 2]; z = pts[i, 3]
            V[i] = Gravitational_potential_cuboid(x, y, z, L, B, D, rho, G)
        end
    end

    return V
end
