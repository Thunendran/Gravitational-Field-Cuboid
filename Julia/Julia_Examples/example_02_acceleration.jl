# Julia_Examples/example_02_acceleration.jl
#
# Validates gx, gy, gz using numerical derivatives of the potential.

include(joinpath(@__DIR__, "..", "utils_Julia", "safe_math.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "surfaces.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "potential.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "acceleration.jl"))

using Printf

# ---------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------
a = 1.0
b = 1.0
c = 1.0

rho = 1.0
G   = 1.0

x = 0.0
y = 0.0
z = 0.0

h = 1e-6
results = Dict{String, Tuple{Float64, Float64, Float64}}()

# ---------------------------------------------------------
# g_x COMPARISON
# ---------------------------------------------------------
gx_analytic = gravitational_gx_cuboid(x, y, z, a, b, c, rho, G)

V_plus_x  = Gravitational_potential_cuboid(x + h, y, z, a, b, c, rho, G)
V_minus_x = Gravitational_potential_cuboid(x - h, y, z, a, b, c, rho, G)

gx_numeric = -(V_plus_x - V_minus_x) / (2h)

results["gx"] = (gx_analytic, gx_numeric, gx_analytic - gx_numeric)

# ---------------------------------------------------------
# g_y COMPARISON
# ---------------------------------------------------------
gy_analytic = gravitational_gy_cuboid(x, y, z, a, b, c, rho, G)

V_plus_y  = Gravitational_potential_cuboid(x, y + h, z, a, b, c, rho, G)
V_minus_y = Gravitational_potential_cuboid(x, y - h, z, a, b, c, rho, G)

gy_numeric = -(V_plus_y - V_minus_y) / (2h)

results["gy"] = (gy_analytic, gy_numeric, gy_analytic - gy_numeric)

# ---------------------------------------------------------
# g_z COMPARISON
# ---------------------------------------------------------
gz_analytic = gravitational_gz_cuboid(x, y, z, a, b, c, rho, G)

V_plus_z  = Gravitational_potential_cuboid(x, y, z + h, a, b, c, rho, G)
V_minus_z = Gravitational_potential_cuboid(x, y, z - h, a, b, c, rho, G)

gz_numeric = -(V_plus_z - V_minus_z) / (2h)

results["gz"] = (gz_analytic, gz_numeric, gx_analytic - gx_numeric)

# ---------------------------------------------------------
# PRINT RESULTS
# ---------------------------------------------------------

println("=====================================================")
println("  Analytical vs Numerical Gravitational Acceleration")
println("               at (0.0, 0.0, 0.0) for a Cube")
println("=====================================================")

println("\n--- g_x Component ---")
@printf("Analytical g_x   = %.15e\n", results["gx"][1])
@printf("Numerical g_x    = %.15e\n", results["gx"][2])
@printf("Difference       = %.3e\n",   results["gx"][3])

println("\n--- g_y Component ---")
@printf("Analytical g_y   = %.15e\n", results["gy"][1])
@printf("Numerical g_y    = %.15e\n", results["gy"][2])
@printf("Difference       = %.3e\n",   results["gy"][3])

println("\n--- g_z Component ---")
@printf("Analytical g_z   = %.15e\n", results["gz"][1])
@printf("Numerical g_z    = %.15e\n", results["gz"][2])
@printf("Difference       = %.3e\n",   results["gz"][3])

println("=====================================================")
