# Julia_Examples/example_01_potential.jl

include(joinpath(@__DIR__, "..", "utils_Julia", "safe_math.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "surfaces.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "potential.jl"))

using Printf

# Constants and point of evaluation
G   = 1.0
rho = 1.0

X_0, Y_0, Z_0 = 0.0, 0.0, 0.0
L, B, D       = 1.0, 1.0, 1.0     # half-dimensions, same as Python
x, y, z       = 0.0, 0.0, 0.0

V = Gravitational_potential_cuboid(
        x - X_0,
        y - Y_0,
        z - Z_0,
        L/2, B/2, D/2,
        rho, G
    )

@printf("The gravitational potential at point (%.3f, %.3f, %.3f) is %.13e m^2/s^2\n",
        x, y, z, V)
