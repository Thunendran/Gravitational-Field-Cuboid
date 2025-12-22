###############################################################################
# Example 03 — Gravitational Acceleration Field (z = 0 Plane)
# Matplotlib-style plot using PyPlot.jl 
###############################################################################

include(joinpath(@__DIR__, "..", "utils_Julia", "safe_math.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "surfaces.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "potential.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "acceleration.jl"))

using PyPlot

# ---------------------------------------------------------
# Matplotlib 
# ---------------------------------------------------------
PyPlot.rc("text", usetex=true)
PyPlot.rc("font", family="serif", serif=["Computer Modern Roman"])
PyPlot.rc("axes", labelsize=16)
PyPlot.rc("font", size=16)
PyPlot.rc("legend", fontsize=16)
PyPlot.rc("xtick", labelsize=15)
PyPlot.rc("ytick", labelsize=15)

# ---------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------
G   = 1.0
rho = 1.0

X0, Y0, Z0 = 0.0, 0.0, 0.0
L, B, D   = 2.0, 2.0, 2.0    # full lengths

a = L/2    # semi-dimensions (expected by analytic functions)
b = B/2
c = D/2

# ---------------------------------------------------------
# GRID IN PLANE Z = 0
# ---------------------------------------------------------
x_vals = LinRange(-2, 2, 13)
y_vals = LinRange(-2, 2, 13)

X = [x for x in x_vals, y in y_vals]
Y = [y for x in x_vals, y in y_vals]

U = zeros(size(X))
V = zeros(size(Y))

# ---------------------------------------------------------
# ANALYTICAL ACCELERATION FIELD
# ---------------------------------------------------------
for i in axes(X,1)
    for j in axes(X,2)
        xq = X[i,j]
        yq = Y[i,j]
        zq = 0.0

        gx = gravitational_gx_cuboid(xq, yq, zq, a, b, c, rho, G)
        gy = gravitational_gy_cuboid(xq, yq, zq, a, b, c, rho, G)

        U[i,j] = gx
        V[i,j] = gy
    end
end

# ---------------------------------------------------------
# PLOT
# ---------------------------------------------------------
fig, ax = subplots(figsize=(8,7))

ax.quiver(X, Y, U, V, color="red", scale=65, pivot="middle", width=0.003)

# Cuboid boundary
x_min, x_max = X0 - L/2, X0 + L/2
y_min, y_max = Y0 - B/2, Y0 + B/2

ax.plot(
    [x_min, x_max, x_max, x_min, x_min],
    [y_min, y_min, y_max, y_max, y_min],
    color="black", linewidth=1, label="Cube Boundary"
)

legend = ax.legend(loc="lower right", frameon=false)
for txt in legend.get_texts()
    txt.set_color("black")
    txt.set_fontsize(16)
end

ax.set_title("Gravitational Acceleration in the Plane z = 0",
             fontsize=17, color="black")
ax.set_xlabel("x", fontsize=16, color="black")
ax.set_ylabel("y", fontsize=16, color="black")
ax.tick_params(axis="both", colors="black")

ax.set_xlim([-2.1, 2.1])
ax.set_ylim([-2.1, 2.1])
ax.set_xticks(-2:1:2)
ax.set_yticks(-2:1:2)
ax.set_aspect("equal")
ax.grid(false)

savefig("Gravitational_Acceleration_quiver_z0_julia.png", dpi=300,
        bbox_inches="tight")
show()
