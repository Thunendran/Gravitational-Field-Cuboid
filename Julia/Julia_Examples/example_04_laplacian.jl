###############################################################################
# Example 04 — Laplacian Test (∇²V = Vxx + Vyy + Vzz)
# Replicates Python Matplotlib version exactly using PyPlot.jl
###############################################################################

include(joinpath(@__DIR__, "..", "utils_Julia", "safe_math.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "surfaces.jl"))
include(joinpath(@__DIR__, "..", "utils_Julia", "tensor.jl"))

using PyPlot

# --------------------------------------------------------------------
# Matplotlib STYLE SETTINGS
# --------------------------------------------------------------------
PyPlot.rc("text", usetex=true)
PyPlot.rc("font", family="serif", serif=["Computer Modern Roman"])
PyPlot.rc("axes", labelsize=16)
PyPlot.rc("font", size=16)
PyPlot.rc("legend", fontsize=16)
PyPlot.rc("xtick", labelsize=15)
PyPlot.rc("ytick", labelsize=15)

# --------------------------------------------------------------------
# PARAMETERS
# --------------------------------------------------------------------
G   = 1.0
rho = 1.0

L, B, D = 2.0, 2.0, 2.0        # full dimensions
X0, Y0, Z0 = 0.0, 0.0, 0.0     # cube center

a, b, c = L/2, B/2, D/2        # semi-dimensions for analytic functions

z_plane = Z0                   # z = 0 plane

# --------------------------------------------------------------------
# ANALYTIC LAPLACIAN FUNCTION
# --------------------------------------------------------------------
function laplacian_analytic(x, y, z)
    x_rel = x - X0
    y_rel = y - Y0
    z_rel = z - Z0

    Vxx = gravitational_Vxx_cuboid(x_rel, y_rel, z_rel, a, b, c, rho, G)
    Vyy = gravitational_Vyy_cuboid(x_rel, y_rel, z_rel, a, b, c, rho, G)
    Vzz = gravitational_Vzz_cuboid(x_rel, y_rel, z_rel, a, b, c, rho, G)

    return Vxx + Vyy + Vzz
end

# --------------------------------------------------------------------
# GRID (300 × 300)
# --------------------------------------------------------------------
x_vals = LinRange(X0 - L*0.75, X0 + L*0.75, 300)
y_vals = LinRange(Y0 - B*0.75, Y0 + B*0.75, 300)

X = [x for x in x_vals, y in y_vals]
Y = [y for x in x_vals, y in y_vals]
laplace_vals = zeros(size(X))

# --------------------------------------------------------------------
# COMPUTE LAPLACIAN AT GRID POINTS
# --------------------------------------------------------------------
for i in axes(X,1)
    for j in axes(X,2)
        laplace_vals[i,j] = laplacian_analytic(X[i,j], Y[i,j], z_plane)
    end
end

# --------------------------------------------------------------------
# DISCRETE COLORMAP 
# --------------------------------------------------------------------
N = 50
colors = Array{Tuple{Float64,Float64,Float64}}(undef, N)

colors[1] = (0.1, 0.4, 1.0)   # blue start

for i in 2:N-1
    ratio = (i-2) / (N-3)
    r = ratio
    g = 1.0
    b = 1.0 - ratio
    colors[i] = (r, g, b)
end

colors[end] = (1.0, 0.0, 0.0) # red end

cmap_discrete = PyPlot.matplotlib.colors.ListedColormap(colors)

vmin = -4 * pi
vmax = 0.0
bounds = range(vmin, vmax, length=N+1)
norm = PyPlot.matplotlib.colors.BoundaryNorm(bounds, N)

# --------------------------------------------------------------------
# PLOTTING 
# --------------------------------------------------------------------
fig, ax = subplots(figsize=(10,8))

im = ax.imshow(
    laplace_vals,
    extent=(x_vals[1], x_vals[end], y_vals[1], y_vals[end]),
    origin="lower",
    cmap=cmap_discrete,
    norm=norm,
    interpolation="nearest"
)

# Colorbar

cbar = plt.colorbar(im)

# Resize colorbar height to match plot height
pos = ax.get_position()
cbar.ax.set_position([
    pos.x1 + 0.02,   # x-position (slightly right of main plot)
    pos.y0,          # bottom aligned
    0.02,            # width
    pos.height       # EXACT same height as main plot
])

#cbar = colorbar(im, ax=ax, boundaries=bounds)
tick_positions = [-12, -10, -8, -6, -4, -2, 0]
cbar.set_ticks(tick_positions)
cbar.set_ticklabels(string.(tick_positions))
cbar.set_label("\$\\nabla^2 V\$", fontsize=14)

# ---------------------------------------------------------
# DRAW CUBE OUTLINE (white)
# ---------------------------------------------------------
xmin, xmax = X0 - L/2, X0 + L/2
ymin, ymax = Y0 - B/2, Y0 + B/2

ax.plot(
    [xmin, xmax, xmax, xmin, xmin],
    [ymin, ymin, ymax, ymax, ymin],
    color="white",
    linewidth=2,
    label="Cube Boundary"
)

legend = ax.legend(loc="lower right", frameon=false)
for txt in legend.get_texts()
    txt.set_color("white")
    txt.set_fontsize(16)
end

# Text labels
ax.text(0, 0, "\$\\nabla^2 V = -4\\pi\$", color="white",
        fontsize=20, ha="center", va="center", fontweight="bold")

outside_y = Y0 + B * 0.6
ax.text(0, outside_y, "\$\\nabla^2 V = 0\$", color="white",
        fontsize=20, ha="center", va="center", fontweight="bold")

ax.set_title("Laplacian Test in the x-y Plane Through the Cube", fontsize=17)
ax.set_xlabel("x", fontsize=16)
ax.set_ylabel("y", fontsize=16)
ax.grid(false)

# Save
savefig("laplacian_plot_Julia.png", dpi=300, bbox_inches="tight")
show()
