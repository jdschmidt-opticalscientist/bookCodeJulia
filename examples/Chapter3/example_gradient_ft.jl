using Plots
using OpticalWavePropSim

# --- Setup ---
N = 64
L = 6.0
delta = L / N
x_vec = ((-N/2):(N/2-1)) * delta
x = repeat(x_vec', N, 1)
y = repeat(x_vec, 1, N)

# 2D Gaussian Function: g = exp(-(x^2 + y^2))
g = exp.(-(x .^ 2 + y .^ 2))

# --- Gradient Calculation ---
# 1. Numerical (FT-based)
gx_num, gy_num = gradient_ft(g, delta)
gx_num, gy_num = real.(gx_num), real.(gy_num)

# 2. Analytic: ∇g = [-2x * g, -2y * g]
gx_ana = -2 .* x .* g
gy_ana = -2 .* y .* g

# Function to prepare downsampled, flattened vectors for quiver
function get_quiver_data(qx, qy, step, x_mat, y_mat)
    idx = 1:step:size(qx, 1)
    return vec(x_mat[idx, idx]), vec(y_mat[idx, idx]), vec(qx[idx, idx]), vec(qy[idx, idx])
end

step = 4
px, py, pu_num, pv_num = get_quiver_data(gx_num, gy_num, step, x, y)
_, _, pu_ana, pv_ana = get_quiver_data(gx_ana, gy_ana, step, x, y)

# Plot 1: Numerical Result
p1 = heatmap(x_vec, x_vec, g, c = :viridis, title = "Numerical (FT)", aspect_ratio = :equal)
quiver!(p1, px, py, quiver = (pu_num, pv_num), color = :white, lw = 1.0)

# Plot 2: Analytic Result
p2 = heatmap(
    x_vec,
    x_vec,
    g,
    c = :viridis,
    title = "Analytic (-2x·g, -2y·g)",
    aspect_ratio = :equal,
)
quiver!(p2, px, py, quiver = (pu_ana, pv_ana), color = :white, lw = 1.0)

# Combine into 1x2 layout
p_final = plot(p1, p2, layout = (1, 2), size = (1100, 500), margin = 5Plots.mm)

display(p_final)
