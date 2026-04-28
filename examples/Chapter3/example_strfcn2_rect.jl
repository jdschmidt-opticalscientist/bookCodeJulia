using Plots
using OpticalWavePropSim

# --- Parameters ---
N = 256         # number of samples
L = 16.0        # grid size [m]
delta = L / N   # sample spacing [m]
F = 1 / L         # frequency-domain grid spacing [1/m]

x_vec = ((-N/2):(N/2-1)) * delta
# Create 2D grids
x = repeat(x_vec', N, 1)
y = repeat(x_vec, 1, N)

w = 2.0         # width of rectangle
A = rect.(x / w) .* rect.(y / w)  # signal
mask = ones(N, N)

# --- Structure Function Calculation ---
# Discrete 2D structure function calculation
# str_fcn2_ft returns (D_real, idx)
D_samp, idx = str_fcn2_ft(A, mask, delta)

# Normalize by delta^2 as per the original MATLAB logic
C = D_samp / delta^2

# Analytic structure function: 2 * w^2 * (1 - tri(x/w) * tri(y/w))
C_cont = 2 * (w^2) * (1 .- tri.(x / w) .* tri.(y / w))

# --- Visualization ---

# Top Left: Analytic Heatmap
p1 = heatmap(
    x_vec,
    x_vec,
    C_cont,
    aspect_ratio = :equal,
    title = "Analytic",
    xlabel = "x [m]",
    ylabel = "y [m]",
    xlims = (-L / 2, L / 2),
    ylims = (-L / 2, L / 2),
)

# Top Right: Numerical Heatmap
p2 = heatmap(
    x_vec,
    x_vec,
    real.(C),
    aspect_ratio = :equal,
    title = "Numerical",
    xlabel = "x [m]",
    ylabel = "y [m]",
    xlims = (-L / 2, L / 2),
    ylims = (-L / 2, L / 2),
)

# Bottom: y=0 Slice Comparison
mid = Int(N / 2) + 1
slice_num = real.(C[mid, :])
slice_ana = C_cont[mid, :]

p3 = plot(
    x_vec,
    slice_ana,
    color = :red,
    marker = :square,
    ms = 3,
    lw = 1.5,
    label = "Analytic",
    title = "Cross-section comparison at y=0",
    xlabel = "x [m]",
    ylabel = "Amplitude",
    grid = true,
)

scatter!(p3, x_vec, slice_num, marker = :x, mc = :blue, ms = 4, label = "Numerical")

# Layout: Two heatmaps on top, one wide plot on bottom
l = @layout [grid(1, 2); c]
p_final = plot(p1, p2, p3, layout = l, size = (900, 700), margin = 5Plots.mm)

display(p_final)
