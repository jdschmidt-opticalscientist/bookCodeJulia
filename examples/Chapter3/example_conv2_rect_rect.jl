# example_conv2_rect_rect.jl -------------------------------------

using Plots
using OpticalWavePropSim

N = 256         # number of samples
L = 16.0        # grid size [m]
delta = L / N   # sample spacing [m]

x_vec = ((-N/2):(N/2-1)) * delta
# Create 2D grids (Idiomatic Julia meshgrid equivalent)
x = repeat(x_vec', N, 1)
y = repeat(x_vec, 1, N)

w = 2.0         # width of rectangle
# Define 2D signals
A = rect.(x / w) .* rect.(y / w)
B = A

# Perform discrete 2D convolution
C = myconv2(A, B, delta)

# Analytic 2D convolution
C_cont = (w^2) * tri.(x / w) .* tri.(y / w)

# Plot 1: Analytic Heatmap
p1 = heatmap(
    x_vec,
    x_vec,
    C_cont,
    aspect_ratio = :equal,
    title = "Analytic",
    xlabel = "x [m]",
    ylabel = "y [m]",
    colorbar_title = "Amplitude",
    xlims = (-L / 2, L / 2),
    ylims = (-L / 2, L / 2),
)

# Plot 2: Numerical Heatmap
p2 = heatmap(
    x_vec,
    x_vec,
    real.(C),
    aspect_ratio = :equal,
    title = "Numerical",
    xlabel = "x [m]",
    ylabel = "y [m]",
    colorbar_title = "Amplitude",
    xlims = (-L / 2, L / 2),
    ylims = (-L / 2, L / 2),
)

# Plot 3: y=0 Slice Comparison
# Julia is 1-indexed, so N/2 + 1 is the center row
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
    linestyle = :dash,
)

scatter!(p3, x_vec, slice_num, marker = :x, mc = :blue, ms = 4, label = "Numerical")

# Layout: Two heatmaps on top, slice on bottom
l = @layout [grid(1, 2); c]
p_final = plot(p1, p2, p3, layout = l, size = (900, 700), margin = 5Plots.mm)

display(p_final)
