# example_ft_gaussian.jl -----------------------------------------
using Plots
using OpticalWavePropSim

# function values to be used in DFT
L = 5.0             # spatial extent of the grid
N = 32              # number of samples
delta = L / N       # sample spacing
x = ((-N/2):(N/2-1)) * delta
f = ((-N/2):(N/2-1)) / (N * delta)
a = 1.0

# sampled function & its DFT
g_samp = exp.(-π * a * x .^ 2)   # function samples
g_dft = ft(g_samp, delta)      # DFT

# analytic function & its continuous FT
M = 1024
x_cont = range(x[1], stop = x[end], length = M)
f_cont = range(f[1], stop = f[end], length = M)
g_cont = exp.(-π * a * x_cont .^ 2)
g_ft_cont = exp.(-π * f_cont .^ 2 / a) / a

# --- Plotting (Matching Python Style) ---

# Initialize the plot with the continuous analytic FT
p = plot(
    f_cont,
    g_ft_cont,
    linecolor = :red,
    linewidth = 1.5,
    label = "Continuous FT (Analytic)",
    title = "Fourier Transform of a Gaussian: Numerical vs. Analytic",
    xlabel = "Frequency f [cycles/m]",
    ylabel = "G(f)",
    grid = true,
    linestyle = :dash,
    gridalpha = 0.7,
    legend = :topright,
    size = (800, 500),
)

# Add the DFT result as discrete points (matching Python's 'bx')
scatter!(
    p,
    f,
    real.(g_dft),
    marker = :x,
    markercolor = :blue,
    markersize = 4,
    label = "DFT (Numerical)",
)

# Explicitly display the plot object
display(p)
