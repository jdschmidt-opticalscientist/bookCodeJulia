# example_ft_gaussian_shift.jl -----------------------------------

using Plots
using OpticalWavePropSim

L = 10.0      # spatial extent of the grid
N = 64        # number of samples
delta = L / N # sample spacing
x = ((-N / 2):(N / 2 - 1)) * delta
x0 = 5 * delta
f = ((-N / 2):(N / 2 - 1)) / (N * delta)
a = 1.0

# sampled function & its DFT
g_samp = exp.(-π * a * (x .- x0) .^ 2) # function samples
g_dft = ft(g_samp, delta)            # DFT

# analytic function and its continuous FT
M = 1024
x_cont = range(x[1]; stop=x[end], length=M)
f_cont = range(f[1]; stop=f[end], length=M)
g_cont = exp.(-π * a * (x_cont .- x0) .^ 2)
g_ft_cont = exp.(-im * 2 * π * x0 * f_cont) .* exp.(-π * f_cont .^ 2 / a) / a

# Initialize plot with the first series and Python-style formatting
p = plot(
    f_cont,
    real.(g_ft_cont);
    linecolor=:red,
    linewidth=1.5,
    label="Analytic FT (Real)",
    title="Fourier Transform of a Gaussian: Numerical vs. Analytic",
    xlabel="Frequency f [cycles/m]",
    ylabel="G(f)",
    grid=true,
    linestyle=:dash,      # Matches Python's '--'
    gridalpha=0.7,        # Matches Python's alpha=0.7
    size=(800, 500),      # Matches Python's figsize=(10, 6) aspect
    legend=:topright,
)

# Add remaining analytic lines
plot!(p, f_cont, imag.(g_ft_cont); linecolor=:green, lw=1.5, label="Analytic FT (Imag)")
plot!(p, f_cont, abs.(g_ft_cont); linecolor=:blue, lw=1.5, label="Analytic FT (Abs)")

# Add DFT result as discrete points (matching Python's 'rx', 'gx', 'bx')
scatter!(p, f, real.(g_dft); marker=:x, markercolor=:red, markersize=4, label="DFT (Real)")
scatter!(
    p, f, imag.(g_dft); marker=:x, markercolor=:green, markersize=4, label="DFT (Imag)"
)
scatter!(p, f, abs.(g_dft); marker=:x, markercolor=:blue, markersize=4, label="DFT (Abs)")

# Explicitly display the plot
display(p)
