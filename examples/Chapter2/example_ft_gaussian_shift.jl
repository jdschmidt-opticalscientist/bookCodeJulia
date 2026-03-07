# example_ft_gaussian_shift.jl -----------------------------------

using Plots
using OpticalWavePropSim

L = 10.0  # spatial extent of the grid
N = 64    # number of samples
delta = L / N  # sample spacing
x = collect(-N/2:N/2-1) .* delta
f = collect(-N/2:N/2-1) ./ (N * delta)
a = 1.0
x0 = 5 * delta
# sampled function & its DFT
g_samp = exp.(-π * a .* (x .- x0) .^ 2) # function samples
g_dft = ft(g_samp, delta) # DFT
# analytic function and its continuous FT
M = 1024
x_cont = range(x[1], x[end], length=M)
f_cont = range(f[1], f[end], length=M)
g_ft_cont = exp.(-im * 2 * π * x0 .* f_cont) .* exp.(-π .* f_cont .^ 2 ./ a) ./ a

# --- Plotting ---
p = plot(f_cont, real.(g_ft_cont), label="Analytic FT (Real)", color=:red, lw=2)
plot!(p, f_cont, imag.(g_ft_cont), label="Analytic FT (Imag)", color=:green, lw=2)
plot!(p, f_cont, abs.(g_ft_cont), label="Analytic FT (Abs)", color=:blue, lw=2)
scatter!(p, f, real.(g_dft), label="DFT (Real)", marker=:x, markercolor=:red, markersize=3)
scatter!(p, f, imag.(g_dft), label="DFT (Imag)", marker=:x, markercolor=:green, markersize=3)
scatter!(p, f, abs.(g_dft), label="DFT (Abs)", marker=:x, markercolor=:blue, markersize=3)
xlabel!("Frequency f [cyc/m]")
ylabel!("G(f)")
title!("Fourier Transform of a Gaussian: Numerical vs. Analytic", titlefont=font(12))
display(p)