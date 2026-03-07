# example_ft_gaussian.jl -----------------------------------------

using Plots
using OpticalWavePropSim

L = 5.0
N = 32
delta = L / N
x = collect(-N/2:N/2-1) .* delta
f = collect(-N/2:N/2-1) ./ (N * delta)
a = 1.0
# sampled function & its DFT
g_samp = exp.(-π * a .* x .^ 2)
g_dft = ft(g_samp, delta)
# analytic comparison
M = 1024
x_cont = range(x[1], x[end], length=M)
f_cont = range(f[1], f[end], length=M)
g_cont = exp.(-π * a .* x_cont .^ 2)
g_ft_cont = exp.(-π .* f_cont .^ 2 ./ a) ./ a

# --- Plotting ---
p = plot(f_cont, g_ft_cont, label="Continuous FT", linecolor=:red, lw=2)
scatter!(p, f, real.(g_dft), label="DFT", markercolor=:blue, markersize=3)

xlabel!("Frequency f")
ylabel!("G(f)")
title!("Gaussian FT Comparison")
display(p)