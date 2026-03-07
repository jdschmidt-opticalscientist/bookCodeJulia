# using .OpticalWavePropSim
# using Plots

# example_derivative_ft.jl
N = 64
L = 6.0
delta = L / N
x = collect(-N/2 : N/2-1) .* delta

w = 3.0
window = rect(x ./ w)
g = (x.^5) .* window

# Discrete derivatives
gp_samp = real.(derivative_ft(g, delta, 1)) .* window
gpp_samp = real.(derivative_ft(g, delta, 2)) .* window

# Analytic derivatives
gp = 5 .* (x.^4) .* window
gpp = 20 .* (x.^3) .* window

# Visualization
# p1 = plot(x, gp_samp, seriestype=:scatter, label="FT")
# plot!(p1, x, gp, label="Analytic", title="1st Deriv")
# p2 = plot(x, gpp_samp, seriestype=:scatter, label="FT")
# plot!(p2, x, gpp, label="Analytic", title="2nd Deriv")
# plot(p1, p2, layout=(1,2))