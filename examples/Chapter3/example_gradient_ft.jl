# using .OpticalWavePropSim
# using Plots

# example_gradient_ft.jl
N = 64
L = 6.0
delta = L / N
vec = collect(-N/2 : N/2-1) .* delta
x = ones(N) .* vec'
y = vec .* ones(1, N)

# 2D Gaussian
g = exp.(-(x.^2 .+ y.^2))

# Discrete derivatives
gx_samp, gy_samp = gradient_ft(g, delta)
gx_samp = real.(gx_samp)
gy_samp = real.(gy_samp)

# Analytic derivatives
gx = -2 .* x .* exp.(-(x.^2 .+ y.^2))
gy = -2 .* y .* exp.(-(x.^2 .+ y.^2))

# Visualization (using Plots.jl)
# p1 = heatmap(vec, vec, gx_samp, title="Numerical d/dx", aspect_ratio=:equal)
# p2 = heatmap(vec, vec, gx, title="Analytic d/dx", aspect_ratio=:equal)
# plot(p1, p2, layout=(1,2))