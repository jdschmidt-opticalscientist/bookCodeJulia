# using .OpticalWavePropSim
# using Plots

# example_corr2_rect_rect.jl
N = 256
L = 16.0
delta = L / N

vec = collect(-N/2 : N/2-1) .* delta
x = ones(N) .* vec'
y = vec .* ones(1, N)

w = 2.0
# Define 2D signal
A = rect(x ./ w) .* rect(y ./ w)
mask = ones(N, N)

# Perform discrete 2D correlation
C, idx = corr2_ft(A, A, mask, delta)

# Analytic correlation
C_cont = (w^2) .* tri(x ./ w) .* tri(y ./ w)

# Visualization
# heatmap(vec, vec, real(C), title="2D Correlation Result", aspect_ratio=:equal)