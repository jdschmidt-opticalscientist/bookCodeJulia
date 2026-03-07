# using .OpticalWavePropSim
# using Plots

# example_strfcn2_rect.jl
N = 256
L = 16.0
delta = L / N

vec = collect(-N/2 : N/2-1) .* delta
x = ones(N) .* vec'
y = vec .* ones(1, N)

w = 2.0
A = rect(x ./ w) .* rect(y ./ w)
mask = ones(N, N)

# Perform discrete structure function
D_samp, idx = str_fcn2_ft(A, mask, delta)
C = D_samp ./ delta^2

# Analytic structure function
C_cont = 2 * w^2 .* (1 .- tri(x ./ w) .* tri(y ./ w))

# Visualization
# p1 = heatmap(vec, vec, C, title="Numerical D(r)", aspect_ratio=:equal)
# p2 = heatmap(vec, vec, C_cont, title="Analytic D(r)", aspect_ratio=:equal)
# plot(p1, p2, layout=(1,2))