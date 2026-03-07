# using .OpticalWavePropSim
# using Plots

# example_conv2_rect_rect.jl
N = 256
L = 16.0
delta = L / N

vec = collect(-N/2 : N/2-1) .* delta
x = ones(N) .* vec'
y = vec .* ones(1, N)

w = 2.0
# Define 2D signals
A = rect(x ./ w) .* rect(y ./ w)
B = A

# Perform discrete 2D convolution
C = myconv2(A, B, delta)

# Analytic 2D convolution
C_cont = (w^2) .* tri(x ./ w) .* tri(y ./ w)

# Optional: Visualization
# heatmap(vec, vec, real(C), title="2D Convolution Result", aspect_ratio=:equal)