# using .OpticalWavePropSim
# using Plots

# example_conv_rect_rect.jl
N = 64
L = 8.0
delta = L / N
x = collect(-N/2 : N/2-1) .* delta

w = 2.0
A = rect(x ./ w)
B = A

# Perform discrete convolution
C = myconv(A, B, delta)

# Continuous convolution
C_cont = w .* tri(x ./ w)

# Visualization
# plot(x, real(C), seriestype=:scatter, label="Discrete", title="Rect-Rect Convolution")
# plot!(x, C_cont, label="Analytic")