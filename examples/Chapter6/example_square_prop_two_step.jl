using .OpticalWavePropSim
using Plots
gr()

# example_square_prop_two_step.jl
N = 1024
L = 1e-2
delta1 = L / N
D = 2e-3
wvl = 1e-6
k = 2 * π / wvl
Dz = 1.0

# Source plane setup
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)

# Square aperture
ap = rect.(x1 ./ D) .* rect.(y1 ./ D)

# Output grid spacing
delta2 = (wvl * Dz) / (N * delta1)

# Numerical propagation
x2, y2, Uout = two_step_prop(ap, wvl, delta1, delta2, Dz)

# Analytic result for y2=0 slice
mid_idx = Int(N/2) + 1
x2_slice = x2[mid_idx, :]
Uout_an = fresnel_prop_square_ap(x2_slice, 0.0, D, wvl, Dz)

# Plotting
p1 = heatmap(x2[1,:], y2[:,1], abs2.(Uout), 
             aspect_ratio=:equal, c=:inferno, 
             title="2D Irradiance")

p2 = plot(x2_slice, abs2.(Uout_an), label="Analytic", lw=2, color=:red)
scatter!(p2, x2_slice, abs2.(Uout[mid_idx, :]), 
         label="Two-Step", markersize=2, markerstrokewidth=0)
xlabel!(p2, "x [m]")
title!(p2, "1D Profile")

plot(p1, p2, layout=(1,2), size=(950, 400))