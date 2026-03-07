# using .OpticalWavePropSim

# example_square_prop_one_step.jl
N = 1024
L = 1e-2
delta1 = L / N
D = 2e-3
wvl = 1e-6
k = 2 * π / wvl
Dz = 1.0

# Source plane
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)

# Square aperture
ap = rect.(x1 ./ D) .* rect.(y1 ./ D)

# Numerical propagation
x2, y2, Uout = one_step_prop(ap, wvl, delta1, Dz)

# Analytic result for y2=0 slice
mid_idx = Int(N/2) + 1
x2_slice = x2[mid_idx, :]
Uout_an = fresnel_prop_square_ap(x2_slice, 0.0, D, wvl, Dz)

# Plotting (requires Plots.jl)
# plot(x2_slice, abs2.(Uout[mid_idx, :]), seriestype=:scatter, label="One-Step")
# plot!(x2_slice, abs2.(Uout_an), label="Analytic")