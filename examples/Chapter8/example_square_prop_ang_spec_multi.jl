using .OpticalWavePropSim
using Plots
gr()

# example_square_prop_ang_spec_multi.jl
D1 = 2e-3
D2 = 6e-3
wvl = 1e-6
z_total = 2.0
delta1 = D1 / 30
deltan = D2 / 30
N = 128
n_steps = 5

# Array of distances
z_steps = collect(1:n_steps) .* z_total ./ n_steps

# Source plane
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)
ap = rect.(x1 ./ D1) .* rect.(y1 ./ D1)

# Multi-plane propagation
x2, y2, Uout = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z_steps)

# Analytic result
mid = Int(N/2) + 1
x2_slice = x2[mid, :]
Uout_an = fresnel_prop_square_ap(x2_slice, 0.0, D1, wvl, z_total)

# Plotting
p1 = heatmap(x2[1,:], y2[:,1], abs2.(Uout), 
             aspect_ratio=:equal, c=:viridis, title="Multi-plane ASM")

p2 = plot(x2_slice, abs2.(Uout_an), label="Analytic", lw=2, color=:red)
scatter!(p2, x2_slice, abs2.(Uout[mid, :]), label="Numerical ($n_steps steps)", ms=3)
xlabel!(p2, "x [m]")

plot(p1, p2, layout=(1,2), size=(900, 400))