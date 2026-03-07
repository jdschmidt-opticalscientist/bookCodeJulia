using .OpticalWavePropSim
using Plots
gr()

# example_square_one_step_prop_samp.jl
D1 = 2e-3
D2 = 3e-3
delta1 = D1 / 50
wvl = 1e-6
Dz = 0.5

# Calculate minimum N
Nmin = (D1 * wvl * Dz) / (delta1 * (wvl * Dz - D2 * delta1))
N = Int(2^ceil(log2(abs(Nmin))))

# Source plane
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)
ap = rect.(x1 ./ D1) .* rect.(y1 ./ D1)

# Propagation
x2, y2, Uout = one_step_prop(ap, wvl, delta1, Dz)

# Analytic slice
mid = Int(N/2) + 1
x2_slice = x2[mid, :]
Uout_an = fresnel_prop_square_ap(x2_slice, 0.0, D1, wvl, Dz)

# Plotting
p1 = heatmap(x2[1,:], y2[:,1], abs2.(Uout), 
             aspect_ratio=:equal, c=:viridis, title="N=$N Propagation")

p2 = plot(x2_slice, abs2.(Uout_an), label="Analytic", color=:red, lw=2)
scatter!(p2, x2_slice, abs2.(Uout[mid, :]), label="Numerical", ms=2)
xlabel!(p2, "x [m]")

plot(p1, p2, layout=(1,2), size=(900, 400))