using .OpticalWavePropSim
using Plots
gr()

# example_square_prop_ang_spec.jl
D1 = 2e-3
D2 = 4e-3
wvl = 1e-6
Dz = 0.1

delta1 = 9.4848e-6
delta2 = 28.1212e-6

# Minimum N calculation
Nmin = D1/(2*delta1) + D2/(2*delta2) + (wvl*Dz)/(2*delta1*delta2)
N = Int(2^ceil(log2(Nmin)))

# Source plane
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)
ap = rect.(x1 ./ D1) .* rect.(y1 ./ D1)

# Propagation
x2, y2, Uout = ang_spec_prop(ap, wvl, delta1, delta2, Dz)

# Analytic slice
mid = Int(N/2) + 1
x2_slice = x2[mid, :]
Uout_an = fresnel_prop_square_ap(x2_slice, 0.0, D1, wvl, Dz)

# Plotting
p1 = heatmap(x2[1,:], y2[:,1], abs2.(Uout), 
             aspect_ratio=:equal, c=:inferno, title="ASM N=$N")

p2 = plot(x2_slice, abs2.(Uout_an), label="Analytic", lw=2, color=:red)
scatter!(p2, x2_slice, abs2.(Uout[mid, :]), label="ASM", ms=2, color=:blue)
xlabel!(p2, "x [m]")

plot(p1, p2, layout=(1,2), size=(900, 400))