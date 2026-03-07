using Plots
gr()

# analysis_pt_source_atmos_samp.jl

# Effective diameters
c = 2.0
D1p = D1 + c * wvl * Dz / r0sw
D2p = D2 + c * wvl * Dz / r0sw

d1_vec = range(1e-6, 1.1 * wvl * Dz / D2p, length=100)
dn_vec = range(1e-6, 1.1 * wvl * Dz / D1p, length=100)

# Constraint boundaries
dn_max1 = -D2p / D1p .* d1_vec .+ wvl * Dz / D1p
dn_min3 = (1 + Dz / R) .* d1_vec .- wvl * Dz / D1p
dn_max3 = (1 + Dz / R) .* d1_vec .+ wvl * Dz / D1p

# Constraint 4: Numerical check
d1, d2 = 10e-3, 10e-3
N = 512
zmax = min(d1, d2)^2 * N / wvl
nmin = ceil(Dz / zmax) + 1

println("Minimum screens needed for N=512: $nmin")

# Plotting the sampling window
p = plot(d1_vec .* 1000, dn_max1 .* 1000, label="Boundary 1", color=:red, lw=2)
plot!(p, d1_vec .* 1000, dn_min3 .* 1000, label="Nyquist Min", color=:green, style=:dash)
plot!(p, d1_vec .* 1000, dn_max3 .* 1000, label="Nyquist Max", color=:green)
xlabel!("delta 1 [mm]")
ylabel!("delta n [mm]")
title!("Atmospheric Propagation Sampling Constraints")

display(p)