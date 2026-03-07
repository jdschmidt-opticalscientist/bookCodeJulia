# using .OpticalWavePropSim

# example_pt_source.jl
D = 8e-3
wvl = 1e-6
k = 2 * π / wvl
Dz = 1.0

arg = D / (wvl * Dz)
delta1 = 1 / (10 * arg)
delta2 = D / 100
N = 1024

# Source-plane coordinates
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)
r1 = sqrt.(x1.^2 .+ y1.^2)

# Band-limited point source
# Note: Julia's sinc(x) is also sin(πx)/(πx)
pt = exp.(-im * k / (2 * Dz) .* r1.^2) .* (arg^2) .* sinc.(arg .* x1) .* sinc.(arg .* y1)

# Propagate
x2, y2, Uout = ang_spec_prop(pt, wvl, delta1, delta2, Dz)