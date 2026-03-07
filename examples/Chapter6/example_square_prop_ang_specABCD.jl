# using .OpticalWavePropSim

# example_square_prop_ang_specABCD.jl
N = 1024
L = 1e-2
delta1 = L / N
D = 2e-3
wvl = 1e-6
k = 2 * π / wvl
Dz = 1.0
f = Inf

# Source plane
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)

# Square aperture
ap = rect.(x1 ./ D) .* rect.(y1 ./ D)

# Observation scaling
delta2 = (wvl * Dz) / (N * delta1)

# Matrix multiplication in Julia
M_prop = [1.0 Dz; 0.0 1.0]
M_lens = [1.0 0.0; -1/f 1.0]
ABCD = M_prop * M_lens

# Numerical propagation
x2, y2, Uout = ang_spec_propABCD(ap, wvl, delta1, delta2, ABCD)