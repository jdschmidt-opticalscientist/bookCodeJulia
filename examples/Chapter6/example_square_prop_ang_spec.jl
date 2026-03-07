# using .OpticalWavePropSim

# example_square_prop_ang_spec.jl
N = 1024
L = 1e-2
delta1 = L / N
D = 2e-3
wavelength = 1e-6
k = 2 * π / wavelength
Dz = 1.0

# Source plane
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)

# Square aperture
ap = rect.(x1 ./ D) .* rect.(y1 ./ D)

# Observation scaling
delta2 = (wavelength * Dz) / (N * delta1)

# Numerical propagation
x2, y2, Uout = ang_spec_prop(ap, wavelength, delta1, delta2, Dz)

# Analytic result for y2=0 slice
# In Julia, N/2 + 1 is the 1-based center index
mid_idx = Int(N / 2) + 1
x2_slice = x2[mid_idx, :]
Uout_an = fresnel_prop_square_ap(x2_slice, 0.0, D, wavelength, Dz)