using .OpticalWavePropSim
using Plots
gr()

# example_pt_source_vac_prop.jl

delta1 = d1
deltan = d2
n = nscr

# Coordinates
vec1 = collect(-N/2 : N/2-1) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)
r1 = sqrt.(x1.^2 .+ y1.^2)

# Point Source Field
pt = exp.(-im * k / (2 * R) .* r1.^2) ./ D1^2 .* sinc.(x1 ./ D1) .* sinc.(y1 ./ D1) .* exp.(-(r1 ./ (4 * D1)).^2)

# Planes
z_planes = collect(1:n-1) .* Dz ./ (n - 1)

# Super-Gaussian Absorber
sg = exp.(-(x1 ./ (0.47 * N * d1)).^16) .* exp.(-(y1 ./ (0.47 * N * d1)).^16)

# Construct mask stack (3D Array)
t = zeros(ComplexF64, N, N, n)
for i in 1:n
    t[:, :, i] = sg
end

# Multi-plane propagation
xn, yn, Uvac = ang_spec_multi_prop(pt, wvl, delta1, deltan, z_planes, t)

# Collimation
Uvac = Uvac .* exp.(-im * π / (wvl * R) .* (xn.^2 .+ yn.^2))