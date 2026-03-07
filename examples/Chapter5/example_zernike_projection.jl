# using .OpticalWavePropSim
# using LinearAlgebra

# example_zernike_projection.jl
N = 32
L = 2.0
delta = L / N

vec = collect(-N/2 : N/2-1) .* delta
x = ones(N) .* vec'
y = vec .* ones(1, N)
r = sqrt.(x.^2 .+ y.^2)
theta = atan.(y, x)

# Unit circle aperture
ap = circ(x, y, 2.0)

# Index map for Zernike (n, m)
z_idx = Dict(2 => (1, 1), 4 => (2, 0), 21 => (5, 1))

# 3 Zernike modes
z2 = zernike(2, r, theta, z_idx) .* ap
z4 = zernike(4, r, theta, z_idx) .* ap
z21 = zernike(21, r, theta, z_idx) .* ap

# Create aberration
W_full = 0.5 .* z2 .+ 0.25 .* z4 .- 0.6 .* z21

# Find indices inside aperture
idx = findall(x -> x > 0, ap)

# Flatten and construct matrix
W = W_full[idx]
Z = hcat(z2[idx], z4[idx], z21[idx])

# Solve the system to compute coefficients
A = Z \ W

println("Recovered Coefficients: ", A)