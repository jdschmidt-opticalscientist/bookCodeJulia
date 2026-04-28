using OpticalWavePropSim
using LinearAlgebra
using Printf
# example_zernike_projection.jl
N = 32    # number of grid points per side
L = 2.0      # total size of the grid [m]
delta = L / N  # grid spacing [m]
# cartesian & polar coordinates
x_vec = ((-N / 2):(N / 2 - 1)) * delta
x = repeat(x_vec', N, 1)
y = repeat(x_vec, 1, N)
r = sqrt.(x .^ 2 .+ y .^ 2)
theta = atan.(y, x)
# unit circle aperture
ap = circ.(x, y, 2.0)
# 3 Zernike modes
z2 = zernike(2, r, theta) .* ap
z4 = zernike(4, r, theta) .* ap
z21 = zernike(21, r, theta) .* ap
# create the aberration
W = 0.5 * z2 + 0.25 * z4 - 0.6 * z21
# find only grid points within the aperture
idx = findall(ap .> 0)
# perform linear indexing in column-major order
W = W[idx]
Z = hcat(z2[idx], z4[idx], z21[idx])
# solve the system of equations to compute coefficients
A = Z \ W
@printf("A = [%.4f, %.4f, %.4f]\n", A...)
