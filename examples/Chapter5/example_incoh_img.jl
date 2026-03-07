# using .OpticalWavePropSim
# using Plots

# example_incoh_img.jl
N = 256
L = 0.1
D = 0.07
delta = L / N
wvl = 1e-6
z = 0.25

# Pupil-plane coordinates
vec = collect(-N/2 : N/2-1) .* delta
x = ones(N) .* vec'
y = vec .* ones(1, N)
r = sqrt.(x.^2 .+ y.^2)
theta = atan.(y, x)

# Zernike index 4 map
z_idx = zeros(4, 2); z_idx[4, :] = [2, 0]
W = 0.05 .* zernike(4, 2 .* r ./ D, theta, z_idx)

# Pupil and ASF
P = circ.(x, y, D) .* exp.(im * 2 * π .* W)
h = ft2(P, delta)
psf = abs2.(h)

# Image-plane coordinates
U_grid = wvl * z / (N * delta)
u_vec = collect(-N/2 : N/2-1) .* U_grid
u = ones(N) .* u_vec'
v = u_vec .* ones(1, N)

# Object intensity
obj_int = (rect.((u .- 1.4e-4) ./ 5e-5) .+ 
           rect.(u ./ 5e-5) .+ 
           rect.((u .+ 1.4e-4) ./ 5e-5)) .* rect.(v ./ 2e-4)

# Convolve intensities
img_incoh = myconv2(obj_int, psf, 1.0)

# heatmap(u_vec, u_vec, real(img_incoh), cmap=:grays, title="Incoherent Image")