# using .OpticalWavePropSim
# using Plots

# example_coh_img.jl
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

# Zernike index 4 (Defocus)
z_idx = zeros(4, 2)
z_idx[4, :] = [2, 0]

# Wavefront aberration
W = 0.05 .* zernike(4, 2 .* r ./ D, theta, z_idx)

# Complex pupil function
P = circ(x, y, D) .* exp.(im * 2 * π .* W)

# Amplitude spread function
h = ft2(P, delta)

# Image-plane coordinates
delta_u = wvl * z / (N * delta)
u_vec = collect(-N/2 : N/2-1) .* delta_u
u = ones(N) .* u_vec'
v = u_vec .* ones(1, N)

# Object
obj = (rect.((u .- 1.4e-4) ./ 5e-5) .+ 
       rect.(u ./ 5e-5) .+ 
       rect.((u .+ 1.4e-4) ./ 5e-5)) .* rect.(v ./ 2e-4)

# Convolve
img = myconv2(obj, h, 1.0)

# Intensity Visualization
# heatmap(u_vec, u_vec, abs2.(img), aspect_ratio=:equal, title="Coherent Image Intensity")