# using .OpticalWavePropSim
# using LinearAlgebra

# example_zernike_synthesis.jl
N = 40
L = 2.0
delta = L / N

vec = collect(-N/2 : N/2-1) .* delta
x = ones(N) .* vec'
y = vec .* ones(1, N)
r = sqrt.(x.^2 .+ y.^2)
theta = atan.(y, x)

ap = circ(x, y, 2.0)
idx_ap = findall(x -> x > 0, ap)

# Create atmospheric phase screen
r0 = L / 20
# Assuming ft_phase_screen returns a real array
phz = ft_phase_screen(r0, N, delta, Inf, 0.0)
screen_full = (phz ./ (2 * π)) .* ap
W = screen_full[idx_ap]

# Analyze Screen
n_modes = 100
Z = zeros(length(W), n_modes)

for idx in 1:n_modes
    # z_idx_map should be pre-defined
    temp = zernike(idx, r, theta, z_idx_map)
    Z[:, idx] = temp[idx_ap]
end

# Compute mode coefficients (Backslash operator for least squares)
A = Z \ W

# Reconstruct
W_prime = Z * A

# Map back to 2D grid
scr = zeros(N, N)
scr[idx_ap] = W_prime