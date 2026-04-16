using OpticalWavePropSim
using LinearAlgebra
using Plots

# example_zernike_synthesis.jl
N = 40;     # number of grid points per side
L = 2.0;      # total size of the grid [m]
delta = L / N;  # grid spacing [m]
# cartesian & polar coordinates
x_vec = (-N/2:N/2-1) * delta
x = repeat(x_vec', N, 1)
y = repeat(x_vec, 1, N)
r = sqrt.(x .^ 2 .+ y .^ 2)
theta = atan.(y, x)
# unit circle aperture
ap = circ.(x, y, 2.0)
# indices of grid points in aperture
idxAp = findall(ap .> 0)
# create atmospheric phase screen
r0 = L / 20;
screen = ft_phase_screen(r0, N, delta, Inf, 0.0) ./ (2 * π) .* ap;
W = screen[idxAp];   # perform linear indexing

# analyze screen
nModes = 100;   # number of Zernike modes
# create matrix of Zernike polynomial values
Z = zeros(length(W), nModes);
for idx = 1:nModes
    temp = zernike(idx, r, theta)
    Z[:, idx] = temp[idxAp]
end
# compute mode coefficients
A = Z \ W;
# synthesize mode-limited screen
W_prime = Z * A;
# reshape mode-limited screen into 2-D for display
scr = zeros(N, N);
scr[idxAp] = W_prime;

# --- Plots ---

# Common color limits for comparison
c_min = minimum(W)
c_max = maximum(W)

# 1. Original Screen
p_orig = heatmap(x_vec, x_vec, screen,
    aspect_ratio=:equal, c=:jet, clims=(c_min, c_max),
    title="Original Screen", xlabel="m", ylabel="m")

# 2. Cumulative Reconstructions
mode_list = [3, 16, 36, 100]
p_recons = []

for n in mode_list
    # Sum modes 1 through n
    W_n = Z[:, 1:n] * A[1:n]

    scr_n = zeros(N, N)
    scr_n[idxAp] = W_n

    p = heatmap(x_vec, x_vec, scr_n,
        aspect_ratio=:equal,
        c=:jet,
        clims=(c_min, c_max),
        title="Modes 1 to $n",
        xlabel="x [m]",
        ylabel="y [m]",
        colorbar=:none)
    push!(p_recons, p)
end

# Layout: Original on top, 2x2 grid of reconstructions below
l = @layout [
    a{0.3h}
    grid(2, 2)
]

final_plot = plot(p_orig, p_recons..., layout=l, size=(800, 1000), margin=5Plots.mm)
display(final_plot)