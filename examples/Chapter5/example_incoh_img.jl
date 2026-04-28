using Plots
using OpticalWavePropSim

N = 256             # number of grid points per side
L = 0.1             # total size of the grid [m]
D = 0.07            # diameter of pupil [m]
delta = L / N       # grid spacing [m]
wvl = 1e-6          # optical wavelength [m]
z = 0.25            # image distance [m]
# pupil-plane coordinates
vec = ((-N / 2):(N / 2 - 1)) * delta
x = repeat(vec', N, 1)
y = repeat(vec, 1, N)
r = sqrt.(x .^ 2 .+ y .^ 2)
theta = atan.(y, x)
# wavefront aberration
z_idx = Dict(4 => [2, 0])
W = 0.05 .* zernike(4, 2 .* r ./ D, theta, z_idx)
# complex pupil function
P = circ.(x, y, D) .* exp.(im * 2 * π .* W)
# amplitude spread function & point spread function
h = ft2(P, delta)
psf = abs2.(h)      # Incoherent Point Spread Function (Intensity)
# image-plane coordinates
U_grid = wvl * z / (N * delta)
u_vec = ((-N / 2):(N / 2 - 1)) * U_grid
u = repeat(u_vec', N, 1)
v = repeat(u_vec, 1, N)
# object (same coordinates as h)
obj =
    (rect.((u .- 1.4e-4) ./ 5e-5) .+ rect.(u ./ 5e-5) .+ rect.((u .+ 1.4e-4) ./ 5e-5)) .*
    rect.(v ./ 2e-4)
# convolve the object with the PSF to simulate imaging
img = myconv2(obj, psf, 1.0)

# Convert to mm for plot axes
u_mm = u_vec * 1e3

# Plot 1: Object Intensity
p1 = heatmap(
    u_mm,
    u_mm,
    obj;
    aspect_ratio=:equal,
    c=:viridis,
    title="Object Irradiance",
    xlabel="x [mm]",
    ylabel="y [mm]",
)

# Plot 2: Point Spread Function (PSF)
p2 = heatmap(
    u_mm,
    u_mm,
    psf;
    aspect_ratio=:equal,
    c=:viridis,
    title="PSF",
    xlabel="u [mm]",
    ylabel="v [mm]",
)

# Plot 3: Incoherent Image Irradiance
p3 = heatmap(
    u_mm,
    u_mm,
    real.(img);
    aspect_ratio=:equal,
    c=:viridis,
    title="Image Irradiance",
    xlabel="u [mm]",
    ylabel="v [mm]",
)

# Final Layout
p_final = plot(p1, p2, p3; layout=(1, 3), size=(1300, 450), margin=5Plots.mm)

display(p_final)
