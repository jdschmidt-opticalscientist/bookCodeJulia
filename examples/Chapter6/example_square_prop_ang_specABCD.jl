using OpticalWavePropSim
using Plots

# example_square_prop_ang_specABCD.jl
N = 1024         # number of grid points per side
L = 1e-2         # total size of the grid [m]
delta1 = L / N   # grid spacing [m]
D = 2e-3         # diameter of the aperture [m]
wvl = 1e-6       # optical wavelength [m]
k = 2 * π / wvl
Dz = 1.0         # propagation distance [m]
f = Inf          # source field radius of curvature [m]

vec1 = collect((-N/2):(N/2-1)) .* delta1
x1 = repeat(vec1', N, 1)
y1 = repeat(vec1, 1, N)
ap = rect.(x1 ./ D) .* rect.(y1 ./ D)
delta2 = (wvl * Dz) / (N * delta1)

M_prop = [1.0 Dz; 0.0 1.0]
M_lens = [1.0 0.0; -1/f 1.0]
ABCD = M_prop * M_lens
x2, y2, Uout = ang_spec_propABCD(ap, wvl, delta1, delta2, ABCD)

# Analytic result for y2=0 slice
mid_idx = Int(N / 2) + 1
x2_slice = x2[mid_idx, :]
y2_val = 0.0
Uout_an = fresnel_prop_square_ap(x2_slice, y2_val, D, wvl, Dz)

# --- Visualization ---
x2_slice_mm = x2_slice .* 1e3

# Irradiance Plot
p1 = plot(
    x2_slice_mm,
    abs2.(Uout_an),
    seriestype = :scatter,
    markershape = :square,
    color = :red,
    label = "Analytic",
    alpha = 0.6,
    markersize = 3,
)
plot!(
    p1,
    x2_slice_mm,
    abs2.(Uout[mid_idx, :]),
    color = :blue,
    label = "Numerical",
    title = "Square Aperture Diffraction Irradiance\n(y=0 slice at z=1m)",
    xlabel = "x₂ [mm]",
    ylabel = "Irradiance [W/m²]",
    xlims = (-5, 5),
    grid = true,
)

# Phase Plot
p2 = plot(
    x2_slice_mm,
    angle.(Uout_an),
    seriestype = :scatter,
    markershape = :square,
    color = :red,
    label = "Analytic",
    alpha = 0.6,
    markersize = 3,
)
plot!(
    p2,
    x2_slice_mm,
    angle.(Uout[mid_idx, :]),
    color = :blue,
    label = "Numerical",
    title = "Square Aperture Diffraction Phase\n(y=0 slice at z=1m)",
    xlabel = "x₂ [mm]",
    ylabel = "Phase [rad]",
    xlims = (-5, 5),
    grid = true,
)

# Combine and display
final_plot = plot(p1, p2, layout = (1, 2), size = (1000, 450), margin = 5Plots.mm)
display(final_plot)
