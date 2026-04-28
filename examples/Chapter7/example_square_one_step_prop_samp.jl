using OpticalWavePropSim
using Plots

# example_square_one_step_prop_samp.jl
D1 = 2e-3;  # diam of the source aperture [m]
D2 = 3e-3;  # diam of the obs-plane region of interest [m]
delta1 = D1 / 50;   # want at least 50 grid pts across ap
wvl = 1e-6;  # optical wavelength [m]
k = 2 * π / wvl;
Dz = 0.5;        # propagation distance [m]
# minimum number of grid points
Nmin = D1 * wvl * Dz / (delta1 * (wvl * Dz - D2 * delta1));
N = Int(2^ceil(log2(abs(Nmin))));    # number of grid pts per side
# source plane
vec1 = collect((-N / 2):(N / 2 - 1)) .* delta1
x1 = repeat(vec1', N, 1)
y1 = repeat(vec1, 1, N)
ap = rect.(x1 ./ D1) .* rect.(y1 ./ D1);
# simulate the propagation
x2, y2, Uout = one_step_prop(ap, wvl, delta1, Dz);

# analytic result for y2=0 slice
mid_idx = Int(N / 2) + 1
x2_slice = x2[mid_idx, :]
Uout_an = fresnel_prop_square_ap(x2_slice, 0.0, D1, wvl, Dz);

# --- Visualization ---
x2_slice_mm = x2_slice .* 1e3

# Irradiance Plot
p1 = plot(
    x2_slice_mm,
    abs2.(Uout_an);
    seriestype=:scatter,
    markershape=:square,
    color=:red,
    label="Analytic",
    alpha=0.6,
    markersize=3,
)
plot!(
    p1,
    x2_slice_mm,
    abs2.(Uout[mid_idx, :]);
    color=:blue,
    label="Numerical",
    title="Irradiance (N=$N)\n(y=0 slice at z=$Dz m)",
    xlabel="x₂ [mm]",
    ylabel="Irradiance [W/m²]",
    xlims=(-5, 5),
    grid=true,
)

# Phase Plot
p2 = plot(
    x2_slice_mm,
    angle.(Uout_an);
    seriestype=:scatter,
    markershape=:square,
    color=:red,
    label="Analytic",
    alpha=0.6,
    markersize=3,
)
plot!(
    p2,
    x2_slice_mm,
    angle.(Uout[mid_idx, :]);
    color=:blue,
    label="Numerical",
    title="Phase (N=$N)\n(y=0 slice at z=$Dz m)",
    xlabel="x₂ [mm]",
    ylabel="Phase [rad]",
    xlims=(-5, 5),
    grid=true,
)

# Combine and display
final_plot = plot(p1, p2; layout=(1, 2), size=(1000, 450), margin=5Plots.mm)
display(final_plot)
