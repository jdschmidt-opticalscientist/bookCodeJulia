using OpticalWavePropSim
using Plots

# example_square_prop_ang_spec_multi.jl
D1 = 2e-3;   # diameter of the source aperture [m]
D2 = 6e-3;   # diameter of the observation aperture [m]
wvl = 1e-6;  # optical wavelength [m]
k = 2 * π / wvl; # optical wavenumber [rad/m]
z = 2;     # propagation distance [m]
delta1 = D1 / 30; # source-plane grid spacing [m]
deltan = D2 / 30; # observation-plane grid spacing [m]
N = 128;        # number of grid points
n = 5;          # number of partial propagations
# switch from total distance to individual distances
# collect() is used to create a concrete array for the propagation function
z_array = collect(1:n) .* z ./ n;
# source-plane coordinates
vec1 = collect((-N / 2):(N / 2 - 1)) .* delta1
x1 = repeat(vec1', N, 1)
y1 = repeat(vec1, 1, N)
ap = rect.(x1 ./ D1) .* rect.(y1 ./ D1);    # source aperture
x2, y2, Uout = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z_array);

# analytic result for y2=0 slice
mid_idx = Int(N / 2) + 1
x2_slice = x2[mid_idx, :]
Dz = z_array[end]; # switch back to total distance
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
    title="Irradiance (Multi-step ASM)\n(y=0 slice at z=$Dz m)",
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
    title="Phase (Multi-step ASM)\n(y=0 slice at z=$Dz m)",
    xlabel="x₂ [mm]",
    ylabel="Phase [rad]",
    xlims=(-5, 5),
    grid=true,
)

# Combine and display
final_plot = plot(p1, p2; layout=(1, 2), size=(1000, 450), margin=5Plots.mm)
display(final_plot)
