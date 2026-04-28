using OpticalWavePropSim
using Plots
using Printf

# example_pt_source.jl
D = 8e-3            # diameter of observation aperture [m]
wvl = 1e-6          # wavelength [m]
k = 2 * π / wvl     # optical wavenumber [rad/m]
Dz = 1.0            # propagation distance [m]
arg = D / (wvl * Dz)
delta1 = 1 / (10 * arg)   # source-plane grid spacing [m]
delta2 = D / 100          # observation-plane grid spacing [m]
N = 1024                  # number of grid points
# Source-plane coordinates
vec1 = collect((-N / 2):(N / 2 - 1)) .* delta1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)
r1 = sqrt.(x1 .^ 2 .+ y1 .^ 2)
pt = exp.(-im * k / (2 * Dz) .* r1 .^ 2) .* (arg^2) .* sinc.(arg .* x1) .* sinc.(arg .* y1)
x2, y2, Uout = ang_spec_prop(pt, wvl, delta1, delta2, Dz)

# --- Visualization ---

# Scale coordinates to mm
mid = Int(N / 2) + 1
x_slice = x2[mid, :]
x2_mm = x2 .* 1e3
y2_mm = y2 .* 1e3
x_slice_mm = x_slice .* 1e3

# Propagated Irradiance (2D)
I_out = abs2.(Uout)
p1 = heatmap(
    x2[1, :] .* 1e3,
    y2[:, 1] .* 1e3,
    I_out / 1e12;
    aspect_ratio=:equal,
    c=:grays,
    title="Numerically Propagated\nPoint-Source Irradiance",
    xlabel="x₂ [mm]",
    ylabel="y₂ [mm]",
    cb_formatter=x -> @sprintf("%.1e", x),
    colorbar_title="Irradiance / 1e12",
)

# y=0 Slice of Irradiance
I_slice = I_out[mid, :]
p2 = plot(
    x_slice_mm,
    I_slice / 1e12;
    title="Numerically Propagated\nPoint-Source Irradiance Slice",
    xlabel="x₂ [mm]",
    ylabel="Irradiance / 1e12",
    grid=true,
    legend=false,
)

# y=0 Slice of Unwrapped Phase
wrapped_phase = angle.(Uout[mid, :])
# Simple unwrap helper
function simple_unwrap!(p)
    idxs = eachindex(p)
    for i in idxs[nextind(idxs, first(idxs)):end]
        d = p[i] - p[i - 1]
        p[i] -= 2π * round(d / 2π)
    end
    return p
end
unwrapped_phase = simple_unwrap!(copy(wrapped_phase))
unwrapped_shifted = unwrapped_phase .- unwrapped_phase[mid]

expected_phase = (k / (2 * Dz) .* x_slice .^ 2)

p3 = plot(
    x_slice_mm,
    expected_phase;
    line=(:solid, 2, 0.8),
    color=:red,
    label="Analytic",
    title="Numerically Propagated\nPoint-Source Phase",
    xlabel="x₂ [mm]",
    ylabel="Phase [rad]",
    grid=true,
)

plot!(p3, x_slice_mm, unwrapped_shifted; line=(:dashdot), color=:blue, label="Numerical")

# Combine plots into a 1x3 layout
final_plot = plot(p1, p2, p3; layout=(1, 3), size=(1200, 400), margin=5Plots.mm)

display(final_plot)
