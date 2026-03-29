using Plots
using OpticalWavePropSim

N = 512              # number of grid points per side
L = 7.5e-3           # total size of the grid [m]
d1 = L / N           # source-plane grid spacing [m]
D = 1e-3             # diameter of the aperture [m]
wvl = 1e-6           # optical wavelength [m]
k = 2 * π / wvl
Dz = 20.0            # propagation distance [m]

vec1 = (-N/2:N/2-1) * d1
x1 = repeat(vec1', N, 1)
y1 = repeat(vec1, 1, N)
Uin = circ.(x1, y1, D)
Uout, x2, y2 = fraunhofer_prop(Uin, wvl, d1, Dz)

# analytic Result
Uout_th = (exp.(im * k / (2 * Dz) .* (x2 .^ 2 .+ y2 .^ 2))
           /
           (im * wvl * Dz) * (D^2 * π / 4)
           .*
           jinc.(D .* sqrt.(x2 .^ 2 .+ y2 .^ 2) ./ (wvl * Dz)))


# 1D slice extraction at y = 0
mid = Int(N / 2) + 1
x2_vec = x2[mid, :]
slice_num = abs2.(Uout[mid, :]) * 1e3     # numerical irradiance [mW/m^2]
slice_ana = abs2.(Uout_th[mid, :]) * 1e3  # analytic irradiance [mW/m^2]

# Plot 1: Numerical Irradiance Heatmap
# x2[1,:] and y2[:,1] provide the axis vectors for the heatmap
p1 = heatmap(x2[1, :], y2[:, 1], abs2.(Uout) * 1e3,
    aspect_ratio=:equal,
    c=:viridis,
    title="Numerical Irradiance",
    xlabel="x₂ [m]", ylabel="y₂ [m]",
    colorbar_title="mW/m²")

# Plot 2: y=0 Slice Comparison
p2 = plot(x2_vec, slice_ana,
    color=:red, lw=2,
    label="Analytic",
    title="y₂=0 Slice",
    xlabel="x₂ [m]", ylabel="Irradiance [mW/m²]",
    grid=true)

scatter!(p2, x2_vec[1:8:end], slice_num[1:8:end], # Downsample scatter for clarity
    marker=:x, mc=:blue, ms=4,
    label="Numerical")

# Final Layout
p_final = plot(p1, p2, layout=(1, 2), size=(1100, 500), margin=5Plots.mm)

display(p_final)