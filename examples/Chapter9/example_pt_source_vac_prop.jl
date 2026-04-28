using OpticalWavePropSim, Plots, Statistics

# example_pt_source_vac_prop.m
delta1 = d1;    # source-plane grid spacing [m]
deltan = d2;    # observation-plane grid spacing [m]
n = nscr;       # number of planes

# coordinates
v_coord = ((-N/2):(N/2-1)) * delta1;
x1 = repeat(v_coord', N, 1);
y1 = repeat(v_coord, 1, N);
r1 = sqrt.(x1 .^ 2 .+ y1 .^ 2);

# point source
# Equation: $$ U(r) = \frac{1}{D_1^2} \exp\left(\frac{-ik r^2}{2R}\right) \text{sinc}\left(\frac{x}{D_1}\right) \text{sinc}\left(\frac{y}{D_1}\right) \exp\left(-\left(\frac{r}{4D_1}\right)^2\right) $$
pt =
    exp.(-im * k / (2 * R) .* r1 .^ 2) ./ D1^2 .* sinc.(x1 / D1) .* sinc.(y1 / D1) .*
    exp.(-(r1 / (4 * D1)) .^ 2);

# partial prop planes
z = collect(1:(n-1)) * Dz / (n - 1);

# simulate vacuum propagation
sg = exp.(-(x1 / (0.47 * N * d1)) .^ 16) .* exp.(-(y1 / (0.47 * N * d1)) .^ 16);

# Stack the super-gaussian mask for each plane
t = repeat(sg, 1, 1, n);

xn, yn, Uvac = ang_spec_multi_prop(pt, wvl, delta1, deltan, z, t);

# collimate the beam
Uvac = Uvac .* exp.(-im * pi / (wvl * R) .* (xn .^ 2 .+ yn .^ 2));
# --- Visualization ---

# Calculate irradiance
Ivac = abs2.(Uvac);

p1 = heatmap(
    xn[1, :] .* 1e3,
    yn[:, 1] .* 1e3,
    Ivac,
    aspect_ratio = :equal,
    color = :inferno,
    title = "Vacuum Irradiance",
    xlabel = "x [mm]",
    ylabel = "y [mm]",
)

p2 = heatmap(
    xn[1, :] .* 1e3,
    yn[:, 1] .* 1e3,
    angle.(Uvac),
    aspect_ratio = :equal,
    color = :phase,
    title = "Vacuum Phase",
    xlabel = "x [mm]",
    ylabel = "y [mm]",
)

l = @layout [a b]
p_final = plot(p1, p2, layout = l, size = (900, 400))
display(p_final)
