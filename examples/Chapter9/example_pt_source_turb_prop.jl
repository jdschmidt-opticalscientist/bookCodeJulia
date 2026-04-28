using OpticalWavePropSim, Plots, Printf, Statistics

# example_pt_source_turb_prop.m
l0 = 0.0;     # inner scale [m]
L0 = Inf;     # outer scale [m]
zt = vcat(0.0, z);  # propagation plane locations
Delta_z = diff(zt);    # propagation distances

# grid spacings
alpha_path = zt / zt[end];
delta_path = (1 .- alpha_path) * delta1 + alpha_path * deltan;

# initialize array for phase screens
phz = zeros(Float64, N, N, n);
nreals = 50;    # number of random realizations (increased for better stats)

# initialize arrays for propagated fields,
# aperture mask, and MCF
Uout = zeros(ComplexF64, N, N);
mask = circ.(xn / D2, yn / D2, 1);
MCF2 = zeros(ComplexF64, N, N);

# sg from vacuum script is reused here
for idxreal in 1:nreals     # loop over realizations
    if mod(idxreal, 10) == 0 || idxreal == 1
        @printf("Realization %d of %d\n", idxreal, nreals)
    end

    # loop over screens
    for idxscr in 1:n
        phz_lo, phz_hi = ft_sh_phase_screen(r0scrn[idxscr], N, delta_path[idxscr], L0, l0)
        phz[:, :, idxscr] = phz_lo .+ phz_hi
    end

    # Disambiguate scope for the variables returned by the propagator
    global xn, yn, Uout

    # simulate turbulent propagation
    xn, yn, Uout = ang_spec_multi_prop(pt, wvl, delta1, deltan, z, t .* exp.(im .* phz))

    # collimate the beam
    Uout = Uout .* exp.(-im * pi / (wvl * R) * (xn .^ 2 .+ yn .^ 2))

    # accumulate realizations of the MCF
    # Note: Using the exported function from your package
    global MCF2
    MCF2_current, _ = corr2_ft(Uout, Uout, mask, deltan)
    MCF2 .+= MCF2_current
end

# modulus of the complex degree of coherence
mid = Int(N / 2) + 1;
MCDOC2 = abs.(MCF2) ./ abs(MCF2[mid, mid]);

# --- Final Comparison with Theory ---

# Theoretical coherence for Spherical Wave (SW)
# Equation: $$ \text{MCDOC}_{th} = \exp\left[-3.44 \left(\frac{\rho}{r_0}\right)^{5/3}\right] $$
rho = abs.(xn[mid, mid:(end - 1)]);
MCDOC_th = exp.(-3.44 * (rho ./ r0sw) .^ (5 / 3));

# Plotting the results
p_mcf = plot(
    rho ./ r0sw,
    MCDOC2[mid, mid:(end - 1)];
    label="Simulation (Julia)",
    color=:orange,
    lw=2,
    title="Coherence (MCDOC) with Atmosphere",
    xlabel="Separation / r0",
    ylabel="Spatial Coherence Factor",
    xlims=(0, 4),
    ylims=(0, 1.1),
    grid=true,
)

plot!(
    p_mcf,
    rho ./ r0sw,
    MCDOC_th;
    label="Theory (SW)",
    color=:black,
    linestyle=:dashdot,
    lw=1.5,
)

display(p_mcf)

# Show last realization irradiance
p_irr = heatmap(
    xn[1, :] .* 1e3,
    yn[:, 1] .* 1e3,
    abs2.(Uout);
    aspect_ratio=:equal,
    color=:inferno,
    title="Last Realization Irradiance",
    xlabel="x [mm]",
    ylabel="y [mm]",
)
display(p_irr)

p_phase = heatmap(
    xn[1, :] .* 1e3,
    yn[:, 1] .* 1e3,
    angle.(Uout);
    aspect_ratio=:equal,
    color=:phase, # A cyclic colormap is best for wrapped phase
    title="Turbulent Phase (Last Realization)",
    xlabel="x [mm]",
    ylabel="y [mm]",
    colorbar_title="Phase [rad]",
)

display(p_phase)
