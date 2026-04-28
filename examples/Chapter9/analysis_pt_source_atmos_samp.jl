using Plots, Printf

# analysis_pt_source_atmos_samp.m
c = 2;
D1p = D1 + c * wvl * Dz / r0sw;
D2p = D2 + c * wvl * Dz / r0sw;
delta1_vec = range(1e-6, 1.1 * wvl * Dz / D2p, length = 100);
deltan_vec = range(1e-6, 1.1 * wvl * Dz / D1p, length = 100);

# constraint 1
deltan_max = -D2p / D1p .* delta1_vec .+ wvl * Dz / D1p;

# constraint 3
d2min3 = (1 + Dz / R) .* delta1_vec .- wvl * Dz / D1p;
d2max3 = (1 + Dz / R) .* delta1_vec .+ wvl * Dz / D1p;

# [delta1 deltan] = meshgrid(delta1, deltan);
# In Julia, we can calculate N2 across the grid using broadcasting
N2 = [
    (wvl * Dz + D1p * dn + D2p * d1) / (2 * d1 * dn) for dn in deltan_vec, d1 in delta1_vec
];

# constraint 4
d1 = 10e-3;
d2 = 10e-3;
N = 512;
# d1*d2 * N / wvl (Matlab comment)
zmax = min(d1, d2)^2 * N / wvl;
nmin = ceil(Dz / zmax) + 1;

# --- Plotting Constraints 1, 2, & 3 ---
# Define the explicit levels we want to track (e.g., N=2^k to N=2^14)
# Make sure the values match the levels used in the contour
levels_to_track = 1:14

# Plot the filled contour with optimized colorbar
p = contourf(
    delta1_vec .* 1e3,
    deltan_vec .* 1e3,
    log2.(N2),
    levels = levels_to_track,
    color = :jet,
    clims = (minimum(levels_to_track), maximum(levels_to_track)),
    title = "Constraints 1, 2, & 3",
    xlabel = "delta_1 [mm]",
    ylabel = "delta_n [mm]",
    colorbar_title = "Required log2(N)",
    # --- Colorbar Customization to replace inline labels ---
    right_margin = 10Plots.mm, # Make space for the labels
    colorbar_ticks = (levels_to_track, ["$l" for l in levels_to_track]), # Force labels on k=1, 2, ... 14
)

# Add black lines on the contour levels for clarity
contour!(
    p,
    delta1_vec .* 1e3,
    deltan_vec .* 1e3,
    log2.(N2),
    levels = levels_to_track,
    color = :black,
    lw = 0.5,
)

# Overlay Constraint 1
plot!(
    p,
    delta1_vec .* 1e3,
    deltan_max .* 1e3,
    label = "Constraint 1",
    color = :black,
    lw = 2,
    linestyle = :dash,
)

# Overlay Constraint 3
plot!(
    p,
    delta1_vec .* 1e3,
    d2min3 .* 1e3,
    label = "Constraint 3 Min",
    color = :green,
    lw = 2,
    linestyle = :dashdot,
)
plot!(
    p,
    delta1_vec .* 1e3,
    d2max3 .* 1e3,
    label = "Constraint 3 Max",
    color = :blue,
    lw = 2,
    linestyle = :dashdot,
)

scatter!(
    p,
    [d1 * 1e3],
    [d2 * 1e3],
    marker = :x,
    ms = 8,
    mc = :white,
    msw = 3,
    label = "Chosen (d1, dn)",
)

# Set axis limits
xlims!(p, (0, maximum(delta1_vec) * 1e3))
ylims!(p, (0, maximum(deltan_vec) * 1e3))

display(p)

# --- Print key results to command line ---
println("-"^45)
println("Sampling Analysis Results")
println("-"^45)
@printf("Max propagation step (zmax): %.2f m\n", zmax)
@printf("Minimum screens needed:      %d\n", nmin)
@printf("Actual screens planned:      %d\n", nscr)
println("-"^45)

if nscr < nmin
    @warn "The number of screens (nscr=$nscr) is less than the required nmin ($nmin). Consider increasing nscr in the setup script."
end
