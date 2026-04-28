using OpticalWavePropSim
using Plots

# example_ft_sh_phase_screen.jl
D = 2;    # length of one side of square phase screen [m]
r0 = 0.1; # coherence diameter [m]
N = 256;  # number of grid points per side
L0 = 100; # outer scale [m]
l0 = 0.01;# inner scale [m]

delta = D / N;    # grid spacing [m]
# spatial grid
x = collect((-N / 2):(N / 2 - 1)) .* delta;
y = x;
# generate a random draw of an atmospheric phase screen
phz_lo, phz_hi = ft_sh_phase_screen(r0, N, delta, L0, l0);
phz = phz_lo .+ phz_hi;

# --- Visualization ---
p1 = heatmap(
    x,
    y,
    phz_hi;
    aspect_ratio=:equal,
    c=:viridis,
    title="High Freq (FFT only)",
    xlabel="x [m]",
    ylabel="y [m]",
)

p2 = heatmap(
    x,
    y,
    phz_lo;
    aspect_ratio=:equal,
    c=:viridis,
    title="Low Freq (Subharmonics)",
    xlabel="x [m]",
    ylabel="y [m]",
)

p3 = heatmap(
    x,
    y,
    phz;
    aspect_ratio=:equal,
    c=:viridis,
    title="Combined Phase Screen",
    xlabel="x [m]",
    ylabel="y [m]",
)

# Combine and display
final_plot = plot(p1, p2, p3; layout=(1, 3), size=(1200, 400), margin=5Plots.mm)
display(final_plot)
