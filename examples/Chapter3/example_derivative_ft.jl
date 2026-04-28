# example_derivative_ft.jl ---------------------------------------

using Plots
using OpticalWavePropSim

N = 64          # number of samples
L = 6.0         # grid size [m]
delta = L / N   # grid spacing [m]
x = ((-N / 2):(N / 2 - 1)) * delta

w = 3.0         # size of window (or region of interest) [m]
window = rect.(x / w)   # window function
g = (x .^ 5) .* window    # function

# discrete derivatives using FT property
gp_samp = real.(derivative_ft(g, delta, 1)) .* window
gpp_samp = real.(derivative_ft(g, delta, 2)) .* window

# analytic derivatives
gp = 5 * (x .^ 4) .* window
gpp = 20 * (x .^ 3) .* window

# Left Plot: First Derivative
p1 = plot(
    x,
    gp;
    linecolor=:red,
    lw=1.5,
    label="Analytic",
    title="First Derivative (5x⁴)",
    xlabel="x [m]",
    grid=true,
    linestyle=:dash,
)

scatter!(p1, x, gp_samp; marker=:x, mc=:blue, ms=4, label="FT Derivative")

# Right Plot: Second Derivative
p2 = plot(
    x,
    gpp;
    linecolor=:red,
    lw=1.5,
    label="Analytic",
    title="Second Derivative (20x³)",
    xlabel="x [m]",
    grid=true,
    linestyle=:dash,
)

scatter!(p2, x, gpp_samp; marker=:+, mc=:blue, ms=5, label="FT Derivative")

# Combine into a 1x2 layout
p_final = plot(p1, p2; layout=grid(1, 2), size=(1000, 450), margin=5Plots.mm)

display(p_final)
