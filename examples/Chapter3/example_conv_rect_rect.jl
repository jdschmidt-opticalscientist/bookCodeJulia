# example_conv_rect_rect.jl --------------------------------------

using Plots
using OpticalWavePropSim

N = 64          # number of samples
L = 8.0         # grid size [m]
delta = L / N   # sample spacing [m]
x = (-N/2:N/2-1) * delta

w = 2.0         # width of rectangle
A = rect.(x / w) # signal A
B = A           # signal B

# Perform discrete convolution
C = myconv(A, B, delta)

# Continuous convolution
C_cont = w * tri.(x / w)

# Plot (a): Signal A
p1 = plot(x, A,
    seriestype=:scatter,
    marker=:x, mc=:blue, ms=4,
    title="(a) A(x)",
    ylabel="Amplitude",
    grid=true, label=false)

# Plot (b): Signal B
p2 = plot(x, B,
    seriestype=:scatter,
    marker=:x, mc=:green, ms=4,
    title="(b) B(x)",
    grid=true, label=false)

# Plot (c): Convolution Result (Analytic vs Numerical)
p3 = plot(x, C_cont,
    color=:red, marker=:square, ms=3, lw=1,
    label="Analytic",
    title="(c) A(x) * B(x)",
    xlabel="x [m]", ylabel="Amplitude",
    grid=true)

scatter!(p3, x, real.(C),
    marker=:x, mc=:blue, ms=4,
    label="Numerical")

# Combine into a layout: 2 rows. Top row has 2 plots, bottom row has 1 spanning both.
# @layout [a b; c] creates exactly the Python subplot(2,2,1), (2,2,2), and (2,1,2)
l = @layout [grid(1, 2); c]
p_final = plot(p1, p2, p3, layout=l, size=(800, 600), margin=5Plots.mm)

display(p_final)