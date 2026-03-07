using .OpticalWavePropSim
using Plots
gr()

# example_ft_sh_phase_screen.jl
D = 2.0
r0 = 0.1
N = 256
L0 = 100.0
l0 = 0.01

delta = D / N
x_vec = collect(-N/2 : N/2-1) .* delta

# Generate phase screen components
phz_lo, phz_hi = ft_sh_phase_screen(r0, N, delta, L0, l0)
phz = phz_lo .+ phz_hi

# Plotting
p1 = heatmap(x_vec, x_vec, phz_hi, title="FFT High", aspect_ratio=:equal, c=:viridis)
p2 = heatmap(x_vec, x_vec, phz_lo, title="SH Low", aspect_ratio=:equal, c=:viridis)
p3 = heatmap(x_vec, x_vec, phz, title="Total Phase", aspect_ratio=:equal, c=:viridis)

plot(p1, p2, p3, layout=(1,3), size=(1200, 400))