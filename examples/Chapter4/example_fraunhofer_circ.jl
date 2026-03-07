# using .OpticalWavePropSim
# using Plots

# example_fraunhofer_circ.jl
N = 512
L = 7.5e-3
d1 = L / N
D = 1e-3
wvl = 1e-6
k = 2 * π / wvl
Dz = 20.0

# Setup source plane
vec1 = collect(-N/2 : N/2-1) .* d1
x1 = ones(N) .* vec1'
y1 = vec1 .* ones(1, N)
Uin = circ(x1, y1, D)

# Fraunhofer propagation
Uout, x2, y2 = fraunhofer_prop(Uin, wvl, d1, Dz)

# Analytic result
Uout_th = (exp.(im * k / (2 * Dz) .* (x2.^2 .+ y2.^2)) 
           ./ (im * wvl * Dz) .* (D^2 * π / 4) 
           .* jinc.(D .* sqrt.(x2.^2 .+ y2.^2) ./ (wvl * Dz)))

# Visualization
# p1 = heatmap(x2[1,:], y2[:,1], abs2.(Uout), title="Fraunhofer Pattern")
# p2 = plot(x2[Int(N/2)+1, :], abs2.(Uout[Int(N/2)+1, :]), label="Numerical")
# plot!(p2, x2[Int(N/2)+1, :], abs2.(Uout_th[Int(N/2)+1, :]), label="Analytic")