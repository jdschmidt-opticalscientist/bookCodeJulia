using .OpticalWavePropSim
using Plots

# 1. Parameters
l0 = 0.0
L0 = Inf
nreals = 20

zt = vcat(0.0, z_planes)
Delta_z = diff(zt)

alpha = zt ./ zt[end]
delta_path = (1 .- alpha) .* delta1 .+ alpha .* deltan

# 2. Initialization
phz = zeros(N, N, nscr)
MCF2 = zeros(ComplexF64, N, N)
mask = circ.(xn, yn, D2)

# 3. Monte Carlo Loop
for idxreal in 1:nreals
    println("Realization $idxreal")
    
    for idxscr in 1:nscr
        phz_lo, phz_hi = ft_sh_phase_screen(r0scrn[idxscr], N, delta_path[idxscr], L0, l0)
        phz[:, :, idxscr] = phz_lo .+ phz_hi
    end
    
    # Applying phase screens and super-gaussian absorber
    # t is the sg_stack initialized in the vacuum script
    complex_screens = t .* exp.(im .* phz)
    
    # Propagate
    xn, yn, Uout = ang_spec_multi_prop(pt, wvl, delta1, deltan, z_planes, complex_screens)
    
    # Collimate
    Uout .*= exp.(-im * π / (wvl * R) .* (xn.^2 .+ yn.^2))
    
    # MCF accumulation
    MCF2 .+= corr2_ft(Uout, Uout, mask, deltan)
end

# 4. Final Calculation
mid = Int(N/2) + 1
MCDOC2 = abs.(MCF2) ./ abs(MCF2[mid, mid])

# Plotting the result
# heatmap(MCDOC2, title="Complex Degree of Coherence", c=:viridis)