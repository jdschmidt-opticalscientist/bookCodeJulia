using JuMP, Ipopt, Printf

# example_pt_source_atmos_setup.jl
# determine geometry
D2 = 0.5;   # diameter of the observation aperture [m]
wvl = 1e-6;  # optical wavelength [m]
k = 2 * π / wvl; # optical wavenumber [rad/m]
Dz = 50e3;     # propagation distance [m]

# use sinc to model pt source
DROI = 4 * D2;  # diam of obs-plane region of interest [m]
D1 = wvl * Dz / DROI;    # width of central lobe [m]
R = Dz; # wavefront radius of curvature [m]

# atmospheric properties
Cn2 = 1e-16;    # structure parameter [m^-2/3], constant
# SW and PW coherence diameters [m]
r0sw = (0.423 * k^2 * Cn2 * 3 / 8 * Dz)^(-3 / 5);
r0pw = (0.423 * k^2 * Cn2 * Dz)^(-3 / 5);
p_vec = range(0, Dz, length = 1000);
dp = step(p_vec);
# log-amplitude variance
rytov =
    0.563 *
    k^(7 / 6) *
    sum(Cn2 .* max.(0, 1 .- p_vec ./ Dz) .^ (5 / 6) .* p_vec .^ (5 / 6) .* dp);

# screen properties
nscr = 11; # number of screens
A = zeros(2, nscr); # matrix
alpha = range(0, 1, length = nscr);
A[1, :] = collect(alpha) .^ (5 / 3);
A[2, :] = max.(0, 1 .- alpha) .^ (5 / 6) .* collect(alpha) .^ (5 / 6);
b = [r0sw^(-5 / 3), rytov / 1.33 * (k / Dz)^(5 / 6)];
# --- Optimization using JuMP & Ipopt ---
model = Model(Ipopt.Optimizer)
set_silent(model) # Keep the output clean
@variable(model, X[1:nscr] >= 0) # X >= 0 constraint
# Set upper bounds
rmax = 0.1; # maximum Rytov number per partial prop
for i = 1:nscr
    if A[2, i] == 0
        set_upper_bound(X[i], 50.0^(-5 / 3))
    else
        set_upper_bound(X[i], (rmax / 1.33 * (k / Dz)^(5 / 6)) / A[2, i])
    end
end
# Objective: Minimize sum of squared residuals
@objective(model, Min, sum((sum(A[j, i] * X[i] for i = 1:nscr) - b[j])^2 for j = 1:2))
optimize!(model)
X_opt = value.(X); # Extract optimized values
# check screen r0s
r0scrn = X_opt .^ (-3 / 5);
replace!(r0scrn, Inf => 1e6);

# check resulting r0sw & rytov
bp = A * X_opt;
sim_r0 = bp[1]^(-3 / 5);
sim_rytov = bp[2] * 1.33 * (Dz / k)^(5 / 6);

# Calculate Error
error_r0 = abs(r0sw - sim_r0) / r0sw;

# --- Formatted Console Output ---

# Individual Phase Screen Strengths Table
println("Individual Phase Screen Strengths:")
@printf("%-10s | %10s\n", "Screen #", "r0 [m]")
println("-"^23)
for i = 1:nscr
    # Print as 0-indexed to match textbook notation, using %g for formatting
    @printf("Screen %-3d | %10.4g\n", i - 1, r0scrn[i])
end
println("-"^23)
println()

# Summary Metrics Table
println("="^45)
@printf("%-15s | %12s | %12s\n", "Metric", "Theoretical", "Simulated")
println("-"^45)
@printf("%-15s | %12.4f | %12.4f\n", "r0 (SW)", r0sw, sim_r0)
@printf("%-15s | %12.4e | %12.4e\n", "Rytov Var.", rytov, sim_rytov)
# Using %% to print the literal percent sign
@printf("%-15s | %12s | %11.2f%%\n", "r0 Error", "-", error_r0 * 100)
println("="^45)
println()
