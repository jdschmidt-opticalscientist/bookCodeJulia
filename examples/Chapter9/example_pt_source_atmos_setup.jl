using JuMP, Ipopt

# example_pt_source_atmos_setup.jl
D2 = 0.5
wvl = 1e-6
k = 2 * π / wvl
Dz = 50e3

DROI = 4 * D2
D1 = wvl * Dz / DROI

Cn2 = 1e-16
r0sw = (0.423 * k^2 * Cn2 * 3/8 * Dz)^(-3/5)

p = range(0, Dz, length=1000)
dp = step(p)
rytov = 0.563 * k^(7/6) * sum(Cn2 .* (1 .- p./Dz).^(5/6) .* p.^(5/6) .* dp)

nscr = 11
alpha = range(0, 1, length=nscr)
A = zeros(2, nscr)
A[1, :] = collect(alpha).^(5/3)
A[2, :] = (1 .- alpha).^(5/6) .* collect(alpha).^(5/6)
b = [r0sw^(-5/3), rytov / 1.33 * (k/Dz)^(5/6)]

# Optimization using JuMP
model = Model(Ipopt.Optimizer)
set_silent(model)

@variable(model, X[i=1:nscr] >= 0)

# Set upper bounds
rmax = 0.1
for i in 1:nscr
    if A[2, i] == 0
        set_upper_bound(X[i], 50.0^(-5/3))
    else
        set_upper_bound(X[i], (rmax / 1.33 * (k/Dz)^(5/6)) / A[2, i])
    end
end

@objective(model, Min, sum((A * X - b).^2))
optimize!(model)

X_opt = value.(X)
r0scrn = X_opt.^(-3/5)
replace!(r0scrn, Inf => 1e6)

# Check results
bp = A * X_opt
println("Target [r0sw, rytov]: ", [r0sw, rytov])
println("Actual [r0sw, rytov]: ", [bp[1]^(-3/5), bp[2]*1.33*(Dz/k)^(5/6)])