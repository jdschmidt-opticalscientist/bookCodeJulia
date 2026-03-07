# OpticalWavePropSim.jl
# Copyright (c) 2026, Jason D. Schmidt
# This source code is licensed under the BSD 3-Clause License 
# found in the LICENSE file in the root directory of this source tree.

module OpticalWavePropSim

using FFTW
using SpecialFunctions # For besselj, fresnelc, fresnels, and gamma
using Statistics       # For mean

export ft, ift, ft2, ift2, rect, tri, jinc, circ,
    ang_spec_multi_prop, ang_spec_prop, corr2_ft,
    ft_phase_screen, ft_sh_phase_screen, zernike

# --- Fourier Transform Utilities ---

function ft(g, delta)
    return fftshift(fft(fftshift(g))) * delta
end

function ift(G, delta_f)
    return ifftshift(ifft(ifftshift(G))) * (length(G) * delta_f)
end

function ft2(g, delta)
    return fftshift(fft2(fftshift(g))) * (delta^2)
end

function ift2(G, delta_f)
    N = size(G, 1)
    return ifftshift(ifft2(ifftshift(G))) * (N * delta_f)^2
end

# --- Basic Optical Functions ---

function rect(x, D=1.0)
    x_abs = abs.(x)
    y = float.(x_abs .< D / 2)
    y[x_abs.==D/2] .= 0.5
    return y
end

function tri(t)
    t_abs = abs.(t)
    y = zeros(size(t))
    idx = t_abs .< 1.0
    y[idx] = 1.0 .- t_abs[idx]
    return y
end

function jinc(x)
    y = ones(size(x))
    idx = x .!= 0
    y[idx] = 2.0 .* besselj.(1, π .* x[idx]) ./ (π .* x[idx])
    return y
end

function circ(x, y, D=1.0)
    r = sqrt.(x .^ 2 .+ y .^ 2)
    return float.(r .<= D / 2)
end

# --- Propagation Algorithms ---

function ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t)
    N = size(Uin, 1)
    vec = collect(-N/2:N/2-1)
    nx = vec' .* ones(N)
    ny = ones(N, 1) .* vec
    k = 2 * π / wvl

    nsq = nx .^ 2 .+ ny .^ 2
    sg = exp.(-nsq .^ 8 ./ (0.47 * N)^16)

    z_arr = vcat(0.0, z)
    n = length(z_arr)
    Delta_z = diff(z_arr)

    alpha = z_arr ./ z_arr[end]
    delta = (1 .- alpha) .* delta1 .+ alpha .* deltan
    m = delta[2:end] ./ delta[1:end-1]

    r1sq = (nx .* delta[1]) .^ 2 .+ (ny .* delta[1]) .^ 2
    Q1 = exp.(im * k / 2 * (1 - m[1]) / Delta_z[1] .* r1sq)
    Uin = Uin .* Q1 .* t[:, :, 1]

    curr_U = Uin
    local dz = 0.0
    for idx in 1:(n-1)
        deltaf = 1 / (N * delta[idx])
        fsq = (nx .* deltaf) .^ 2 .+ (ny .* deltaf) .^ 2
        dz = Delta_z[idx]

        Q2 = exp.(-im * π^2 * 2 * dz / m[idx] / k .* fsq)
        # Proper propagation step
        curr_U = sg .* t[:, :, idx+1] .* ift2(Q2 .* ft2(curr_U ./ m[idx], delta[idx]), deltaf)
    end

    xn, yn = nx .* delta[end], ny .* delta[end]
    rnsq = xn .^ 2 .+ yn .^ 2
    Q3 = exp.(im * k / 2 * (m[end] - 1) / (m[end] * dz) .* rnsq)
    return xn, yn, Q3 .* curr_U
end

function ang_spec_prop(Uin, wvl, d1, d2, Dz)
    N = size(Uin, 1)
    k = 2 * π / wvl
    vec = collect(-N/2:N/2-1)
    nx = vec' .* ones(N)
    ny = ones(N, 1) .* vec

    df1 = 1 / (N * d1)
    m = d2 / d1

    Q1 = exp.(im * k / 2 * (1 - m) / Dz .* ((nx .* d1) .^ 2 .+ (ny .* d1) .^ 2))
    Q2 = exp.(-im * π^2 * 2 * Dz / m / k .* ((nx .* df1) .^ 2 .+ (ny .* df1) .^ 2))
    Q3 = exp.(im * k / 2 * (m - 1) / (m * Dz) .* ((nx .* d2) .^ 2 .+ (ny .* d2) .^ 2))

    Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin ./ m, d1), df1)
    return nx .* d2, ny .* d2, Uout
end

# --- Turbulence & Analysis ---

function ft_phase_screen(r0, N, delta, L0, l0)
    del_f = 1 / (N * delta)
    vec = collect(-N/2:N/2-1) .* del_f
    fx = vec' .* ones(N)
    fy = ones(N, 1) .* vec

    f = sqrt.(fx .^ 2 .+ fy .^ 2)
    fm = 5.92 / (l0 * 2 * π)
    f0 = 1 / L0

    PSD_phi = 0.023 * r0^(-5 / 3) .* exp.(-(f ./ fm) .^ 2) ./ (f .^ 2 .+ f0^2) .^ (11 / 6)
    PSD_phi[Int(N / 2)+1, Int(N / 2)+1] = 0

    cn = (randn(N, N) .+ im .* randn(N, N)) .* sqrt.(PSD_phi) .* del_f
    return real(ift2(cn, 1.0))
end

function ft_sh_phase_screen(r0, N, delta, L0, l0)
    D = N * delta
    phz_hi = ft_phase_screen(r0, N, delta, L0, l0)

    vec = collect(-N/2:N/2-1) .* delta
    x = vec' .* ones(N)
    y = ones(N, 1) .* vec
    phz_lo = zeros(ComplexF64, N, N)

    for p in 1:3
        del_f = 1 / (3^p * D)
        vec_f = collect(-1:1) .* del_f
        fx = vec_f' .* ones(3)
        fy = ones(3, 1) .* vec_f
        f = sqrt.(fx .^ 2 .+ fy .^ 2)

        fm = 5.92 / (l0 * 2 * π)
        f0 = 1 / L0
        PSD_phi = 0.023 * r0^(-5 / 3) .* exp.(-(f ./ fm) .^ 2) ./ (f .^ 2 .+ f0^2) .^ (11 / 6)
        PSD_phi[2, 2] = 0

        cn = (randn(3, 3) .+ im .* randn(3, 3)) .* sqrt.(PSD_phi) .* del_f

        SH = zeros(ComplexF64, N, N)
        for ii in 1:9
            SH .+= cn[ii] .* exp.(im * 2 * π .* (fx[ii] .* x .+ fy[ii] .* y))
        end
        phz_lo .+= SH
    end

    phz_lo_real = real(phz_lo)
    phz_lo_real .-= mean(phz_lo_real)
    return phz_lo_real, phz_hi
end

function corr2_ft(u1, u2, mask, delta)
    N = size(u1, 1)
    delta_f = 1 / (N * delta)

    U1 = ft2(u1 .* mask, delta)
    U2 = ft2(u2 .* mask, delta)
    U12corr = ift2(conj.(U1) .* U2, delta_f)

    areamask = sum(mask) * delta^2
    maskcorr = ift2(abs.(ft2(mask, delta)) .^ 2, delta_f) / areamask

    idx = maskcorr .>= (delta^2 / areamask)
    c = zeros(ComplexF64, N, N)
    c[idx] = U12corr[idx] ./ maskcorr[idx]

    return c, idx
end

function _zrf(n, m, r)
    R = zeros(size(r))
    for s in 0:Int((n - m) / 2)
        num = (-1)^s * gamma(n - s + 1)
        denom = gamma(s + 1) * gamma((n + m) / 2 - s + 1) * gamma((n - m) / 2 - s + 1)
        R .+= (num / denom) .* r .^ (n - 2 * s)
    end
    return R
end

function zernike(i, r, theta, zernike_index)
    n = Int(zernike_index[i, 1])
    m = Int(zernike_index[i, 2])
    if m == 0
        return sqrt(n + 1) .* _zrf(n, 0, r)
    else
        if i % 2 == 0
            return sqrt(2 * (n + 1)) .* _zrf(n, m, r) .* cos.(m .* theta)
        else
            return sqrt(2 * (n + 1)) .* _zrf(n, m, r) .* sin.(m .* theta)
        end
    end
end

end # module