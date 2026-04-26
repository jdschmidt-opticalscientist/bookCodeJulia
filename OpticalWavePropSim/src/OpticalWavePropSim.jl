# OpticalWavePropSim.jl
# Copyright (c) 2026, Jason D. Schmidt
# This source code is licensed under the BSD 3-Clause License 
# found in the LICENSE file in the root directory of this source tree.

module OpticalWavePropSim

using SpecialFunctions
using FFTW
using Statistics: mean

export ft, ift, ft2, ift2, myconv, myconv2, rect, tri, jinc, circ,
    derivative_ft, gradient_ft, corr2_ft, str_fcn2_ft, zernike,
    fraunhofer_prop, ang_spec_multi_prop, ang_spec_prop, ang_spec_propABCD,
    one_step_prop, two_step_prop, ang_spec_multi_prop_vac,
    ft_phase_screen, ft_sh_phase_screen, fresnel_prop_square_ap


# --- Fourier Transform Utilities ---

function ft(g, delta)
    return fftshift(fft(fftshift(g))) * delta
end

function ift(G, delta_f)
    return ifftshift(ifft(ifftshift(G))) * (length(G) * delta_f)
end

function ft2(g, delta)
    G = fftshift(fft(fftshift(g))) * delta^2
    return G
end

function ift2(G, delta_f)
    N = size(G, 1)
    g = ifftshift(ifft(ifftshift(G))) * (N * delta_f)^2
    return g
end

function myconv(A, B, delta)
    N = length(A)
    C = ift(ft(A, delta) .* ft(B, delta), 1 / (N * delta))
    return C
end

function myconv2(A, B, delta)
    N = size(A, 1)
    df = 1 / (N * delta)
    C = ift2(ft2(A, delta) .* ft2(B, delta), df)
    return C
end
# --- Basic Optical Functions ---

function rect(x, D=1.0)
    # logic for a single scalar value
    x_abs = abs(x)
    if x_abs < D / 2
        return 1.0
    elseif x_abs == D / 2
        return 0.5
    else
        return 0.0
    end
end

function tri(t)
    t_abs = abs(t)
    if t_abs < 1.0
        return 1.0 - t_abs
    else
        return 0.0
    end
end

function jinc(x)
    # Handle the limit at x = 0 to avoid division by zero
    if x == 0
        return 1.0
    else
        # Bessel function of the first kind, order 1: J1(pi*x) / (2*pi*x)
        # Note: Using SpecialFunctions.besselj or similar
        return 2 * besselj1(π * x) / (π * x)
    end
end

function circ(x, y, D=1.0)
    r = sqrt.(x .^ 2 .+ y .^ 2)
    return float.(r .<= D / 2)
end

# --- Propagation Algorithms ---

function fraunhofer_prop(Uin, wvl, d1, Dz)
    # function [Uout x2 y2] = fraunhofer_prop(Uin, wvl, d1, Dz)
    N = size(Uin, 1)    # assume square grid
    k = 2 * π / wvl     # optical wavevector

    # Frequency-domain grid spacing
    f_vec = (-N/2:N/2-1) / (N * d1)

    # Observation-plane coordinates: x2 = wvl * Dz * fx
    v = wvl * Dz * f_vec
    x2 = repeat(v', N, 1)
    y2 = repeat(v, 1, N)

    # Fraunhofer propagation formula:
    # Uout(x2,y2) = exp(i*k/(2*z)*(x2^2+y2^2)) / (i*wvl*z) * FT{Uin}
    # Note: ft2 handles the delta^2 scaling internally
    Uout = (exp.(im * k / (2 * Dz) .* (x2 .^ 2 .+ y2 .^ 2))
            /
            (im * wvl * Dz) .* ft2(Uin, d1))

    return Uout, x2, y2
end

# Use the Faddeeva function from SpecialFunctions to define the Fresnel integrals
function fresnel_integrals(z)
    # The Faddeeva function w(z) relates to Fresnel integrals
    # Fresnel integrals S(z) and C(z)
    # Note: SpecialFunctions.faddeeva(z) is available
    w = faddeeva((1 + im) * sqrt(π) * z / 2)
    val = (1 + im) / 2 * (1 - exp(-im * π * z^2 / 2) * w)
    return real(val), imag(val) # Returns (C, S)
end

# Inside OpticalWavePropSim.jl (below your imports)
import SpecialFunctions: erf

# Fresnel integrals defined via the complex Error Function
function fresnelc(z)
    val = (1 + im) / 2 * erf((1 - im) * sqrt(π) / 2 * z)
    return real(val)
end

function fresnels(z)
    val = (1 + im) / 2 * erf((1 - im) * sqrt(π) / 2 * z)
    return imag(val)
end

function fresnel_prop_square_ap(x2, y2, D1, wvl, Dz)
    # function U = fresnel_prop_square_ap(x2, y2, D1, wvl, Dz)
    N_F = (D1 / 2)^2 / (wvl * Dz) # Fresnel number

    # substitutions
    bigX = x2 ./ sqrt(wvl * Dz)
    bigY = y2 ./ sqrt(wvl * Dz)

    alpha1 = -sqrt(2) .* (sqrt(N_F) .+ bigX)
    alpha2 = sqrt(2) .* (sqrt(N_F) .- bigX)
    beta1 = -sqrt(2) .* (sqrt(N_F) .+ bigY)
    beta2 = sqrt(2) .* (sqrt(N_F) .- bigY)

    # Fresnel sine and cosine integrals
    ca1 = fresnelc.(alpha1)
    sa1 = fresnels.(alpha1)
    ca2 = fresnelc.(alpha2)
    sa2 = fresnels.(alpha2)

    cb1 = fresnelc.(beta1)
    sb1 = fresnels.(beta1)
    cb2 = fresnelc.(beta2)
    sb2 = fresnels.(beta2)

    # observation-plane field
    U = 1 / (2 * im) .* ((ca2 .- ca1) .+ im .* (sa2 .- sa1)) .* ((cb2 .- cb1) .+ im .* (sb2 .- sb1))

    return U
end

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

function ang_spec_propABCD(Uin, wvl, d1, d2, ABCD)
    # function [x2, y2, Uout] = ang_spec_propABCD(Uin, wvl, d1, d2, ABCD)

    N = size(Uin, 1)   # assume square grid
    # source-plane coordinates
    vec1 = collect(-N/2:N/2-1) .* d1
    x1 = repeat(vec1', N, 1)
    y1 = repeat(vec1, 1, N)
    r1sq = x1 .^ 2 .+ y1 .^ 2
    # spatial frequencies (of source plane)
    df1 = 1 / (N * d1)
    vecf = collect(-N/2:N/2-1) .* df1
    fX = repeat(vecf', N, 1)
    fY = repeat(vecf, 1, N)
    fsq = fX .^ 2 .+ fY .^ 2
    # scaling parameter
    m = d2 / d1
    # observation-plane coordinates
    vec2 = collect(-N/2:N/2-1) .* d2
    x2 = repeat(vec2', N, 1)
    y2 = repeat(vec2, 1, N)
    r2sq = x2 .^ 2 .+ y2 .^ 2
    # optical system matrix (1-based indexing)
    A = ABCD[1, 1]
    B = ABCD[1, 2]
    D = ABCD[2, 2]
    # quadratic phase factors
    Q1 = exp.(im * π / (wvl * B) * (A - m) .* r1sq)
    Q2 = exp.(-im * π * wvl * B / m .* fsq)
    Q3 = exp.(im * π / (wvl * B) * (D - 1 / m) .* r2sq)
    # compute the propagated field
    Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin ./ m, d1), df1)
    return x2, y2, Uout
end

function one_step_prop(Uin, wvl, d1, Dz)
    # function [x2 y2 Uout] ...
    #     = one_step_prop(Uin, wvl, d1, Dz)
    N = size(Uin, 1)   # assume square grid
    k = 2 * π / wvl    # optical wavevector
    # source-plane coordinates
    vec1 = collect(-N/2:1:N/2-1) .* d1
    x1 = repeat(vec1', N, 1)
    y1 = repeat(vec1, 1, N)
    # observation-plane coordinates
    vec2 = collect(-N/2:N/2-1) ./ (N * d1) .* wvl .* Dz
    x2 = repeat(vec2', N, 1)
    y2 = repeat(vec2, 1, N)
    # evaluate the Fresnel-Kirchhoff integral
    Uout = 1 / (im * wvl * Dz) .* exp.(im * k / (2 * Dz) .* (x2 .^ 2 .+ y2 .^ 2)) .* ft2(Uin .* exp.(im * k / (2 * Dz) .* (x1 .^ 2 .+ y1 .^ 2)), d1)

    return x2, y2, Uout
end

function two_step_prop(Uin, wvl, d1, d2, Dz)
    # function [x2 y2 Uout] ...
    #     = two_step_prop(Uin, wvl, d1, d2, Dz)
    N = size(Uin, 1)   # number of grid points
    k = 2 * π / wvl    # optical wavevector
    # source-plane coordinates
    vec1 = collect(-N/2:1:N/2-1) .* d1
    x1 = repeat(vec1', N, 1)
    y1 = repeat(vec1, 1, N)
    # magnification
    m = d2 / d1
    # intermediate plane
    Dz1 = Dz / (1 - m)  # propagation distance
    d1a = wvl * abs(Dz1) / (N * d1)    # coordinates
    vec1a = collect(-N/2:N/2-1) .* d1a
    x1a = repeat(vec1a', N, 1)
    y1a = repeat(vec1a, 1, N)
    # evaluate the Fresnel-Kirchhoff integral for the intermediate plane
    Uitm = 1 / (im * wvl * Dz1) .* exp.(im * k / (2 * Dz1) .* (x1a .^ 2 .+ y1a .^ 2)) .* ft2(Uin .* exp.(im * k / (2 * Dz1) .* (x1 .^ 2 .+ y1 .^ 2)), d1)
    # observation plane
    Dz2 = Dz - Dz1  # propagation distance
    # coordinates
    vec2 = collect(-N/2:N/2-1) .* d2
    x2 = repeat(vec2', N, 1)
    y2 = repeat(vec2, 1, N)
    # evaluate the Fresnel diffraction integral for the observation plane
    Uout = 1 / (im * wvl * Dz2) .* exp.(im * k / (2 * Dz2) .* (x2 .^ 2 .+ y2 .^ 2)) .* ft2(Uitm .* exp.(im * k / (2 * Dz2) .* (x1a .^ 2 .+ y1a .^ 2)), d1a)

    return x2, y2, Uout
end

function ang_spec_multi_prop_vac(Uin, wvl, delta1, deltan, z)
    # function [xn yn Uout] = ang_spec_multi_prop_vac ...
    #     (Uin, wvl, delta1, deltan, z)
    N = size(Uin, 1)   # number of grid points
    v = collect(-N/2:1:N/2-1)
    nx = repeat(v', N, 1)
    ny = repeat(v, 1, N)
    k = 2 * π / wvl    # optical wavevector
    # super-Gaussian absorbing boundary
    nsq = nx .^ 2 .+ ny .^ 2
    w = 0.47 * N
    sg = exp.(-nsq .^ 8 ./ w^16) # element-wise exponentiation

    # propagation plane locations
    z_planes = [0.0; vec(collect(z))]
    n = length(z_planes)
    # propagation distances
    Delta_z = diff(z_planes)
    # grid spacings
    alpha = z_planes ./ z_planes[n]
    delta = (1 .- alpha) .* delta1 .+ alpha .* deltan
    m = delta[2:n] ./ delta[1:n-1]
    x1 = nx .* delta[1]
    y1 = ny .* delta[1]
    r1sq = x1 .^ 2 .+ y1 .^ 2

    Q1 = exp.(im * k / 2 * (1 - m[1]) / Delta_z[1] .* r1sq)
    Uin = Uin .* Q1
    local Z = 0.0
    for idx = 1:n-1
        # spatial frequencies (of i^th plane)
        deltaf = 1 / (N * delta[idx])
        fX = nx .* deltaf
        fY = ny .* deltaf
        fsq = fX .^ 2 .+ fY .^ 2
        Z = Delta_z[idx]   # propagation distance
        # quadratic phase factor
        Q2 = exp.(-im * π^2 * 2 * Z / m[idx] / k .* fsq)
        # compute the propagated field
        Uin = sg .* ift2(Q2 .* ft2(Uin ./ m[idx], delta[idx]), deltaf)
    end

    # observation-plane coordinates
    xn = nx .* delta[n]
    yn = ny .* delta[n]
    rnsq = xn .^ 2 .+ yn .^ 2
    Q3 = exp.(im * k / 2 * (m[n-1] - 1) / (m[n-1] * Z) .* rnsq)
    Uout = Q3 .* Uin

    return xn, yn, Uout
end

# --- Turbulence & Analysis ---

function ft_phase_screen(r0, N, delta, L0, l0)

    # setup the PSD
    del_f = 1 / (N * delta)   # frequency grid spacing [1/m]
    v_f = collect(-N/2:N/2-1) .* del_f
    # frequency grid [1/m]
    fx = repeat(v_f', N, 1)
    fy = repeat(v_f, 1, N)
    f = sqrt.(fx .^ 2 .+ fy .^ 2)  # polar grid
    fm = 5.92 / (l0 * 2 * π) # inner scale frequency [1/m]
    f0 = 1 / L0           # outer scale frequency [1/m]
    # modified von Karman atmospheric phase PSD
    PSD_phi = 0.023 * r0^(-5 / 3) .* exp.(-(f ./ fm) .^ 2) ./ (f .^ 2 .+ f0^2) .^ (11 / 6)
    PSD_phi[Int(N / 2)+1, Int(N / 2)+1] = 0
    # random draws of Fourier coefficients
    cn = (randn(N, N) .+ im .* randn(N, N)) .* sqrt.(PSD_phi) .* del_f
    # synthesize the phase screen
    return real(ift2(cn, 1.0))
end

function ft_sh_phase_screen(r0, N, delta, L0, l0)
    D = N * delta
    # high-frequency screen from FFT method
    phz_hi = ft_phase_screen(r0, N, delta, L0, l0)
    # spatial grid [m]
    v_space = collect(-N/2:N/2-1) .* delta
    x = repeat(v_space', N, 1)
    y = repeat(v_space, 1, N)
    # initialize low-freq screen
    phz_lo = zeros(ComplexF64, N, N)
    # loop over frequency grids with spacing 1/(3^p*L)
    for p in 1:3
        # setup the PSD
        del_f = 1 / (3^p * D)
        v_f = collect(-1:1) .* del_f
        # frequency grid [1/m]
        fx = repeat(v_f', 3, 1)
        fy = repeat(v_f, 1, 3)
        f = sqrt.(fx .^ 2 .+ fy .^ 2)  # polar grid
        fm = 5.92 / (l0 * 2 * π) # inner scale frequency [1/m]
        f0 = 1 / L0           # outer scale frequency [1/m]
        # modified von Karman atmospheric phase PSD
        PSD_phi = 0.023 * r0^(-5 / 3) .* exp.(-(f ./ fm) .^ 2) ./ (f .^ 2 .+ f0^2) .^ (11 / 6)
        PSD_phi[2, 2] = 0
        # random draws of Fourier coefficients
        cn = (randn(3, 3) .+ im .* randn(3, 3)) .* sqrt.(PSD_phi) .* del_f
        SH = zeros(ComplexF64, N, N)
        # loop over frequencies on this grid
        for ii in 1:9
            # Indexing into 3x3 grids fx and fy is now safe
            SH .+= cn[ii] .* exp.(im * 2 * π .* (fx[ii] .* x .+ fy[ii] .* y))
        end
        phz_lo .+= SH   # accumulate subharmonics
    end
    phz_lo_real = real(phz_lo)
    phz_lo_real .-= mean(phz_lo_real)
    return phz_lo_real, phz_hi
end

function corr2_ft(u1, u2, mask, delta)
    N = size(u1, 1)
    delta_f = 1 / (N * delta)

    # DFTs of masked signals
    U1 = ft2(u1 .* mask, delta)
    U2 = ft2(u2 .* mask, delta)

    # Cross-correlation in frequency domain
    U12corr = ift2(conj.(U1) .* U2, delta_f)

    # Area and mask normalization
    areamask = sum(mask) * delta^2
    # Mask autocorrelation
    maskcorr = ift2(abs.(ft2(mask, delta)) .^ 2, delta_f) / areamask

    # Logical indexing for valid mask overlap
    idx = real.(maskcorr) .>= (delta^2 / areamask)

    # Pre-allocate complex result array
    c = zeros(ComplexF64, N, N)

    # Perform division only where mask overlap is sufficient
    # We use .= to broadcast the assignment into the subset of c
    c[idx] .= U12corr[idx] ./ maskcorr[idx]

    return c, idx
end

function derivative_ft(g, delta, n)
    # function der = derivative_ft(g, delta, n)
    N = length(g)
    # frequency domain grid spacing [1/m]
    F = 1 / (N * delta)
    # frequency values range
    f_X = (-N/2:N/2-1) * F

    # Fourier derivative property: F{g^(n)} = (i*2*pi*f)^n * F{g}
    # Note the use of .^ and .* for element-wise operations
    der = ift((im * 2 * π * f_X) .^ n .* ft(g, delta), F)

    return der
end

function gradient_ft(g, delta)
    # function [gx gy] = gradient_ft(g, delta)
    N = size(g, 1)   # number of samples per side
    # grid spacing in the frequency domain [1/m]
    F = 1 / (N * delta)

    # Create frequency vectors and 2D grids
    f_vec = (-N/2:N/2-1) * F
    fX = repeat(f_vec', N, 1)
    fY = repeat(f_vec, 1, N)

    # Pre-calculate the FT of the input to avoid redundant transforms
    G = ft2(g, delta)

    # Fourier Gradient Property: F{∂g/∂x} = (i*2*pi*fX) * F{g}
    gx = ift2(im * 2 * π * fX .* G, F)
    gy = ift2(im * 2 * π * fY .* G, F)

    return gx, gy
end

function str_fcn2_ft(ph, mask, delta)
    # ph: phase screen matrix
    # mask: aperture mask matrix
    # delta: grid spacing [m]

    N = size(ph, 1)
    ph_masked = ph .* mask

    # Fourier Transforms
    P = ft2(ph_masked, delta)
    S = ft2(ph_masked .^ 2, delta)
    W = ft2(mask, delta)

    delta_f = 1 / (N * delta)

    # Weighting and intermediate terms
    # Using conj. for element-wise conjugation
    w2 = ift2(W .* conj.(W), delta_f)

    # Phase structure function calculation
    # We take the real part of the numerator terms as per the mathematical definition
    num = real.(S .* conj.(W)) .- abs.(P) .^ 2
    D = 2 * ift2(num, delta_f) ./ w2

    # Mask normalization and safety indexing (from your Python version)
    areamask = sum(mask) * delta^2
    maskcorr = ift2(abs.(ft2(mask, delta)) .^ 2, delta_f) / areamask

    # Avoid division by zero/noise in regions outside the mask overlap
    idx = real.(maskcorr) .>= (delta^2 / areamask)

    # Final result: ensure it's real and zeroed out where invalid
    D_real = zeros(Float64, N, N)
    D_real[idx] .= real.(D[idx])

    return D_real, idx
end

function _zrf(n, m, r)
    R = zeros(size(r))
    # n-m must be even and non-negative
    if (n - m) % 2 != 0 || n < m
        return R
    end

    for s in 0:Int((n - m) / 2)
        num = (-1)^s * gamma(n - s + 1)
        denom = (gamma(s + 1) * gamma((n + m) / 2 - s + 1) * gamma((n - m) / 2 - s + 1))
        # Use .^ for element-wise power on the radius matrix
        R .+= (num / denom) .* r .^ (n - 2 * s)
    end
    return R
end

function noll_to_nm(j)
    # Corrected Noll to (n, m) mapping
    n = Int(floor((sqrt(8 * j - 7) - 1) / 2))
    m_abs = j - n * (n + 1) ÷ 2 - 1

    if n % 2 == 0
        m = 2 * (m_abs ÷ 2)
    else
        m = 2 * ((m_abs + 1) ÷ 2) - 1
    end

    if j % 2 > 0
        m = -m
    end
    return n, m
end

function zernike(j, r, theta, z_map=nothing)
    # Use mapping if provided, otherwise default to Noll
    if !isnothing(z_map) && haskey(z_map, j)
        n, m_val = z_map[j]
    else
        n, m_val = noll_to_nm(j)
    end

    m = abs(m_val)
    # Epsilon boundary check using Julia's ternary broadcasting
    r_safe = [(val <= 1.000001) ? val : 0.0 for val in r]

    if m == 0
        return sqrt(n + 1) .* _zrf(n, 0, r_safe)
    elseif j % 2 == 0
        return sqrt(2 * (n + 1)) .* _zrf(n, m, r_safe) .* cos.(m .* theta)
    else
        return sqrt(2 * (n + 1)) .* _zrf(n, m, r_safe) .* sin.(m .* theta)
    end
end

end # module