using FFTW, Plots, SpecialFunctions

include("propagar.jl")

# ******* Definicion de parametros fisicos ********
w0 = 0.5e-3          # Cintura del haz Gaussiano [m]
λ = 633e-9           # longitud de onda [m]
k = 2π/λ             # numero de onda [1/m]
zR = k * w0^2 / 2    # Distancia de Rayleigh [m]
z = 1 * zR           # Distancia de propagacion maxima

# ******* Definicion de parametros numericos *******
N = 2^9
NV = (-N/2:1:N/2-1) * 2/N

L = 4 * w0
dx = 2 * L / N
kmax = π / dx
X, Y = meshgrid(NV .* L, NV .* L)
Kx, Ky = meshgrid(NV .* kmax, NV .* kmax)

r = sqrt.(X.^2 .+ Y.^2)
kt = sqrt.(Kx.^2 .+ Ky.^2)

nz = 300
dz = z / nz

# ******************** Perfil Inicial ********************
flens = 0.4 * zR
Tlens = 1   #exp(-1im * k / (2 * flens) * r.^2)
k1 = 10000
f = exp.(-r.^2/w0^2) .* besselj0.(k1 * r) .* Tlens

U0 = f
Ur = zeros(Complex{Float64}, N, nz+1)
Ur[:, 1] = U0[:, Int(N/2)+1]

# Plotting the initial field
heatmap(X[Int(N/2)+1, :]./w0, Y[:, Int(N/2)+1]./w0, abs.(U0), title="Campo U₀", xlabel="x/w₀", ylabel="y/w₀")

# ****************************** Propagador ******************************
Prop = exp.(-1im * 0.5 * dz * (kt.^2) / k)

F = fftshift(fft(U0))
for ii=1:nz
    F .*= Prop
    U = ifft(F)
    Ur[:, ii+1] = U[:, Int(N/2)+1]
end

# Plotting the propagated field
heatmap(0:dz:z, Y[:, Int(N/2)+1]./w0, abs.(Ur), title="Propagacion del campo en z", xlabel="z/zR", ylabel="y/w₀")

