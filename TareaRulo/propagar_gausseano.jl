# Assuming propagar.jl is in the same directory as this script
include("propagar.jl")

using Plots
using SpecialFunctions

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

# ******************** Perfil Inicial ********************
flens = 0.4 * zR
Tlens = 1#exp.(-1im * k / (2 * flens) * r.^2)
k1 = 10000
f = exp.(-r.^2/w0^2) .* Tlens
#f = exp.(-r.^2/w0^2) .* besselj0.(k1 * r) .* Tlens

U0 = f
nz = 300
dz = z / nz

# Call the propagar function to get propagated fields
propagated_fields = propagar(U0, z, dz, Kx, Ky, k)

Ur = zeros(Complex{Float64}, N, nz+1)

# Extract the relevant slices
Ur = propagated_fields[:, Int(N/2)+1, :]


# Plot the propagated fields
heatmap(0:dz:z, Y[:, Int(N/2)+1]./w0, abs.(Ur), title="Propagacion del campo en z", xlabel="z/zR", ylabel="y/w₀", color=:viridis)

# Save as a PNG file
savefig("propagacionSinLente.png")

