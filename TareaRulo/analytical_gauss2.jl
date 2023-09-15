

include("propagar.jl")

# Propagador de Fourier

# Parametros iniciales
N = 2^8  # Numero de puntos a evaluar
w_0 = 3e-4  # Diametro de un haz en m
Lambda = 632.8e-9  # Longitud de onda
k0 = 2 * π / Lambda  # Número de onda
Lx = 12 * w_0  # Límite espacial en terminos
Ly = 12 * w_0
Zr = k0 * w_0^2 / 2  # Distancia de Rayleigh

# Vector de ventana
x = Lx * (2/N) * (-N/2:N/2-1)
y = Ly * (2/N) * (-N/2:N/2-1)
X, Y = meshgrid(x, y)  # Malla de trabajo
#phi, r = cart2pol.(X, Y)  # Se define en coordenadas polares


Mz = 1.5 * Zr       # Distancia maxima de propagación
dz = Mz / N         # Delta de propagación en z
nz = Int(Mz / dz)   # Número de puntos en z
Z = 0:dz:Mz         # Vector de z

# Límites espectrales
dx = 2 * Lx / N
dy = 2 * Ly / N
kx = (2/N) * (-N/2:N/2-1) * (π / dx)  # Vectores auxiliares de frecuencia
ky = (2/N) * (-N/2:N/2-1) * (π / dy)
Kx, Ky = meshgrid(kx, ky)             # Ventana espectral
Kt = sqrt.(Kx.^2 + Ky.^2)             # K transversal
Kz = sqrt.(k0^2 .- abs.(Kt).^2)

u0 = exp.((-X.^2 .- Y.^2)./w_0.^2)  # Campo inicial

kap1 = 8665
m = [0 1]

ii = 2
U0_propagated = zeros(Complex{Float64}, N, N, nz)  # Inicialización de U
U0_propagated[:, :, 1] = u0
A = 1
B = Z
C = 0
D = 1

using SpecialFunctions: besselj0, besselj1, besselj
using Plots

# for del cálculo analítico
for zz in 2:nz
    Rz = Z[zz] * (1 + (Zr / Z[zz])^2)  # Cálculo de la parte real
    Wz = w_0^2 * (1 + (Z[zz] / Zr)^2)  # Parte imaginaria
    qinv = 1 / Rz + (2im / (k0 * Wz))
    q = 1 / qinv
    q2 = (q * A + B[zz]) / (q * C + D)  # Cálculo de q2
    kap2 = kap1 / (A + B[zz] / q)

    U0_propagated[:, :, zz] = exp.((kap1 * kap2 * B[zz]) / (2im * k0)) * exp.(1im * k0 * Z[zz]) / (A + B[zz] / q) * exp.(1im * k0 * (X.^2 + Y.^2) / (2 * Rz)) .* besselj0.(kap2 * sqrt.(X.^2 + Y.^2)) .* exp.(1im * m[ii] * atan.(Y, X))  # Campo con parametro q1

end


# Corte para graficar en las coordenadas 'y'
Ua0_y = U0_propagated[Int(end/2), :, :]

# Asumiendo que la primera dimensión de Ua0_y corresponde a ys
ys_corrected = y[1:size(Ua0_y, 1)]

# La segunda dimensión de Ua0_y todavía corresponde a Z
Z_corrected = Z[2:2+size(Ua0_y, 2)-1]


# Graficando el corte Ua0_y
heatmap(ys_corrected, Z_corrected, abs.(Ua0_y), color=:viridis, xlabel="y (m)", ylabel="Z (m)", title="Magnitude of Ua0 in y-Z plane")


