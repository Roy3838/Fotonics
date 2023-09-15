

# # gamma = (1im* k)/2 * ((q1*A + B)/q1*B)
# # alpha = 1im* k * (2*B) * x2
# # beta  = 1im* k * (2*B) * y2

# # U = k*exp(1im*k*L_0)/(12*pi) * exp((i*k)/(2*B) * (D * r^2)) * (pi)(gamma) * exp(alpha^2 / gamma + beta^2 / gamma)


# using Plots

# function meshgrid(x, y)
#     X = repeat(x', length(y), 1)
#     Y = repeat(y, 1, length(x))
#     return X, Y
# end


# using LinearAlgebra: I

# function U_2(r_2, k, L_0, D, B, q_1, A, x_2, y_2)
#     # Definimos las constantes i y π
#     i = I
#     π = pi
    
#     # Calculamos γ, α, y β
#     γ = i * k / 2 * (q_1 * A + B) / (q_1 * B)
#     α = i * k * 2B * x_2
#     β = i * k * 2B * y_2
    
#     # Finalmente, calculamos y retornamos U_2
#     return (k * exp(i * k * L_0) / (12π)) * exp(i * k * r_2^2 * D / (2B)) * π * γ * exp((α^2 + β^2) / γ)
# end



# # Original values from the code
# lambda = 633e-9
# k_0 = 2*pi/lambda
# N = 512
# L = 20 * 0.5e-3
# dx = L/N
# NV = collect(-N/2:N/2-1)
# xs = NV*dx
# ys = NV*dx
# Xs, Ys = meshgrid(xs, ys)

# # Placeholder ABCD matrix elements for a thin lens (Replace with actual values)
# f = 0.05 # Focal length of the lens
# A = 1
# B = 1
# C = 1
# D = 1
# q1 = 1.0 # Placeholder value, replace with actual value if available
# L_0 = 3 # Placeholder value, replace with actual value if available

# Rz = Z[zz] * (1 + (Zr / Z[zz])^2)  # Cálculo de la parte real
# Wz = w_0^2 * (1 + (Z[zz] / Zr)^2)  # Parte imaginaria
# qinv = 1 / Rz + (2im / (k_0 * Wz))
# q = 1 / qinv

# gamma = (1im* k)/2 * ((q1*A + B)/q1*B)
# alpha = 1im* k * (2*B) * x2
# beta  = 1im* k * (2*B) * y2

# U = k*exp(1im*k*L_0)/(12*pi) * exp((i*k)/(2*B) * (D * r^2)) * (pi)(gamma) * exp(alpha^2 / gamma + beta^2 / gamma)


# # any(isnan, U) || any(isinf, U)

# # Plot the magnitude of U
# heatmap(xs, ys, abs.(U), color=:viridis, xlabel="x2", ylabel="y2", title="Magnitude of U")




using Plots
# using LinearAlgebra: I

# Definición de la función U_2
function U_2(r_2, k, L_0, D, B, q_1, A, x_2, y_2)
    i = I
    π = pi
    
    γ = i * k / 2 * (q_1 * A + B) / (q_1 * B)
    α = i * k * 2B * x_2
    β = i * k * 2B * y_2
    
    return (k * exp(i * k * L_0) / (12π)) * exp(i * k * r_2^2 * D / (2B)) * π * γ * exp((α^2 + β^2) / γ)
end

# Valores originales del código
lambda = 633e-9
k_0 = 2 * pi / lambda
N = 512
L = 20 * 0.5e-3
dx = L / N
NV = collect(-N/2:N/2-1)
xs = NV * dx
ys = NV * dx
Xs = [x for x in xs, y in ys]
Ys = [y for x in xs, y in ys]

# Elementos de la matriz ABCD para una lente delgada
f = 0.05
A = 1
B = 1
C = 1
D = 1
q1 = 1.0
L_0 = 3

# Calculando la función U_2 en todo el plano y
U_2_values = [U_2(sqrt(x^2 + y^2), k_0, L_0, D, B, q1, A, x, y) for x in xs, y in ys]

# Graficando
# heatmap(xs, ys, abs.(U_2_values), color=:viridis, xlabel="x (m)", ylabel="y (m)", title="Magnitude of U_2 in y-plane")




