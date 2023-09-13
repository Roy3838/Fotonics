

# gamma = (1im* k)/2 * ((q1*A + B)/q1*B)
# alpha = 1im* k * (2*B) * x2
# beta  = 1im* k * (2*B) * y2

# U = k*exp(1im*k*L_0)/(12*pi) * exp((i*k)/(2*B) * (D * r^2)) * (pi)(gamma) * exp(alpha^2 / gamma + beta^2 / gamma)


using Plots

function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

# Original values from the code
lambda = 633e-9
k_0 = 2*pi/lambda
N = 512
L = 20 * 0.5e-3
dx = L/N
NV = collect(-N/2:N/2-1)
xs = NV*dx
ys = NV*dx
Xs, Ys = meshgrid(xs, ys)

# Placeholder ABCD matrix elements for a thin lens (Replace with actual values)
f = 0.05 # Focal length of the lens
A = 1
B = 1
C = 1
D = 1
q1 = 1.0 # Placeholder value, replace with actual value if available
L_0 = 3 # Placeholder value, replace with actual value if available


# any(isnan, U) || any(isinf, U)

# Plot the magnitude of U
heatmap(xs, ys, abs.(U), color=:viridis, xlabel="x2", ylabel="y2", title="Magnitude of U")








