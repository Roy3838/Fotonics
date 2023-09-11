using Plots
include("propagar.jl")

N = 512
U_0 = zeros(Complex{Float64},N,N)
lambda = 633e-9
w0 = 0.5e-3
L = 20*w0
dx = L/N
NV = collect(-N/2:N/2-1)
xs = NV*dx
ys = NV*dx

Xs, Ys = meshgrid(xs, ys)
r2 = Xs.^2 + Ys.^2

U_0 = exp.(-r2/w0^2) # Gaussiano inicial
# U_0 es un cuadrado
# for i in 1:N
#     for j in 1:N
#         if abs(xs[i]) < 0.5e-3 && abs(ys[j]) < 0.5e-3
#             U_0[i,j] = 1
#         end
#     end
# end


k_0 = 2*pi/lambda
k_max = pi/dx
k_x = k_max * (2/N) * NV
k_y = k_max * (2/N) * NV
KXs, KYs = meshgrid(k_x, k_y)

zR = pi*w0^2/lambda

# # Aplicar lente delgada
f = 0.05 # Distancia focal de 5 cm
phase = exp.(-1im *( k_0 ./ (2 .* f)) * (Xs.^2 + Ys.^2))
U_0 = U_0 .* phase

z_final = 0.7e-2*zR

nz = 100 # número de pasos en z
dz = z_final/nz
U_z = zeros(Complex{Float64}, N, N, nz) # Matriz tridimensional

U_z[:,:,1] = U_0

# Propagar y almacenar en U_z
for i in 1:nz-1
    U_z[:, :, i+1] = propagar(U_z[:,:,i], dz, KXs, KYs, k_0)
end

#heatmap(abs.(U_z[:,:,1]).^2, color=:viridis)

z_values = range(0, stop=z_final, length=nz)
intensity_plane = abs.(U_z[Int(N/2), :, :]).^2
heatmap(ys, z_values, intensity_plane', color=:viridis, xlabel="y", ylabel="z", title="Propagation in z at x=0")


#heatmap(abs.(U_z[:, :, Int(nz-1)]).^2, color=:viridis)



# U = propagar(U_0, z, KXs, KYs, k_0)
# 





