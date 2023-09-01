using Plots
using LinearAlgebra

using LaTeXStrings
using Gadfly, Compose

n1 = 1.471
n2 = 1.46
n = n2/n1
θc = asin(n)
θ = θc:0.001:π/2

#ΦTE y ΦTM
ΦTEr =  -2*atan.(sqrt.(round.((sin.(θ).^2 .- n.^2), digits = 10)),cos.(θ))
ΦTEt = (1/2).*ΦTEr
plot_n_2_3 = Plots.plot(θ, ΦTEr, linewidth = 3, color = :blue, label = "ΦTEr")
plot!(θ, ΦTEt, linewidth = 3, color = :red, label = "ΦTEt")
xlabel!("Ángulo de Incidencia "*L"(rad)")
ylabel!("Desfase "*L"(rad)")
title!("Desfase para "*L"n_{1} = 1.471, n_{2} = 1.46")
savefig(plot_n_2_3, "Phase_Changes_P2_r_t")

#Modos de propagacion
λ0 = 632.8*10^-9
k1 = n1*2*π/λ0
a = 5*10^-6

α = 0:0.0001:(pi/2-θc)
f1 = 2*k1.*a.*sin.(α)
f2 = -2*atan.(sqrt.((n1^2).*cos.(α).^2 .- n2^2)./(n1.*sin.(α)))
fig = Plots.plot(α, f1, linewidth = 3, color = :blue, label = L"k_{1}d\sin(α)", xlims = (0, 0.14))
plot!(α, f2, linewidth = 3, color = :red, label = L"-2arctan[-]")
plot!(α, f2+f1, linewidth = 3, color = :black, label = L"k_{1}d\sin(\alpha)-2arctan[-]")
for ii in 0:1:6
    display(plot!([0, α[end]], [ii*π, ii*π], label = false, color = :black, linestyle = :dash))
    display(annotate!((α[end]+0.01, ii*π, "m = "*string(ii)), fontsize = 5))
end
title!("Modos de propagación")
xlabel!("Ángulo (rad)")
ylabel!("Desfase (rad)")
savefig(fig, "Modos_de_propagacion")

#Angulos de modos de propagación

for ii in 0:5
    modo = argmin(abs.(f2+f1.-ii*pi))
    angulo_modal = 90-360*α[modo]/(2*π)
    println(angulo_modal)
end

#Monomodo (misma fibra) variando λ0
a = 5*10^-6
λ0 = 4*a√(n1^2-n2^2)
k1 = n1*2*π/λ0


α = 0:0.0001:(pi/2-θc)
f1 = 2*k1.*a.*sin.(α)
f2 = -2*atan.(sqrt.((n1^2).*cos.(α).^2 .- n2^2)./(n1.*sin.(α)))
fig = Plots.plot(α, f1, linewidth = 3, color = :blue, label = L"k_{1}d\sin(α)", xlims = (0, 0.14))
plot!(α, f2, linewidth = 3, color = :red, label = L"-2arctan[-]")
plot!(α, f2+f1, linewidth = 3, color = :black, label = L"k_{1}d\sin(\alpha)-2arctan[-]")
for ii in 0:1:1
    display(plot!([0, α[end]], [ii*π, ii*π], label = false, color = :black, linestyle = :dash))
    display(annotate!((α[end]+0.01, ii*π, "m = "*string(ii)), fontsize = 5))
end
title!("Longitud de corte λ = "*string(round(λ0, digits = 9)))
xlabel!("Ángulo (rad)")
ylabel!("Desfase (rad)")
savefig("Longitud_de_onda_de_corte")


#Monomodo (misma fibra) variando la distancia
a = (7.048/2)*10^-6
λ0 = 632.8*10^-9
k1 = n1*2*π/λ0

α = 0:0.0001:(pi/2-θc)
f1 = 2*k1.*a.*sin.(α)
f2 = -2*atan.(sqrt.((n1^2).*cos.(α).^2 .- n2^2)./(n1.*sin.(α)))
fig = Plots.plot(α, f1, linewidth = 3, color = :blue, label = L"k_{1}d\sin(α)", xlims = (0, 0.14))
plot!(α, f2, linewidth = 3, color = :red, label = L"-2arctan[-]")
plot!(α, f2+f1, linewidth = 3, color = :black, label = L"k_{1}d\sin(\alpha)-2arctan[-]")
for ii in 0:1:5
    display(plot!([0, α[end]], [ii*π, ii*π], label = false, color = :black, linestyle = :dash))
    display(annotate!((α[end]+0.01, ii*π, "m = "*string(ii)), fontsize = 5))
end
title!("Longitud de corte λ = "*string(round(λ0, digits = 9)))
xlabel!("Ángulo (rad)")
ylabel!("Desfase (rad)")
savefig("5Modos_var_distancia.png")

