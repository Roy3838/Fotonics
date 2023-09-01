using Plots
using LinearAlgebra
using LaTeXStrings
using Gadfly, Compose
#Parámetros
N = 1000
n21 = LinRange(0, 1, N)
n21_mat = zeros(Float64, N, N)
for i in 1:N
    n21_mat[i, i:end] = LinRange(n21[i], n21[i], N-i+1)
end
n21_mat
θc = asin.(n21)
θ1 = zeros(Float64, N, N)
θ1[1, 1:end]
θc[1:end]
for i in 1:N
    θ1[i, i:end] = θc[i:end]
end
θ1
#ΦTE Y ΦTM
ΦTE =  -2*atan.(sqrt.(round.((sin.(θ1).^2 .- n21_mat.^2), digits = 10))./cos.(θ1))
plot_ΦTE =  Plots.heatmap(θc, n21, ΦTE, colorbar_title=L"\phi_{TE}")
plot!([0.55, 0.65], [0.7, 0.6], arrow = true, linewidth = 2, color = :black)
annotate!((0.5, 0.75, L"\theta_{C}"))

xlabel!("Ángulo de incidencia (θ)")
ylabel!("n relativo "*L"(n_2/n_1)")
title!("Cambio de fase en RTI")
savefig(plot_ΦTE, "Phase_Change_PhiTE")

((sin.(θ1)./n21_mat).^2 .-1)./n21_mat
pseudo_ΦTM = (sqrt.(round.((sin.(θ1)./n21_mat).^2 .- 1, digits = 10)))./(n21_mat.*cos.(θ1))
pseudo_ΦTM = replace!(pseudo_ΦTM, NaN=>0)
ΦTM =  -2*atan.(pseudo_ΦTM)
plot_ΦTM = Plots.heatmap(θc, n21, ΦTM, colorbar_title=L"\phi_{TM}")

xlabel!("Ángulo de incidencia (θ)")
ylabel!("n relativo "*L"(n_2/n_1)")
title!("Cambio de fase en RTI")
plot!([0.55, 0.65], [0.7, 0.6], arrow = true, linewidth = 2, color = :black)
annotate!((0.5, 0.75, L"\theta_{C}"))
savefig(plot_ΦTM, "Phase_Change_PhiTM")

# n21 = 2/3
n = 2/3
θc = asin(n)
θ = θc:0.001:π/2

#ΦTE y ΦTM
ΦTE =  -2*atan.(sqrt.(round.((sin.(θ).^2 .- n.^2), digits = 10)),cos.(θ))
pseudo_ΦTM = (sqrt.(round.((sin.(θ)./n).^2 .- 1, digits = 10)))
ΦTM =  -2*atan.(pseudo_ΦTM, n.*cos.(θ))
plot_n_2_3 = Plots.plot(θ, ΦTE, linewidth = 3, color = :blue, label = "ΦTE")
plot!(θ, ΦTM, linewidth = 3, color = :red, label = "ΦTM")
xlabel!("Ángulo de Incidencia "*L"(rad)")
ylabel!("Desfase "*L"(rad)")
title!("Desfase para "*L"n_{21}=2/3")
savefig(plot_n_2_3, "Phase_Changes_n_2_3")

#Maximo desfase entre ΦTE y ΦTM para n=0.5
n = 1/2
θc = asin(n)
θ = θc:0.001:π/2
ΦTE =  -2*atan.(sqrt.(round.((sin.(θ).^2 .- n.^2), digits = 10)),cos.(θ))
pseudo_ΦTM = (sqrt.(round.((sin.(θ)./n).^2 .- 1, digits = 10)))
ΦTM =  -2*atan.(pseudo_ΦTM, n*cos.(θ))
plot_n_1_2 = Plots.plot(θ, ΦTE, linewidth = 3, label = L"Φ_{TE}", color = :blue, legend =:right)
plot!(θ, ΦTM, linewidth = 3, label = L"Φ_{TM}", color = :red)
xlabel!("Ángulo de Incidencia "*L"(rad)")
ylabel!("Desfase "*L"(rad)")
title!("Desfase para "*L"n_{21}= 0.5")
plot!(θ, -abs.(ΦTE-ΦTM), linewidth = 3, label = L"ΔΦ = -|Φ_{TE}-Φ_{TM}|", color = :black)
maxi = findmax(abs.(ΦTE-ΦTM))
scatter!([θ[max[2]]], [-max[1]], color = :red, label = "Max ΔΦ")
plot!([θ[max[2]], θ[max[2]]], [-π, 0], color = :black, linestyle = :dash, label = false)
plot!([0.75, θ[max[2]]], [-2.7, -2.7], arrow = true, linewidth = 2, color = :black, label = false)
annotate!((0.87, -2.7, "θ = "*string(round(θ[max[2]], digits = 3))))
savefig(plot_n_1_2, "Desfase_n_1_2")

maxi[1]