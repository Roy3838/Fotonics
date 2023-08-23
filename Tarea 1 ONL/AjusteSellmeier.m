clear; close all; clc;
%% Datos
lmbda = [0.330 0.532 0.6328 1.064 1.57 2.804]; %Longitud de onde
n = [2.47379 2.23421 2.20271 2.15601 2.1375 2.1029]; %indice de refraccion

A_i = [2.5^2 2.5 1 1 0];
sellmeier_i = A_i(1) + lmbda.^2.*(A_i(2)./(lmbda.^2 - A_i(3)))...
            +lmbda.^2.*(A_i(4)./(lmbda.^2 - A_i(5)));
sellmeier0 = sellmeier_i;
n2 = n.^2;
%% Recocido simulado
A_best = A_i;
ti = 3;
tf = 1e-7;
while ti > tf 
    pdsteps = makedist('Uniform',-2*exp(-5+ti),2*exp(-5+ti));%Tamaño de paso
    A_ti = A_best;
    for jj = 0:10000
        r = random(pdsteps,1,5);
        A_t = A_ti + -2*(r < -A_ti).*r + r;
        sellmeier = A_t(1) + lmbda.^2.*(A_t(2)./(lmbda.^2 - A_t(3)))...
            +lmbda.^2.*(A_t(4)./(lmbda.^2 - A_t(5)));
        if norm(sellmeier - n2) < norm(sellmeier0-n2)
            sellmeier0 = sellmeier;
            A_best = A_t;
            A_ti = A_t;
        elseif rand > exp(-ti)
            A_ti = A_t;
        end
    end
    ti = ti*0.98;
end
%% Visualize solution
lmbda_var = linspace(0, 3, 10000);
sellmeier = A_best(1) + lmbda_var.^2.*(A_best(2)./(lmbda_var.^2 - A_best(3)))...
            +lmbda_var.^2.*(A_best(4)./(lmbda_var.^2 - A_best(5)));
plot(lmbda_var, real(sqrt(sellmeier)), 'LineWidth', 3, 'Color' ,'r')
hold on 
scatter(lmbda, n, 'filled', 'MarkerEdgeColor', 'b', "MarkerFaceColor", "b")
xlim([0, 3])
ylim([2, 2.6])
title('Ajuste de curva de Sellmeier')
legend('Curva de Sellmeier', 'Datos Experimentales')
xlabel('Longitud de onda (\mum)')
ylabel('Indice de refraccion (n)')


