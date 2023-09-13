%% ********************* PROPAGADOR DE FOURIER ****************************
% Fisica Experimental 1 (EM2022) - Prof. Raul Hernandez 
clear all; close all;
set(0,'defaultTextInterpreter','latex');

%% Datos propios
lambda = 633e-9;            % longitud de onda [m]
k = 2*pi/lambda;            % numero de onda [1/m]
z0 = 1;

%% *************** Definicion de parametros fisicos ***********************
w0 = 0.5e-3;                % Cintura del haz Gaussiano [m]
zR = k*w0^2/2;              % Distancia de Rayleigh [m]
z = 1*zR;                 % Distancia de propagacion maxima en unidades 
                            % de la distancia de Rayleigh [m]


%% *************** Definicion de parametros numericos *********************
N = 2^9;                % Numero de puntos, escogemos una potencia de 2
                        % para que el calculo de la transformada sea mas
                        % rapido

% Creamos el vector de indices
NV = (-N/2:1:N/2-1)*2/N;    % Con esto nos aseguramos de tener la frecuencia 0
                            % y de que NV tiene N elementos
                       
L = 4*w0;                   % Ventana numerica, las unidades son metros [m]
dx = 2*L/N;                 % Intervalo de muestreo
kmax = pi/dx;               % Frecuencia maxima
[X,Y] = meshgrid(NV*L);     % Creamos la malla de evaluacion en el espacio real
[Kx,Ky] = meshgrid(NV*kmax);    %Creamos la malla en el espacio de las frecuencias



r = sqrt(X.^2+Y.^2);        % Coordenada radial
kt = sqrt(Kx.^2+Ky.^2);     % Coordenada radial en el espacio de frecuencias

nz = 300;                   % Numero de puntos en z
dz = z/nz;                  % Tamaño de paso en z, depende de como se quiera propagar
                            % se puede ir evaluando el propagador en cada
                            % coordenada z, o ir propagando en intervalos
                            % de tamaño dz mientras se va actualizando el
                            % campo.

                            
%% ************************ Perfil Inicial *******************************
flens = 0.4*zR;             % Distancia focal de la lente en unidades de zR
Tlens = 1;%exp(-1i*k/(2*flens)*r.^2);             % Funcion de transmitancia de lente
k1 = 10000;                 % Frecuencia espacial transversal para el haz Bessel
f = exp(-r.^2/w0^2).*besselj(0,k1*r).*Tlens;  % Haz Bessel-Gauss
%f = exp(-r.^2/w0^2).*Tlens;                    % Haz Gaussiano

U0 = f;                     % Asignamos el campo inicial
Ur = zeros(N,nz+1);         % Definimos la matriz para guardar la propagacion
Ur(:,1) = U0(:,N/2+1);      % Asignamos el perfil inicial al primer vector


%% ******************  Grafica del campo inicial *************************
figure(1),subplot(1,2,1),surf(X(N/2+1,:)/w0,Y(:,N/2+1)/w0,abs(U0))
shading interp,lighting phong, view(2)
ejes = gca;
ejes.FontSize = 13;
title('Campo $U_{0}$')
xlabel('$x/w_{0}$','FontSize',20);
ylabel('$y/w_{0}$','FontSize',20);
xlim([-L L]/w0),ylim([-L L]/w0),axis square;


%% ****************************** Propagador ******************************
Prop = (exp(-1i*0.5*dz*(kt.^2)/k));      %Propagador paraxial
% Prop = exp(1i*dz*sqrt(k^2-kt.^2));       %Propagador no-paraxial

F = fftshift(fft2(U0));     % Calculo del espectro angular
for ii=1:nz
    F = F.*Prop;            % Propagacion en intervalos de tamaño dz
    U = ifft2(F);           % Recuperamos el campo en el plano real
    Ur(:,ii+1) = U(:,N/2+1);
end


%% ****************** Esta forma de propagar es mas lenta ***************** 
% zvec = linspace(0,z,nz);
% 
% A = fftshift(fft2(U0));     % Calculo del espectro angular
% for ii=1:length(zvec)
%     Prop = exp(-1i*0.5*zvec(ii)*(kt.^2)/k);
%     F = A.*Prop;            % Propagacion en intervalos de tamaño zvec(ii)
%     U = ifft2(F);           % Recuperamos el campo en el plano real
%     Ur(:,ii+1) = U(:,N/2+1);
% end


%% *******************  Grafica del campo propagado ***********************
figure(1),subplot(1,2,2),surf((0:dz:z)/zR,Y(:,N/2+1)/w0,abs(Ur));
shading interp,lighting phong, view(2);   
ejes = gca;
ejes.FontSize = 13;
title('Propagacion del campo en $z$')
xlabel('$z/z_{R}$','FontSize',20);
ylabel('$y/w_{0}$','FontSize',20);
xlim([0 z]/zR),ylim([-L L]/w0);
axis square;

