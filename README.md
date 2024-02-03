# Generación de triplete de fotones entrelazados utilizando TOSPDC en una guía de onda.

La generación de parejas de fotones es crucial en mecánica cuántica, usándose tanto en computadoras cuánticas como en experimentos para entender la teoría cuántica. 
Estas parejas se crean mediante varios métodos en cristales no lineales, como en PDC (donde un fotón se divide en dos) o en SFWM (donde dos fotones se convierten en otros dos). 

Sin embargo, generar tripletes de fotones entrelazados y multipletes de mayor orden implica desafíos tecnológicos significativos. 
Nos enfocamos en la conversión paramétrica descendente espontánea de tercer orden (TOSPDC) en fibras ópticas de sílice fusionado, donde un solo fotón de bombeo se aniquila para generar un triplete de fotones. 
Este proceso es prometedor para la emisión anunciada de pares de fotones y para la generación directa de estados entrelazados de polarización de Greenberger-Horne-Zeilinger (GHZ), sin necesidad de postselección. 

Este documento explorará la teoría detrás de nuestra propuesta de fuentes de tripletes de fotones TOSPDC, enfocándose en el estado del triplete de fotones y las características de coincidencia de fases de TOSPDC en fibras ópticas delgadas. 

Este documento tiene como objetivo proporcionar una visión general de la generación de pares de fotones en fibras ópticas, incluyendo una descripción del estado de dos fotones, y la optimización de las dimensiones de la guía de onda para TOSPDC.

## Efectos No-Lineales
Existen dos métodos de generacion de parejas de fotones mediante cristales no-lineales, uno que considera el efecto $\chi^{(2)}$ y otro que considera el efecto $\chi^{(3)}$.
El efecto de $\chi^{(2)}$ se llama Spontaneus Parametric Down Conversion, al ser de segundo orden de no-linealidad es el mas estudiado. 
También se ha propuesto recientemente el efecto $\chi^{(3)}$ de Spontaneous Four Wave Mixing, en el cual también se generan dos fotones pero a costa de otros dos fotones de bombeo. 
El fenómeno que estudiaremos es el tercer fenómeno $\chi^{(3)}$ llamado Third Order Spontaneous Parametric Down Conversion, que será referenciado como (TOSPDC). 

## Derivación del estádo cuántico
La generación de fotones mediante TOSPDC sigue el siguiente hamiltoniano

$\hat{H}= \frac{3}{4}\epsilon_0 \chi^{(3)} \int dV \hat{E}_p^{(+)} \hat{E}_r^{(-)} \hat{E}_s^{(-)} \hat{E}_i^{(-)}$

Donde $\hat{E}^{(+)}$ y $\hat{E}^{(-)}$ son las partes de frecuencia positiva y precuencia negativa del operador de campo eléctrico.

La generación de fotones en TOSPDC sigue una formulación donde el campo eléctrico para los modos $r,s,i$ se expresa como una función que incluye el operador de aniquilación dependiente del número de onda asociado al modo de propagación en la fibra.

$\hat{E}^{(+)}(r,t) = iA(x,y)\sqrt{\delta k} \sum_k l(w) \exp{[i(kz-wt)]} \hat{a}(k)$

Y se representa los fotones de bombeo como una onda clásica, que al sustituir (2) en el hamiltoniano (1) permite derivar el estado producido por TOSPDC, que se puede escribir en términos del componente de tres fotones del estado $\ket{\Psi_3}$.

$\ket{\Psi} = \ket{0}_r\ket{0}_s\ket{0}_i + \xi \ket{\Psi_3}$

Donde $\xi$ es la eficiencia de conversión.

$\xi=$
$\frac{3\epsilon_0\chi^{(3)} (2\pi)A_0(\delta k)^{3/2} L}{4\hbar}\times \int dx \int dy A_p(x,y)$

$A_r^*(x,y)$ 

$A_s^*(x,y)$

$A_i^*(x,y)$

Y el estado $\ket{\Psi_3}$ está descrito por.

$\ket{\Psi_3} = \sum_{k_r}\sum_{k_s}\sum_{k_i} G_k (k_r,k_s,k_i) \times \hat{a}^\dag(k_r)\hat{a}^\dag(k_s)\hat{a}^\dag(k_i)\ket{0}_r\ket{0}_s\ket{0}_i$


La función $G_k(k_r,k_s,k_i)$ es la amplitud conjunta de número de onda.
Las propiedades espectrales del triplete de fotones están determinadas por la función de amplitud espectral conjunta $G_k(k_r,k_s,k_i)$, que se simplifica a ser $F(w_r,w_s,w_i)$

$F(w_r,w_s,w_i) = \alpha(w_r+w_s+w_i)\cdot\Phi(w_r,w_s,w_i)$

que se relaciona con la Pump Spectral Amplitude $\alpha$ y la función de Phase Matching $\Phi$.

La función de Phase Matching se define como

$\Phi = sinc[\frac{L}{2}\Delta k]\exp{[\frac{iL}{2}\Delta k]}$

donde 

$\Delta k (w_r,w_s,w_i) = k_p(w_r+w_s+w_i) - k_r(w_r) - k_s(w_s) -k_i(w_i) + NL$

Se puede observar el mecanismo de TOSPDC en la anhiquilación de un fotón de bombeo para la creación de los fotones $k_r, k_s, k_i$.
Finalmente, se aborda la contribución no lineal en la desviación de fase, que para fines de ésta investiagación se considera que son 0.

Sabiendo la función de Phase Matching de los fotones generados, se puede saber la función de amplitud espectral conjunta, que te da información del estado de los tres fotones $\ket{\Psi_3}$.





