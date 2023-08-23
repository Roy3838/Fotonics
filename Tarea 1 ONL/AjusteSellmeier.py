import numpy as np
import matplotlib.pyplot as plt

# Datos
lmbda = np.array([0.330, 0.532, 0.6328, 1.064, 1.57, 2.804])  # Longitud de onde
n = np.array([2.47379, 2.23421, 2.20271, 2.15601, 2.1375, 2.1029])  # indice de refraccion

A_i = np.array([2.5**2, 2.5, 1, 1, 0])
sellmeier_i = A_i[0] + lmbda**2 * (A_i[1] / (lmbda**2 - A_i[2])) + lmbda**2 * (A_i[3] / (lmbda**2 - A_i[4]))
sellmeier0 = sellmeier_i
n2 = n**2

# Recocido simulado
A_best = A_i.copy()
ti = 3
tf = 1e-4

while ti > tf:
    print(ti)
    r = np.random.uniform(-2 * np.exp(-5 + ti), 2 * np.exp(-5 + ti), 5)
    A_ti = A_best.copy()
    for jj in range(10001):
        r_new = np.random.uniform(-2 * np.exp(-5 + ti), 2 * np.exp(-5 + ti), 5)
        A_t = A_ti + -2 * (r_new < -A_ti) * r_new + r_new
        sellmeier = A_t[0] + lmbda**2 * (A_t[1] / (lmbda**2 - A_t[2])) + lmbda**2 * (A_t[3] / (lmbda**2 - A_t[4]))
        if np.linalg.norm(sellmeier - n2) < np.linalg.norm(sellmeier0 - n2):
            sellmeier0 = sellmeier
            A_best = A_t
            A_ti = A_t
        elif np.random.rand() > np.exp(-ti):
            A_ti = A_t
    ti = ti * 0.98

# Visualize solution
lmbda_var = np.linspace(0, 3, 10000)
sellmeier = A_best[0] + lmbda_var**2 * (A_best[1] / (lmbda_var**2 - A_best[2])) + lmbda_var**2 * (A_best[3] / (lmbda_var**2 - A_best[4]))

plt.plot(lmbda_var, np.real(np.sqrt(sellmeier)), linewidth=3, color='r')
plt.scatter(lmbda, n, edgecolor='b', facecolor='b')
plt.xlim([0, 3])
plt.ylim([2, 2.6])
plt.title('Ajuste de curva de Sellmeier')
plt.legend(['Curva de Sellmeier', 'Datos Experimentales'])
plt.xlabel('Longitud de onda (µm)')
plt.ylabel('Indice de refraccion (n)')
plt.show()

print(A_best)