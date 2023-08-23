import numpy as np
import matplotlib.pyplot as plt

def Xchi_w (w, w0, gamma,m,N,q):
    return N*q**2/(m*(w0**2-w**2-1j*w*gamma))


def n_bar(Xchi):
    return np.sqrt(1+Xchi)

# Parámetros
omega_0 = 1  # Dado en el problema
omega_p = 1  # Dado en el problema
gamma_values = [0.06, 0.2, 0.5]  # Dados en el problema
omega = np.linspace(0.1, 3, 500)  # Rango de frecuencias

# Gráficos
plt.figure(figsize=(12, 6))

# Gráfico de n
plt.subplot(1, 2, 1)
for gamma in gamma_values:
    plt.plot(omega, np.real(n_bar(Xchi_w(omega,omega_0,gamma,1,1,1))), label=f"γ = {gamma}")
plt.xlabel("ω")
plt.ylabel("Índice de refracción (n)")
plt.title("Índice de refracción vs ω")
plt.legend()
plt.grid(True)

# Gráfico de k
plt.subplot(1, 2, 2)
for gamma in gamma_values:
    plt.plot(omega, np.imag(n_bar(Xchi_w(omega,omega_0,gamma,1,1,1))), label=f"γ = {gamma}")
plt.xlabel("ω")
plt.ylabel("Coeficiente de extinción (k)")
plt.title("Coeficiente de extinción vs ω")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()


