import numpy as np
import matplotlib.pyplot as plt

# Definición de funciones
def theta_c(n21):
    return np.arcsin(n21)

def phi_TE(theta1, n21):
    return np.where(theta1 <= theta_c(n21), 0, np.pi)

def phi_TM(theta1, n21):
    condition = theta1 <= theta_c(n21)
    term = 2 * n21 * np.cos(theta1) / np.sqrt(np.sin(theta1)**2 - n21**2)
    return np.where(condition, 0, -np.arctan(term))

# Crear una malla de valores
n21_values = np.linspace(0, 1, 400)
theta1_values = np.linspace(0, np.pi/2, 400)
Theta1, N21 = np.meshgrid(theta1_values, n21_values)

Phi_TE = phi_TE(Theta1, N21)
Phi_TM = phi_TM(Theta1, N21)

# Graficar
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Gráfica de phi_TE
c1 = ax1.contourf(Theta1, N21, Phi_TE, levels=50, cmap='viridis')
ax1.set_title(r'$\phi_{TE}$', fontsize=16)
ax1.set_xlabel(r'$\theta_1$', fontsize=14)
ax1.set_ylabel(r'$n_{21}$', fontsize=14)
fig.colorbar(c1, ax=ax1, orientation='vertical')

# Gráfica de phi_TM
c2 = ax2.contourf(Theta1, N21, Phi_TM, levels=50, cmap='viridis')
ax2.set_title(r'$\phi_{TM}$', fontsize=16)
ax2.set_xlabel(r'$\theta_1$', fontsize=14)
ax2.set_ylabel(r'$n_{21}$', fontsize=14)
fig.colorbar(c2, ax=ax2, orientation='vertical')

plt.tight_layout()
plt.show()



# Definición de n21
n21_value = 2/3

# Calcular los cambios de fase para n21 = 2/3
Phi_TE_single = phi_TE(theta1_values, n21_value)
Phi_TM_single = phi_TM(theta1_values, n21_value)

# Graficar
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Gráfica de phi_TE para n21 = 2/3
ax1.plot(theta1_values, Phi_TE_single, label=f'$n_{{21}} = {n21_value}$', color='blue')
ax1.axvline(x=theta_c(n21_value), linestyle='--', color='red', label=r'$\theta_C$')
ax1.set_title(r'$\phi_{TE}$', fontsize=16)
ax1.set_xlabel(r'$\theta_1$', fontsize=14)
ax1.set_ylabel(r'Cambio de fase', fontsize=14)
ax1.legend()

# Gráfica de phi_TM para n21 = 2/3
ax2.plot(theta1_values, Phi_TM_single, label=f'$n_{{21}} = {n21_value}$', color='blue')
ax2.axvline(x=theta_c(n21_value), linestyle='--', color='red', label=r'$\theta_C$')
ax2.set_title(r'$\phi_{TM}$', fontsize=16)
ax2.set_xlabel(r'$\theta_1$', fontsize=14)
ax2.set_ylabel(r'Cambio de fase', fontsize=14)
ax2.legend()

plt.tight_layout()
plt.show()
