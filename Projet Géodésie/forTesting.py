import matplotlib.pyplot as plt
import numpy as np
from math import *
from New_function import *

# Paramètres de l'ellipse
a = 3
b = 2

# Angle φ en radians
Phi = -pi / 3

# Calculer les angles β et ψ
Beta = beta_from_phi(a, b, Phi)
Psi = psi_from_phi(a, b, Phi)

# Création des points de l'ellipse
theta = np.linspace(-np.pi / 2, np.pi / 2, 10000)
x_ellipse = a * np.cos(theta)
y_ellipse = b * np.sin(theta)

# Création des points pour l'angle φ
phi_x = [0, a * np.cos(Phi)]
phi_y = [0, b * np.sin(Phi)]

# Création des points pour l'angle β
beta_x = [0, a * np.cos(Beta)]
beta_y = [0, b * np.sin(Beta)]

# Création des points pour l'angle ψ
psi_x = [0, a * np.cos(Psi)]
psi_y = [0, b * np.sin(Psi)]

# Dessiner l'ellipse
plt.plot(x_ellipse, y_ellipse)

# Dessiner le cercle (en prenant a comme rayon)
plt.plot(x_ellipse, a * np.sin(theta))

# Relier le centre de l'ellipse à l'angle β
plt.plot([0, a * np.cos(Beta)], [0, b * np.sin(Beta)], color="g", linestyle="--")

# Relier le centre de l'ellipse à l'angle ψ
plt.plot([0, a * np.cos(Psi)], [0, b * np.sin(Psi)], color="r", linestyle="--")

# Dessiner les angles φ, β, ψ
plt.plot(phi_x, phi_y, color="b", label=f"Phi= {round((np.degrees(Phi)),0)}°")
plt.plot(beta_x, beta_y, color="g", label=f"Beta= {round(np.degrees(Beta),0)}°")
plt.plot(psi_x, psi_y, color="r", label=f"Psi= {round(np.degrees(Psi),0)}°")

# Configurer le graphique
plt.axhline(0, color="black", linewidth=0.2)
plt.axvline(0, color="black", linewidth=0.2)
plt.gca().set_aspect("equal", adjustable="box")


plt.xticks([])
plt.yticks([])
# Légende
plt.legend()
plt.legend(loc="lower right")

# Afficher le graphique
plt.show()
