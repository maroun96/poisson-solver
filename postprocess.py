import h5py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['savefig.dpi'] = 300

with h5py.File("u.hdf5", "r") as f:
    u_sol = f["u_sol"][:]

# print(u_sol.shape)



x = np.linspace(-1, 1, 129)
y = np.linspace(-1, 1, 129)
z = np.linspace(-1, 1, 129)
X, Y, Z= np.meshgrid(x, y, z)

# u_true = X*X + Y*Y + Z*Z
u_true = np.cos(2*np.pi*X)*np.cos(2*np.pi*Y)*np.cos(2*np.pi*Z)

print(Z[:,:,64])


plt.contour(X[:,:,64], Y[:,:,64], u_sol[:,:,64], colors="blue",linestyles='solid', levels=10)
plt.contour(X[:,:,64], Y[:,:,64], u_true[:,:,64], colors="red", alpha=0.6, linestyles='dashed',levels=10)
plt.xlabel("x")
plt.ylabel("y")
# Creating custom legend entries
legend_elements = [
    Line2D([0], [0], color='blue', lw=2, label='Numerical Solution'),
    Line2D([0], [0], color='red', lw=2, label='Exact Solution')
]

# Adding the legend to the plot
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2, fontsize='small')
plt.tight_layout()
plt.savefig("comparison_plot.png")



# err_norms = np.loadtxt("error_inf.dat")
# n_grid = np.array([3,5,9,17,33,65,129])
# dx = 2/(n_grid-1)
# plt.plot(dx, err_norms, linestyle="--", marker="x", 
#          markerfacecolor='red', markeredgecolor='red')
# plt.grid(True)
# plt.xlabel("dx")
# plt.ylabel("$|u-u_{exact}|_{\infty}$")
# plt.suptitle("Error convergence")
# plt.savefig("error_convergence.png")

# with h5py.File("poisson_mat.hdf5", "r") as f:
#     print(list(f.keys()))

# A = np.fromfile("poisson_mat.dat")
# print(A)



