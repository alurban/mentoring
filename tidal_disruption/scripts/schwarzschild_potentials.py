# Imports.
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib.patheffects as PE
from matplotlib import ticker

# Define the effective potential.
def potential(r, h):
    return -1/r + (h/r)**2/2 - h**2/r**3

# Set array of orbital separation values.
r = np.linspace(1e-4, 100, 1000)

# Construct a set of arrays corresponding to effective potentials
# with different angular momenta.
h = np.array([3.45, 3.85, 4.3][::-1])
Phi = np.array([potential(r, x) for x in h])

# Construct a figure.
fig = plt.figure( figsize=(6, 7.5) )
ax = [[] for x in range(3)]

# Plot the effective potential as a function of orbital separation.
ax[0] = fig.add_subplot(3, 1, 1)
ax[0].plot([0, 100], [0, 0], 'k--', linewidth=0.5)
ax[0].plot(r, Phi[0], 'k-.', linewidth=2., label='$h = %s \, GM/c$' % h[0])
rC = 0.5 * h[0] * (h[0] + np.sqrt(h[0]**2 - 12))
ax[0].annotate('circular', xy=(rC, -0.034), xycoords='data', size=10,
    ha="center", va="top", color='DodgerBlue')
ax[0].scatter([rC], [potential(rC, h[0])], 25, c='DodgerBlue')
for energy in [-0.02, 0.02]:
    condition = np.logical_and(Phi[0] < energy, r > 5)
    ax[0].plot(r[condition], [energy for x in r[condition]], 'DodgerBlue',
        linestyle='dashed', linewidth=1.5)
ax[0].annotate('elliptical', xy=(20, -0.015), xycoords='data', size=10,
    ha="center", va="center", color='DodgerBlue')
ax[0].annotate('hyperbolic', xy=(50, 0.025), xycoords='data', size=10,
    ha="center", va="center", color='DodgerBlue')
ax[0].set_xlim([0, 100])
ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax[0].set_ylim([-0.05, 0.05])
ax[0].set_ylabel('$\Phi$ $(\mu c^2)$')
ax[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
leg = ax[0].legend(loc=4, fontsize=11, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Plot the effective potential as a function of orbital separation.
ax[1] = fig.add_subplot(3, 1, 2)
ax[1].plot([0, 100], [0, 0], 'k--', linewidth=0.5)
ax[1].plot(r, Phi[1], 'k-.', linewidth=2., label='$h = %s \, GM/c$' % h[1])
rC = 0.5 * h[1] * (h[1] + np.sqrt(h[1]**2 - 12))
ax[1].annotate('circular', xy=(rC, 1.05*potential(rC, h[1])), xycoords='data',
    size=10, ha="center", va="top", color='SandyBrown')
ax[1].scatter([rC], [potential(rC, h[1])], 25, c='SandyBrown')
condition = np.logical_and(Phi[1] < -0.025, r > 5)
ax[1].plot(r[condition], [-0.025 for x in r[condition]], 'SandyBrown',
    linestyle='dashed', linewidth=1.5)
ax[1].annotate('elliptical', xy=(15, -0.021), xycoords='data', size=10,
    ha="center", va="center", color='SandyBrown')
condition = np.logical_and(Phi[1] < -0.01, r > 0)
ax[1].plot(r[condition], [-0.01 for x in r[condition]], 'SandyBrown',
    linestyle='dashed', linewidth=1.5)
ax[1].annotate('capture', xy=(30, -0.006), xycoords='data', size=10,
    ha="center", va="center", color='SandyBrown')
ax[1].set_xlim([0, 100])
ax[1].xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax[1].set_ylim([-0.06, 0.01])
ax[1].set_ylabel('$\Phi$ $(\mu c^2)$')
ax[1].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
leg = ax[1].legend(loc=4, fontsize=11, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Plot the effective potential as a function of orbital separation.
ax[2] = fig.add_subplot(3, 1, 3)
ax[2].plot([0, 100], [0, 0], 'k--', linewidth=0.5)
ax[2].plot(r, Phi[2], 'k-.', linewidth=2., label='$h = %s \, GM/c$' % h[2])
ax[2].annotate('ISCO', xy=(6, -1.07*0.057), xycoords='data',
    size=10, ha="left", va="top", color='OrangeRed')
ax[2].scatter([6], [-0.057], 25, c='OrangeRed')
condition = np.logical_and(Phi[2] < -0.025, r > 0)
ax[2].plot(r[condition], [-0.025 for x in r[condition]], 'OrangeRed',
    linestyle='dashed', linewidth=1.5)
ax[2].annotate('capture', xy=(10, -0.018), xycoords='data', size=10,
    ha="center", va="center", color='OrangeRed')
ax[2].set_xlim([0, 100])
ax[2].set_xlabel('$r$ $(GM/c^2)$')
ax[2].xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax[2].set_ylim([-0.1, 0.01])
ax[2].set_ylabel('$\Phi$ $(\mu c^2)$')
ax[2].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
leg = ax[2].legend(loc=4, fontsize=11, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Save the figure.
fig.tight_layout()
plt.savefig('schwarzschild_potentials.pdf')
