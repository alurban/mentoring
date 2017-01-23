# Imports.
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib.patheffects as PE
from matplotlib import ticker

# Physical constants.
G = 6.67408e-11  # Newton's constant in m^3 / kg / s
MSun = 1.989e30  # Solar mass in kg
M = 1.4 * MSun   # Mass of each neutron star in this example

# Set array of orbital separation values.
a = np.linspace(1e-4, 100, 500) * 1000

# Construct a set of arrays corresponding to effective potentials
# with different angular momenta.
aL = np.array([0, 10, 20, 40, 60])
L = np.sqrt( 0.5 * G * M**3 * aL * 1000 )
Phi = np.array([- G * M**2 / a + x**2 / (M * a**2) for x in L])
colors = ['k--', 'k', 'CornflowerBlue', 'Red', 'DarkSlateGrey']

# Construct a figure.
fig = plt.figure( figsize=(6, 3.25) )

# Plot the total angular momentum as a function of frequency.
ax = fig.add_subplot(1, 1, 1)
ax.plot([0, 100], [0, 0], 'k--', linewidth=0.5)
for x, y, c in zip(aL, Phi, colors):
    ax.plot(a/1000, y/1e45, c, linewidth=2., label='$a_L =$ %s km' % x)
ax.set_xlim([0, 100])
ax.set_xlabel('orbital separation (km)')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax.set_ylim([-30, 10])
ax.set_ylabel('energy (10$^{45}$ J)')
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax.legend(loc=4, fontsize=11, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Save the figure.
fig.tight_layout()
plt.savefig('kepler_potentials.pdf')
