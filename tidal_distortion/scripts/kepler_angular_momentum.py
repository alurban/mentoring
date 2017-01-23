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

# Set array of GW frequency values.
f_GW = np.linspace(1e-4, 4000, 500)

# For each GW frequency, compute the orbital separation in meters,
# the orbital velocity, and the angular momentum of stable circular
# orbits.
a = ( G * (2*M) / (pi * f_GW)**2 )**(1./3)
v = ( 2 * pi * G * M * f_GW )**(1./3) / 299792458.
L = ( G**2 * M**5 / (2 * pi * f_GW) )**(1./3)

# Construct a figure.
fig = plt.figure( figsize=(6, 7.5) )

# Plot the orbital separation as a function of GW frequency.
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(f_GW, a/1e3, 'k', linewidth=2.)
ax1.fill_between(f_GW, 0, 11, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(f_GW, 0, 12, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(f_GW, 0, 13, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.annotate('neutron star radius', xy=(2000, 6), xycoords='data', size=12, ha="center", va="center",
    path_effects=[PE.withStroke(linewidth=3, foreground="w")])
ax1.set_xlim([0, 4000])
ax1.set_ylim([0, 100])
ax1.set_ylabel('orbital separation (km)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax1.get_xticklabels(), visible=False)

# Plot the orbital velocity as a fraction of the speed of light.
ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(f_GW, v, 'k', linewidth=2.)
ax2.set_xlim([0, 4000])
ax2.set_ylim([0, 1])
ax2.set_ylabel('orbital velocity ($c$)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
plt.setp(ax2.get_xticklabels(), visible=False)

# Plot the total angular momentum as a function of frequency.
ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(f_GW, L/1e42, 'k', linewidth=2.)
ax3.set_xlim([0, 4000])
ax3.set_xlabel('gravitational wave frequency (Hz)')
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax3.set_ylim([0, 10])
ax3.set_ylabel('angular momentum (10$^{42}$ J$\cdot$s)')
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

# Save the figure.
plt.savefig('kepler_angular_momentum.pdf')
