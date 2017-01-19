# Imports.
import numpy as np
from numpy import pi
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as PE
from matplotlib import ticker

# Physical constants.
G = 6.67408e-11  # Newton's constant in m^3 / kg / s
MSun = 1.989e30  # Solar mass in kg
M = 1.4 * MSun   # Mass of each neutron star in this example

# Set array of GW frequency values.
f_GW = np.linspace(0, 3000, 500)

# For each GW frequency, compute the orbital separation in meters.
a = ( G * (2*M) / (pi * f_GW)**2 )**(1./3)

# Compute the binding and tidal forces on a neutron star of radius 11 km.
F_bind_11 = G * M**2 / 11e3**2
F_tide_11 = G * M**2 * ( 1 / (a - 10e3)**2 - 1 / (a + 10e3)**2 )

# Do the same for a neurtron star of radius 12 km.
F_bind_12 = G * M**2 / 12e3**2
F_tide_12 = G * M**2 * ( 1 / (a - 12e3)**2 - 1 / (a + 12e3)**2 )

# And for lucky 13 km.
F_bind_13 = G * M**2 / 13e3**2
F_tide_13 = G * M**2 * ( 1 / (a - 13e3)**2 - 1 / (a + 13e3)**2 )

# Construct a figure.
fig = plt.figure( figsize=(12, 10) )

# Plot the orbital separation as a function of GW frequency.
ax1 = fig.add_subplot(2, 1, 1)
ax1.plot(f_GW, a/1e3, 'k', linewidth=2.5)
ax1.fill_between(f_GW, 0, 11, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(f_GW, 0, 12, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(f_GW, 0, 13, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.annotate('neutron star radius', xy=(1500, 6), xycoords='data', size=20, ha="center", va="center",
    path_effects=[PE.withStroke(linewidth=3, foreground="w")])
ax1.set_xlim([0, 3000])
ax1.set_ylim([0, 100])
ax1.set_ylabel('orbital separation (km)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax1.get_xticklabels(), visible=False)

# Plot the tidal and binding forces as a function of frequency.
ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
ax2.plot(f_GW, F_tide_11/1e42, 'k', linewidth=2.5, label='$R = 11$ km')
ax2.plot([0, 3000], [F_bind_11/1e42]*2, 'k')
ax2.plot(f_GW, F_tide_12/1e42, 'k--', linewidth=2.5, label='$R = 12$ km')
ax2.plot([0, 3000], [F_bind_12/1e42]*2, 'k--')
ax2.plot(f_GW, F_tide_13/1e42, 'k-.', linewidth=2.5, label='$R = 13$ km')
ax2.plot([0, 3000], [F_bind_13/1e42]*2, 'k-.')
ax2.set_xlim([0, 3000])
ax2.set_xlabel('gravitational wave frequency (Hz)')
ax2.set_ylim([0, 10])
ax2.set_ylabel('tidal force ($10^{42}$ N)')
leg = ax2.legend(loc=2, fontsize=20, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Save the figure.
plt.savefig('newtonian_tidal_forces.pdf')
