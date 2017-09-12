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

# For each GW frequency, compute the orbital separation in meters
# along with the orbital angular momentum in J*s.
a = ( G * (2*M) / (pi * f_GW)**2 )**(1./3)
L = ( G**2 * M**5 / (2 * pi * f_GW) )**(1./3)

# Compute the binding and tidal forces on a neutron star of radius 11 km.
F_bind_11 = G * M**2 / 11e3**2
F_tide_11 = G * M**2 * ( 1 / (a - 11e3)**2 - 1 / (a + 11e3)**2 )

# Do the same for a neurtron star of radius 12 km.
F_bind_12 = G * M**2 / 12e3**2
F_tide_12 = G * M**2 * ( 1 / (a - 12e3)**2 - 1 / (a + 12e3)**2 )

# And for lucky 13 km.
F_bind_13 = G * M**2 / 13e3**2
F_tide_13 = G * M**2 * ( 1 / (a - 13e3)**2 - 1 / (a + 13e3)**2 )

# Compute the centrifugal force due to orbital motion.
F_c = 2 * L**2 / (M * a**3)

# Construct a figure.
fig = plt.figure( figsize=(6, 5) )

# Plot the tidal and binding forces as a function of frequency.
ax1 = fig.add_subplot(2, 1, 1)
ax1.plot(f_GW, F_tide_11/1e42, 'k', linewidth=2., label='$R = 11$ km')
ax1.plot([0, 4000], [F_bind_11/1e42]*2, 'k')
ax1.plot(f_GW, F_tide_12/1e42, 'k--', linewidth=2., label='$R = 12$ km')
ax1.plot([0, 4000], [F_bind_12/1e42]*2, 'k--')
ax1.plot(f_GW, F_tide_13/1e42, 'k-.', linewidth=2., label='$R = 13$ km')
ax1.plot([0, 4000], [F_bind_13/1e42]*2, 'k-.')
ax1.plot(f_GW, F_c/1e42, 'DodgerBlue', linewidth=2., label='centrifugal')
ax1.set_xlim([0, 4000])
ax1.set_ylim([0, 10])
ax1.set_ylabel('tidal force (10$^{\mathrm{42}}$ N)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax1.get_xticklabels(), visible=False)
leg = ax1.legend(loc=2, fontsize=11, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

ax2 = fig.add_subplot(2, 1, 2)
ax2.plot(f_GW, F_tide_11/1e42 - F_c/1e42, 'k', linewidth=2., label='$R = 11$ km')
ax2.plot([0, 4000], [F_bind_11/1e42]*2, 'k')
ax2.plot(f_GW, F_tide_12/1e42 - F_c/1e42, 'k--', linewidth=2., label='$R = 12$ km')
ax2.plot([0, 4000], [F_bind_12/1e42]*2, 'k--')
ax2.plot(f_GW, F_tide_13/1e42 - F_c/1e42, 'k-.', linewidth=2., label='$R = 13$ km')
ax2.plot([0, 4000], [F_bind_13/1e42]*2, 'k-.')
ax2.set_xlim([0, 4000])
ax2.set_xlabel('gravitational wave frequency (Hz)')
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax2.set_ylim([0, 10])
ax2.set_ylabel('total force (10$^{\mathrm{42}}$ N)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax2.legend(loc=2, fontsize=11, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Save the figure.
plt.savefig('kepler_centrifugal_force.pdf')
