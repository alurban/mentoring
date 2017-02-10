# Imports.
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker


# Physical constants.
G = 6.67408e-11  # Newton's constant in m^3 / kg / s
MSun = 1.989e30  # Solar mass in kg
c = 299792458.   # Speed of light in m/s
m1, m2 = 1.4*MSun, 1.4*MSun  # component masses in kg
mu, M = m1*m2/(m1+m2), m1 + m2  # reduced mass, total mass in kg


# Define the effective potential.
def potential(x, y, a):
    """ Returns the Roche potential evaluated in the equatorial plane
        with the origin at the center of mass, given the parameters:
            x, y:  Cartesian coordinates in meters
            a:     orbital separation in meters
            w:     angular frequency of a circular 2-body orbit """
    a *= 1e3; x *= 1e3; y *= 1e3  # convert distances to km
    dist = m2 / M  # distance from the origin to star 2
    w2 = G*M/a**3  # angular frequency of a circular orbit, squared
    p1 = -G*m1/np.sqrt((x - dist*a)**2 + y**2)
    p2 = -G*m2/np.sqrt((x + (1 - dist)*a)**2 + y**2)
    p3 = -w2 * (x**2 + y**2) / 2
    return p1 + p2 + p3


# Set array of x and y values, in km.
x = np.linspace(-100, 100, 1000)
y = np.linspace(-100, 100, 1000)
a = 75  # orbital separation in km


# Construct a figure.
fig = plt.figure( figsize=(5, 4.75) )
ax = fig.add_subplot(111)

# Construct the 2D potential.
xx, yy = np.meshgrid( x, y )
zz = potential(xx, yy, a)/1e16

# Plot the 11agg rate-redshift posterior.
extent = (x.min(), x.max(), y.min(), y.max())
bound = -G*m2/12e3/1e16  # star 2's self-potential surface
ax.imshow(zz, extent=extent, cmap=cm.YlGnBu, vmin=bound, vmax=0.95*zz.max(), origin='lower')

# Axis 1 labels.
ax.set_xlabel('distance (km)')
ax.set_xlim([x.min(), x.max()])
ax.set_xticks([-100, -50, 0, 50, 100])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))
ax.set_ylabel('distance (km)')
ax.set_ylim([y.min(), y.max()])
ax.set_yticks([-100, -50, 0, 50, 100])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))

# Save the figure.
fig.tight_layout()
plt.savefig('roche_diagram.pdf')


# Construct the first figure.
fig = plt.figure( figsize=(5, 4.75) )
ax = fig.add_subplot(111)

# Plot the 11agg rate-redshift posterior.
levels = np.linspace(0.9*bound, zz.max(), 80)[5::8]
ax.imshow(zz, extent=extent, cmap=cm.YlGnBu, vmin=bound, vmax=0.95*zz.max(), origin='lower')
ax.contour(x, y, zz, levels, colors='MintCream', origin='lower', extent=extent)

# Axis 1 labels.
ax.set_xlabel('distance (km)')
ax.set_xlim([x.min(), x.max()])
ax.set_xticks([-100, -50, 0, 50, 100])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))
ax.set_ylabel('distance (km)')
ax.set_ylim([y.min(), y.max()])
ax.set_yticks([-100, -50, 0, 50, 100])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))

# Save the figure.
fig.tight_layout()
plt.savefig('roche_diagram_with_contours.pdf')
