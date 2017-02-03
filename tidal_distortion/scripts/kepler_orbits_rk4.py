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
c = 299792458.   # Speed of light in m/s


# Function definitions.
def Phi(r, L):
    """ Returns the effective potential energy given the following parameters:
            r: the orbital separation in meters
            L: the prescribed orbital angular momentum """
    return -G*M**2 / r + (L/r)**2 / M

def v0(r, E, L):
    """ Returns the initial radial velocity constrained by the following parameters:
            r: the initial separation in meters
            E: the prescribed total energy of the orbit in Joules
            L: the prescribed orbital angular momentum in J*s """
    return -np.sqrt( (4 / M) * abs(E - Phi(r, L)) )

def g(f, L):
    """ Returns the right-hand side of three equations of motion at a single point,
        given the following parameters:
            f: an array of function values at the same point
            L: the prescribed orbital angular momentum in J*s """
    # Note: we use an effective force to make sure all equations are linear
    #       in the derivatives of r and phi.
    drdt = f[1]
    dvdt = -2*G*M/f[0]**2 + 4*L**2/(M**2 * f[0]**3)
    dphidt = 2 * L / (M * f[0]**2)
    return np.array([drdt, dvdt, dphidt])

def rk4(f, h, L):
    """ Returns the estimated integral at a new point using the RK4 method (an extension
         of Simpson's rule) given the following parameters:
            f: an array of function values at the previous point
            h: the integration step size
            L: the prescribed orbital angular momentum in J*s """
    k1 = g(f, L)
    k2 = g(f + h*k1/2, L)
    k3 = g(f + h*k2/2, L)
    k4 = g(f + h*k3, L)
    return f + (h/6) * (k1 + 2*k2 + 2*k3 + k4)


# Set initial conditions.
L = 4e42  # angular momentum, in J*s
rL = 2 * L**2 / (G * M**3)  # initial separation, in meters
Porb = (pi * M * rL**2) / L  # period of a stable circular orbit
dt = Porb / 100.  # step size, determined as 1% of Porb
t = np.arange(0, 5*1.01*Porb, dt)  # time samples
EL = -(G * M**(5./2) / (2 * L))**2  # total energy of a stable circular orbit
frac = np.array([1, 0.99, 0.8, 0.5, 0]) # fractions of EL to simulate
v0 = np.array([v0(rL, E, L) for E in frac * EL])  # initial radial velocities


# Perform the numerical simulation.
r = [[rL] for x in frac]
vr = [[x] for x in v0]
phi = [[0] for x in frac]
v = [[np.sqrt( y**2 + (2 * L / (M * rL))**2 )] for x, y in zip(frac, v0)]

for i in xrange(len(frac)):
    for j in xrange(1, len(t)):
        # perform the RK4 integration all at once
        f_old = np.array([r[i][j-1], vr[i][j-1], phi[i][j-1]])
        f_new = rk4(f_old, dt, L)

        # append data to arrays
        r[i].append( f_new[0] )
        vr[i].append( f_new[1] )
        phi[i].append( f_new[2] )
        v[i].append( np.sqrt(vr[i][j]**2 + (2 * L / (M * r[i][j]))**2) )


# Construct a figure showing orbital parameters.
fig = plt.figure( figsize=(6, 7.5) )
color = ['k--', 'CornflowerBlue', 'Red', 'Silver', 'k']

# Plot the orbital separation as a function of time.
ax1 = fig.add_subplot(3, 1, 1)
for i in xrange(len(frac)):
    ax1.plot(t/Porb, np.array(r[i])/1000, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax1.fill_between(t/Porb, 0, 11, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(t/Porb, 0, 12, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(t/Porb, 0, 13, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.annotate('neutron star radius', xy=(2.5, 6), xycoords='data', size=12, ha="center", va="center",
    path_effects=[PE.withStroke(linewidth=3, foreground="w")])
ax1.set_xlim([0, 5])
ax1.set_ylim([0, 100])
ax1.set_ylabel('orbital separation (km)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax1.get_xticklabels(), visible=False)

# Plot the orbital separation as a function of time.
ax2 = fig.add_subplot(3, 1, 2)
for i in xrange(len(frac)):
    ax2.plot(t/Porb, np.rad2deg(np.array(phi[i])), color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax2.set_xlim([0, 5])
ax2.set_ylim([0, 1800])
ax2.set_ylabel(r'$\varphi$ (degrees)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax2.get_xticklabels(), visible=False)
leg = ax2.legend(loc=2, fontsize=10, fancybox=True)

# Plot the orbital velocity as a function of time.
ax3 = fig.add_subplot(3, 1, 3)
for i in xrange(len(frac)):
    ax3.plot(t/Porb, np.array(v[i])/c, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax3.set_xlim([0, 5])
ax3.set_xlabel(r'$t/P_{\rm circular}$ (unitless)')
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
ax3.set_ylim([0, 1])
ax3.set_ylabel('velocity ($c$)')
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))

# Save the figure.
fig.tight_layout()
plt.savefig('kepler_orbit_parameters.pdf')


# Construct a figure showing sanity checks.
fig = plt.figure( figsize=(6, 7.5) )

# Plot the simulated energies on a diagram as a function of radius.
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot([0, 100], [0, 0], 'k--', linewidth=0.5)
for i in xrange(1, len(frac)):
    energy = np.array([(M/4)*y**2 + Phi(x, L) for x, y in zip(r[i], vr[i])])
    ax1.plot(np.array(r[i])/1000, energy/1e45, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
a_plot = np.arange(1e-4, 100000, 100)
potential = Phi(a_plot, L)
ax1.plot(a_plot/1000, potential/1e45, 'k-.', linewidth=1.5)
ax1.set_xlim([0, 100])
ax1.set_xlabel('orbital separation (km)')
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax1.set_ylim([-15, 10])
ax1.set_ylabel('energy (10$^{45}$ J)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax1.legend(loc=1, fontsize=10, fancybox=True)

# Plot the simulated energies as a function of time.
ax2 = fig.add_subplot(3, 1, 2)
for i in xrange(len(frac)):
    energy = np.array([(M/4)*y**2 + Phi(x, L) for x, y in zip(r[i], vr[i])])
    ax2.plot(t/Porb, energy/1e45, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax2.set_xlim([0, 5])
ax2.set_ylim([-15, 5])
ax2.set_ylabel('total energy (10$^{45}$ J)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax2.get_xticklabels(), visible=False)

# Plot the simulated angular momenta as a function of time.
ax3 = fig.add_subplot(3, 1, 3)
for i in xrange(len(frac) - 1):
    energy = np.array([(M/4)*y**2 + Phi(x, L) for x, y in zip(r[i], vr[i])])
    error = np.array([abs((E - frac[i]*EL) / (frac[i] * EL)) for E in energy])
    ax3.plot(t/Porb, error*100., color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax3.set_xlim([0, 5])
ax3.set_xlabel(r'$t/P_{\rm circular}$ (unitless)')
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
ax3.set_ylim([1e-15, 1e-1])
ax3.set_yscale('log')
ax3.set_ylabel('relative error (\%)')

# Save the figure.
fig.tight_layout()
plt.savefig('kepler_sanity_checks.pdf')


# Construct a figure illustrating the orbits.
fig = plt.figure( figsize=(6, 6) )

# Plot a radial diagram of the simulated orbits.
ax = fig.add_subplot(1, 1, 1, projection='polar')
for i in xrange(len(frac)-1):
    ax.plot(phi[i], np.array(r[i])/1000, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax.set_rmax(80)
ax.set_rticks([40, 60, 80])
ax.grid(True)
ax.set_xticklabels(['0$^{\circ}$', '45$^{\circ}$', '90$^{\circ}$', '135$^{\circ}$', '180$^{\circ}$',
    '225$^{\circ}$', '270$^{\circ}$', '315$^{\circ}$'])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax.legend(loc=1, fontsize=10, fancybox=True)

# Save the figure.
fig.tight_layout()
plt.savefig('kepler_orbit_diagram.pdf')
