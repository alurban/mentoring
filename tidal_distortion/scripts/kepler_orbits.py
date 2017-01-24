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


# Function definitions.
def dphidt(a, L):
    """ Returns the derivative of phi(t) given orbital separation (a)
        and orbital angular momentum (L). This is based on the definition
        of L, and the fact that L is conserved. """
    return 2 * L / (M * a**2)

def dadt(a, E, L):
    """ Returns the derivative of a(t) given a itself, the total energy (E),
        and orbital angular momentum (L). This is based on the equation
        for total energy, and the fact that energy is conserved. """
    return np.sqrt(abs( (4./M) * (E + (G * M**2 / a) - ((L/a)**2 / M)) ))

def euler(y, h, rhs):
    """ Performs numerical integration using Euler's method on an equation
        with previous value y, step size h and right-hand side rhs. """
    return y + (h * rhs)


# Set initial conditions.
L = 4e42  # angular momentum in J*s
aL = 2 * L**2 / (G * M**3)  # initial separation in m
Porb = (pi * M * aL**2) / L  # period of a stable circular orbit
dt = Porb / 100.  # step size, determined as 1% of Porb
t = np.arange(0, 15*1.01*Porb, dt)  # sample times
E0 = -(G * M**(5./2) / (2 * L))**2  # total energy of a stable circular orbit
frac = np.array([1, 0.975, 0.8, 0.5, 0])
E = frac * E0  # try orbits with these energies


# Perform the numerical simulation.
a = [[aL] for x in frac]
phi = [[0] for x in frac]
v = [[aL * dphidt(aL, L)] for x in frac]

for i in xrange(len(E)):
    for j in xrange(1, len(t)):
        # get the change in a
        rhs_a = -dadt(a[i][j-1], E[i], L)

        # make sure we're going in the right direction
        if j >= 2 and a[i][j-2] < a[i][j-1]:
            rhs_a *= -1

        # check that the kinetic energy is sane
        a_proposed = euler(a[i][j-1], dt, rhs_a)
        if -G*M**2/a_proposed + (L/a_proposed)**2/M > E[i]:
            rhs_a *= -1

        # carry on with the simulation
        a[i].append( euler(a[i][j-1], dt, rhs_a) )
        rhs_p = dphidt(a[i][j-1], L)
        phi[i].append( euler(phi[i][j-1], dt, rhs_p) )
        v[i].append( np.sqrt(rhs_a**2 + (a[i][j] * dphidt(a[i][j], L))**2) )


# Construct a figure showing orbital parameters.
fig = plt.figure( figsize=(6, 7.5) )
c = ['k--', 'CornflowerBlue', 'Red', 'Silver', 'k']

# Plot the orbital separation as a function of time.
ax1 = fig.add_subplot(3, 1, 1)
for i in xrange(len(E)):
    ax1.plot(t/Porb, np.array(a[i])/1000, c[i], linewidth=2., label='$E =$ %s$E_0$' % frac[i])
ax1.fill_between(t/Porb, 0, 11, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(t/Porb, 0, 12, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(t/Porb, 0, 13, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.annotate('neutron star radius', xy=(1.5, 6), xycoords='data', size=12, ha="center", va="center",
    path_effects=[PE.withStroke(linewidth=3, foreground="w")])
ax1.set_xlim([0, 3])
ax1.set_ylim([0, 100])
ax1.set_ylabel('orbital separation (km)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax1.get_xticklabels(), visible=False)

# Plot the orbital separation as a function of time.
ax2 = fig.add_subplot(3, 1, 2)
for i in xrange(len(E)):
    ax2.plot(t/Porb, np.rad2deg(np.array(phi[i])), c[i], linewidth=2., label='$E =$ %s$E_0$' % frac[i])
ax2.set_xlim([0, 3])
ax2.set_ylim([0, 720])
ax2.set_ylabel(r'$\varphi$ (degrees)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax2.get_xticklabels(), visible=False)
leg = ax2.legend(loc=2, fontsize=10, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Plot the orbital velocity as a function of time.
ax3 = fig.add_subplot(3, 1, 3)
for i in xrange(len(E)):
    ax3.plot(t/Porb, np.array(v[i])/299792458., c[i], linewidth=2., label='$E =$ %s$E_0$' % frac[i])
ax3.set_xlim([0, 3])
ax3.set_xlabel(r'$t/P_{\rm circular}$ (unitless)')
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
ax3.set_ylim([0, 1])
ax3.set_ylabel('velocity ($c$)')
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))

# Save the figure.
fig.tight_layout()
plt.savefig('orbit_parameters.pdf')


# Construct a figure showing sanity checks.
fig = plt.figure( figsize=(6, 7.5) )

# Plot the simulated energies on a diagram as a function of radius.
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot([0, 100], [0, 0], 'k--', linewidth=0.5)
for i in xrange(1, len(E)):
    energy = np.array([(M/4)*dadt(x, E[i], L)**2 - G*M**2/x + (L/x)**2/M for x in a[i]])
    ax1.plot(np.array(a[i])/1000, energy/1e45, c[i], linewidth=2., label='$E =$ %s$E_0$' % frac[i])
a_plot = np.arange(1e-4, 100000, 1)
Phi = -G * M**2 / a_plot + (L/a_plot)**2 / M
ax1.plot(a_plot/1000, Phi/1e45, 'k', linewidth=1.5)
ax1.set_xlim([0, 100])
ax1.set_xlabel('orbital separation (km)')
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax1.set_ylim([-15, 10])
ax1.set_ylabel('energy (10$^{45}$ J)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax1.legend(loc=1, fontsize=10, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Plot the simulated energies as a function of time.
ax2 = fig.add_subplot(3, 1, 2)
for i in xrange(len(E)):
    energy = np.array([(M/4)*dadt(x, E[i], L)**2 - G*M**2/x + (L/x)**2/M for x in a[i]])
    ax2.plot(t/Porb, energy/1e45, c[i], linewidth=2., label='$E =$ %s$E_0$' % frac[i])
ax2.set_xlim([0, 3])
ax2.set_ylim([-15, 5])
ax2.set_ylabel('total energy (10$^{45}$ J)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax2.get_xticklabels(), visible=False)

# Plot the simulated angular momenta as a function of time.
ax3 = fig.add_subplot(3, 1, 3)
for i in xrange(len(E)):
    ang = np.array([0.5*M*(x**2)*dphidt(x, L) for x in a[i]])
    ax3.plot(t/Porb, ang/1e42, c[i], linewidth=2., label='$E =$ %s$E_0$' % frac[i])
ax3.set_xlim([0, 3])
ax3.set_xlabel(r'$t/P_{\rm circular}$ (unitless)')
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
ax3.set_ylim([0, 10])
ax3.set_ylabel('angular momentum (10$^{42}$ J$\cdot$s)')
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

# Save the figure.
fig.tight_layout()
plt.savefig('sanity_checks.pdf')


# Construct a figure illustrating the orbits.
fig = plt.figure( figsize=(6, 6) )

# Plot a radial diagram of the simulated orbits.
ax = fig.add_subplot(1, 1, 1, projection='polar')
for i in xrange(1, len(E)-1):
    ax.plot(phi[i], np.array(a[i])/1000, c[i], linewidth=2., label='$E =$ %s$E_0$' % frac[i])
ax.plot(phi[0][:101], np.array(a[0][:101])/1000, c[0], linewidth=2., label='$E = E_0$')
ax.set_rmax(80)
ax.set_rticks([40, 60, 80])
ax.grid(True)
ax.set_xticklabels(['0$^{\circ}$', '45$^{\circ}$', '90$^{\circ}$', '135$^{\circ}$', '180$^{\circ}$',
    '225$^{\circ}$', '270$^{\circ}$', '315$^{\circ}$'])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax.legend(loc=1, fontsize=10, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Save the figure.
fig.tight_layout()
plt.savefig('orbit_diagram.pdf')


# Construct a figure illustrating the orbits.
fig = plt.figure( figsize=(6, 6) )

# Plot a radial diagram of the simulated orbits.
ax = fig.add_subplot(1, 1, 1, projection='polar')
for i in xrange(1, len(E)-1):
    ax.plot(phi[i], np.array(a[i])/1000, c[i], linewidth=0.5)
ax.plot(phi[-1], np.array(a[-1])/1000, c[-1], linewidth=2., label='$E = 0$')
ax.set_rmax(800)
ax.set_rticks([200, 400, 600, 800])
ax.grid(True)
ax.set_xticklabels(['0$^{\circ}$', '45$^{\circ}$', '90$^{\circ}$', '135$^{\circ}$', '180$^{\circ}$', 
    '225$^{\circ}$', '270$^{\circ}$', '315$^{\circ}$'])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax.legend(loc=1, fontsize=10, fancybox=True)
leg.legendPatch.set_path_effects([PE.withSimplePatchShadow()])

# Save the figure.
fig.tight_layout()
plt.savefig('orbit_diagram_unbound.pdf')
