# Imports.
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib.patheffects as PE
from matplotlib import ticker

# The ImageMagick package is needed for gifs.
# Visit http://www.imagemagick.org.
try:
    from matplotlib.animation import FuncAnimation
except:
    print "Couldn't find gif-making packages; continuing without gifs."


# Physical constants.
G = 6.67408e-11  # Newton's constant in m^3 / kg / s
MSun = 1.989e30  # Solar mass in kg
c = 299792458.   # Speed of light in m/s
m1, m2 = 1.4*MSun, 1.4*MSun  # component masses in kg
mu, M = m1*m2/(m1+m2), m1 + m2  # reduced mass, total mass in kg


# Function definitions.
def Phi(r, h):
    """ Returns the effective potential energy given the following parameters:
            r: the orbital separation in meters
            h: the prescribed orbital angular momentum per unit mass in J*s/kg """
    return -G*mu*M/r + (mu/2)*(h/r)**2 - G*M*mu*h**2/(c**2*r**3)

def v0(r, E, h):
    """ Returns the initial radial velocity constrained by the following parameters:
            r: the initial separation in meters
            E: the prescribed total energy of the orbit in Joules
            h: the prescribed orbital angular momentum per unit mass in J*s/kg """
    return -np.sqrt( (2 / mu) * abs(E - Phi(r, h)) )

def g(f, h):
    """ Returns the right-hand side of three equations of motion at a single point,
        given the following parameters:
            f: an array of function values at the same point
            h: the prescribed orbital angular momentum per unit mass in J*s/kg """
    # Note: we use an effective force to make sure all equations are linear
    #       in the derivatives of r and phi.
    drdt = f[1]
    dvdt = -G*M/f[0]**2 + h**2/f[0]**3 - 3*G*M*h**2/(c**2 * f[0]**4)
    dphidt = h / f[0]**2
    return np.array([drdt, dvdt, dphidt])

def rk4(f, dt, h):
    """ Returns the estimated integral at a new point using the RK4 method (an extension
         of Simpson's rule) given the following parameters:
            f: an array of function values at the previous point
            dt: the integration step size
            h: the prescribed orbital angular momentum per unit mass in J*s/kg """
    k1 = g(f, h)
    k2 = g(f + dt*k1/2, h)
    k3 = g(f + dt*k2/2, h)
    k4 = g(f + dt*k3, h)
    return f + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)


# Set initial conditions.
h = 3.82*G*M/c  # angular momentum, in J*s/kg
rC = (h**2 + h*np.sqrt(h**2 - 12*(G*M/c)**2)) / (2*G*M)  # initial separation, in meters
Porb = 2 * pi * np.sqrt(rC**3/(G*M))  # period of a stable circular orbit
dt = Porb / 1000.  # step size, determined as 1% of Porb
t = np.arange(0, 5*1.01*Porb, dt)  # time samples
EL = Phi(rC, h)  # total energy of a stable circular orbit
frac = np.array([1, 0.99, 0.8, 0.5, 0.2]) # fractions of EL to simulate
v0 = np.array([v0(rC, E, h) for E in frac * EL])  # initial radial velocities


# Perform the numerical simulation.
r = [[rC] for x in frac]
vr = [[x] for x in v0]
vr[-2][0] = -vr[-2][0]  # draw the first capture orbit initially going outward
phi = [[0] for x in frac]
v = [[np.sqrt( y**2 + (h / rC)**2 )] for x, y in zip(frac, v0)]

for i in xrange(len(frac)):
    for j in xrange(1, len(t)):
        # perform the RK4 integration all at once
        f_old = np.array([r[i][-1], vr[i][-1], phi[i][-1]])
        f_new = rk4(f_old, dt, h)

        if f_new[0] > 2*G*M/c**2:  # continue only if outside the event horizon
            # append data to arrays
            r[i].append( f_new[0] )
            vr[i].append( f_new[1] )
            phi[i].append( f_new[2] )
            v[i].append( np.sqrt(vr[i][j]**2 + (h / r[i][j])**2) )


# Construct a figure showing orbital parameters.
fig = plt.figure( figsize=(6, 7.5) )
color = ['k--', 'DodgerBlue', 'Red', 'LimeGreen', 'k']

# Plot the orbital separation as a function of time.
ax1 = fig.add_subplot(3, 1, 1)
for i in xrange(len(frac)):
    ax1.plot(t[:len(r[i])]/Porb, np.array(r[i])/1000, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax1.fill_between(t/Porb, 0, 11, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(t/Porb, 0, 12, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.fill_between(t/Porb, 0, 13, facecolor='Tomato', edgecolor='Tomato', alpha=0.5)
ax1.annotate('neutron star radius', xy=(2.5, 6), xycoords='data', size=10, ha="center", va="center",
    path_effects=[PE.withStroke(linewidth=3, foreground="w")])
ax1.set_xlim([0, 5])
ax1.set_ylim([0, 170])
ax1.set_ylabel('orbital separation (km)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax1.get_xticklabels(), visible=False)

# Plot the orbital separation as a function of time.
ax2 = fig.add_subplot(3, 1, 2)
for i in xrange(len(frac)):
    ax2.plot(t[:len(r[i])]/Porb, np.rad2deg(np.array(phi[i])), color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax2.set_xlim([0, 5])
ax2.set_ylim([0, 1800])
ax2.set_ylabel(r'$\varphi$ (degrees)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax2.get_xticklabels(), visible=False)
leg = ax2.legend(loc=2, fontsize=10, fancybox=True)

# Plot the orbital velocity as a function of time.
ax3 = fig.add_subplot(3, 1, 3)
for i in xrange(len(frac)):
    ax3.plot(t[:len(r[i])]/Porb, np.array(v[i])/c, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax3.set_xlim([0, 5])
ax3.set_xlabel(r'$t/P_{\rm circular}$ (unitless)')
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
ax3.set_ylim([0, 1])
ax3.set_ylabel('velocity ($c$)')
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))

# Save the figure.
fig.tight_layout()
plt.savefig('relativistic_orbit_parameters.pdf')


# Construct a figure showing sanity checks.
fig = plt.figure( figsize=(6, 7.5) )

# Plot the simulated energies on a diagram as a function of radius.
ax1 = fig.add_subplot(3, 1, 1)
for i in xrange(1, len(frac)):
    energy = np.array([(mu/2)*y**2 + Phi(x, h) for x, y in zip(r[i], vr[i])])
    ax1.plot(np.array(r[i])/1000, energy/1e45, color[i], linewidth=2.)
a_plot = np.arange(1e-4, 50*G*M/c**2, 100)
potential = Phi(a_plot, h)
ax1.plot(a_plot/1000, potential/1e45, 'k-.', linewidth=1.5)
ax1.set_xlim([0, 200])
ax1.set_xlabel('orbital separation (km)')
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax1.set_ylim([-6, 0])
ax1.set_ylabel('energy (10$^{45}$ J)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

# Plot the simulated energies as a function of time.
ax2 = fig.add_subplot(3, 1, 2)
for i in xrange(len(frac)):
    energy = np.array([(mu/2)*y**2 + Phi(x, h) for x, y in zip(r[i], vr[i])])
    ax2.plot(t[:len(r[i])]/Porb, energy/1e45, color[i], linewidth=2., label='$E =$ %s$E_L$' % frac[i])
ax2.set_xlim([0, 5])
ax2.set_ylim([-6, 0])
ax2.set_ylabel('total energy (10$^{45}$ J)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax2.get_xticklabels(), visible=False)
leg = ax2.legend(loc=1, fontsize=10, fancybox=True)

# Plot the relative error in energy as a function of time.
ax3 = fig.add_subplot(3, 1, 3)
for i in xrange(len(frac)):
    energy = np.array([(mu/2)*y**2 + Phi(x, h) for x, y in zip(r[i], vr[i])])
    error = np.array([abs((E - frac[i]*EL) / (frac[i] * EL)) for E in energy])
    ax3.plot(t[:len(r[i])]/Porb, error*100., color[i], linewidth=2.)
ax3.set_xlim([0, 5])
ax3.set_xlabel(r'$t/P_{\rm circular}$ (unitless)')
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
ax3.set_ylim([1e-15, 2])
ax3.set_yscale('log')
ax3.set_ylabel('relative error (\%)')

# Save the figure.
fig.tight_layout()
plt.savefig('relativistic_sanity_checks.pdf')


# Construct a figure illustrating the orbits.
fig = plt.figure( figsize=(6, 6) )

# Plot a radial diagram of the simulated orbits.
ax = fig.add_subplot(1, 1, 1, projection='polar')
for i in xrange(1, 3):
    ax.plot(phi[i], np.array(r[i])/1000, color[i], linewidth=1.25, label='$E =$ %s$E_L$' % frac[i])
for i in xrange(3, len(frac)):
    ax.plot(phi[i], np.array(r[i])/1000, color[i], linewidth=2.5, label='$E =$ %s$E_L$' % frac[i])
ax.set_rmax(170)
ax.set_rticks([90, 140])
ax.grid(True)
ax.set_xticklabels(['0$^{\circ}$', '45$^{\circ}$', '90$^{\circ}$', '135$^{\circ}$', '180$^{\circ}$',
    '225$^{\circ}$', '270$^{\circ}$', '315$^{\circ}$'])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
leg = ax.legend(loc=1, fontsize=10, fancybox=True)

# Save the figure.
fig.tight_layout()
plt.savefig('relativistic_orbit_diagram.pdf')


# Finally, write a gif for each closed orbit.
try:
    fig = plt.figure( figsize=(6, 6) )
    ax = fig.add_subplot(1, 1, 1)
    star1 = ax.scatter(r[-2][0]/2000, 0, c='Tomato', s=1000, edgecolors='none')
    star2 = ax.scatter(-r[-2][0]/2000, 0, c='Tomato', s=1000, edgecolors='none')
    ax.grid(False)
    ax.set_xlim([-100, 100])
    ax.set_xlabel('$x$ (km)')
    ax.set_ylim([-100, 100])
    ax.set_ylabel('$y$ (km)')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.set_title('marginal capture orbit')
    def update_1(i):
        """ Update the neutron star positions, then return a tuple of
            artists that have to be redrawn for this frame. """
        x1, y1 = r[-2][50*i] * np.cos(phi[-2][50*i]), r[-2][50*i] * np.sin(phi[-2][50*i])
        x2, y2 = -x1, -y1
        star1.set_offsets( (x1/2000, y1/2000) )
        star2.set_offsets( (x2/2000, y2/2000) )
        return star1, star2, ax
    anim = FuncAnimation(fig, update_1, frames=len(r[-2])/50, interval=40)
    anim.save('capture_orbit.gif', dpi=80, writer='imagemagick')

    fig = plt.figure( figsize=(6, 6) )
    ax = fig.add_subplot(1, 1, 1)
    star1 = ax.scatter(r[1][0]/2000, 0, c='Tomato', s=1000, edgecolors='none')
    star2 = ax.scatter(-r[1][0]/2000, 0, c='Tomato', s=1000, edgecolors='none')
    ax.grid(False)
    ax.set_xlim([-100, 100])
    ax.set_xlabel('$x$ (km)')
    ax.set_ylim([-100, 100])
    ax.set_ylabel('$y$ (km)')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.set_title('elliptical orbit')
    def update_2(i):
        """ Update the neutron star positions, then return a tuple of
            artists that have to be redrawn for this frame. """
        x1, y1 = r[1][50*i] * np.cos(phi[1][50*i]), r[1][50*i] * np.sin(phi[1][50*i])
        x2, y2 = -x1, -y1
        star1.set_offsets( (x1/2000, y1/2000) )
        star2.set_offsets( (x2/2000, y2/2000) )
        return star1, star2, ax
    anim = FuncAnimation(fig, update_2, frames=len(r[1])/50, interval=40)
    anim.save('elliptical_orbit.gif', dpi=80, writer='imagemagick')
except:
    import sys; sys.exit(0)
