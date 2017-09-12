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
rISCO = 6*G*M/c**2  # ISCO radius, in meters
ESun = MSun * c**2  # Rest mass-energy of the Sun


# Function definitions.
def time_to_merger(r):
    """ Returns the predicted time until merger assuming a Schwarzschild background and
        given the following parameter:
            r: the initial orbital separation in meters """
    return (5./256) * (c**5/G**3) * 1./(mu*M**2) * (r**4 - rISCO**4)

def drdt(r):
    """ Returns the rate of orbital separation decay in m/s given the following parameter:
            r: the current orbital separation in meters """
    return -(64./5) * (G**3/c**5) * mu*M**2/r**3

def dEdt(r):
    """ Returns the rate of (binding) energy in J decay given the following parameter:
            r: the current orbital separation in meters """
    return -(32./5) * (G**4/c**5) * mu**2*M**3/r**5

def omega(r):
    """ Returns the angular velocity in rad/s given the following parameter:
            r: the current orbital separation in meters """
    return np.sqrt(G*M / r**3)

def g(f):
    """ Returns the right-hand side of three equations of motion at a single point,
        given the following parameters:
            r: the orbital separation in meters evaluated at the same point """
    r = f[0]  # capture the only scalar that affects change
    separation = drdt(r)
    energy = dEdt(r)
    phi = omega(r)
    return np.array([separation, energy, phi])

def rk4(f, dt):
    """ Returns the estimated integral at a new point using the RK4 method (an extension
         of Simpson's rule) given the following parameters:
            f: an array of function values at the previous point
            dt: the integration step size """
    k1 = g(f)
    k2 = g(f + dt*k1/2)
    k3 = g(f + dt*k2/2)
    k4 = g(f + dt*k3)
    return f + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)


# Set initial conditions.
r0 = 300e3  # initial separation, in meters
Tmerge = time_to_merger(r0)  # expected time until merger
fISCO = 1./np.sqrt(216*pi) * c**3/(G*M)  # ISCO frequency (Hz)
dt = min( (1./(2*fISCO), rISCO/(4*c)) )  # step size, determined
# as either twice the Nyquist step size or 1/4 the light-travel time
# across rISCO -- whichever is smaller
t = np.arange(0, Tmerge + dt, dt)  # time samples


# Define initial conditions.
r = [r0]  # orbital separation in meters
phi = [0]  # polar angle (rad)
E = [-G*M*mu/(2*r0)]  # binding energy (J)
v = [np.sqrt( drdt(r0)**2 + (r0*omega(r0))**2 )]  # orbital velocity (m/s)
pot = [-G*M*mu/r0]  # potential energy (J)
kin = [mu * v[0]**2 / 2]  # relativistic kinetic energy (J)
LGW = [-dEdt(r0)]  # GW luminosity (W)
fGW = [omega(r0)/pi]  # GW frequency (Hz)

# Perform the numerical simulation.
for i in xrange(1, len(t)):
    # perform RK4 integration all at once
    f_old = np.array([r[-1], E[-1], phi[-1]])
    f_new = rk4(f_old, dt)

    # append data to arrays
    r.append( f_new[0] )
    E.append( f_new[1] )
    phi.append( f_new[2] )
    v.append( np.sqrt(drdt(r[-1])**2 + (r[-1]*omega(r[-1]))**2) )
    pot.append( -G*M*mu/r[-1] )
    kin.append( mu*v[-1]**2/2 )
    LGW.append( -dEdt(r[-1]) )
    fGW.append( omega(r[-1])/pi )


# Construct a figure showing orbital parameters.
fig = plt.figure( figsize=(6, 7.5) )

# Plot the orbital separation as a function of time.
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(t, np.array(r)/1000, 'k', linewidth=1.75, label='$r(t)$')
ax1.plot([t.min(), 1.05*t.max()], [rISCO/1000]*2, 'k--', linewidth=1.,
    label='ISCO radius')
ax1.set_xlim([t.min(), 1.05*t.max()])
ax1.set_ylim([0, max(r)/1000])
ax1.set_ylabel('orbital separation (km)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax1.get_xticklabels(), visible=False)
leg = ax1.legend(loc=1, fontsize=10, fancybox=True)

# Plot the GW frequency track as a function of time.
ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(t, fGW, 'k', linewidth=1.5, label=r'$f_{\rm GW}$')
ax2.set_xlim([t.min(), 1.05*t.max()])
ax2.set_ylim([10, 2*fISCO])
ax2.set_yscale('log')
ax2.set_ylabel(r'$f_{\rm GW}$ (Hz)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
plt.setp(ax2.get_xticklabels(), visible=False)
leg = ax2.legend(loc=2, fontsize=10, fancybox=True)

# Plot the orbital velocity as a function of time.
ax3 = ax2.twinx()
ax3.plot(t, np.array(v)/c, 'SandyBrown', linewidth=2., label='$v/c$')
ax3.set_xlim([t.min(), 1.05*t.max()])
ax3.set_ylim([0, 1])
ax3.set_ylabel('velocity ($c$)')
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
leg = ax3.legend(loc=1, fontsize=10, fancybox=True)

# Plot a diagram of the gravitational waveform.
ax4 = fig.add_subplot(3, 1, 3)
ax4.plot([t.min(), 1.05*t.max()], [0]*2, 'DarkSlateGray')
ax4.annotate('number of cycles = %.4g' % (phi[-1]/pi), xy=(10, 1.2),
     xycoords='data', size=10, ha="center", va="center")
rmu, rM = 2*G*mu/c**2, 2*G*M/c**2
D = 3.086e22  # 1 Mpc, in meters
A = np.array([(rmu/D) * (rM/x) for x in r])  # GW amplitude
phase = np.array([2*x for x in phi])  # GW phase
h = A * np.cos(phase)  # gravitational waveform h(t)
n = int( 1 / (2 * fISCO * dt) )  # number of samples in 2*fISCO
ax4.plot(t[::n], h[::n]/1e-20, 'DarkSlateGray')
ax4.set_xlim([t.min(), 1.05*t.max()])
ax4.set_xlabel('time (sec)')
ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))
ax4.set_ylim([-1.05*h.max()/1e-20, 1.05*h.max()/1e-20])
ax4.set_ylabel(r'strain (10$^{-20}$)')
ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))

# Save the figure.
fig.tight_layout()
plt.savefig('inspiral_orbit_parameters.pdf')


# Construct a figure showing the evolution of various energies.
fig = plt.figure( figsize=(6, 5) )

# Plot the energies as a function of time.
ax1 = fig.add_subplot(2, 1, 1)
ax1.plot([t.min(), 1.05*t.max()], [0]*2, 'k--', linewidth=1.)
ax1.plot(t, np.array(pot)/ESun, 'DarkTurquoise', linewidth=1.75, label='potential energy')
ax1.plot(t, np.array(kin)/ESun, 'Chocolate', linewidth=1.75, label='kinetic energy')
ax1.plot(t, np.array(E)/ESun, 'k', linewidth=1.75, label='binding energy')
ax1.set_xlim([t.min(), 1.05*t.max()])
ax1.set_ylim([1.05*min(pot)/ESun, 1.05*max(kin)/ESun])
ax1.set_ylabel(r'energy ($M_{\rm \odot}c^2$)')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
plt.setp(ax1.get_xticklabels(), visible=False)
leg = ax1.legend(loc=3, fontsize=10, fancybox=True)

# Plot the GW luminosity as a function of time.
ax2 = fig.add_subplot(2, 1, 2)
ax2.plot(t, np.array(LGW)/ESun, 'k', linewidth=1.5)
ax2.set_xlim([t.min(), 1.05*t.max()])
ax2.set_xlabel('time (sec)')
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))
ax2.set_ylim([1e-5, 1.5*max(LGW)/ESun])
ax2.set_yscale('log')
ax2.set_ylabel(r'luminosity ($M_{\rm \odot}c^2$/s)')
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))

# Save the figure.
fig.tight_layout()
plt.savefig('inspiral_energies.pdf')


# Construct a figure illustrating the orbits.
fig = plt.figure( figsize=(6, 6) )

# Plot a radial diagram of the simulated orbits.
ax = fig.add_subplot(1, 1, 1, projection='polar')
ax.plot(phi, np.array(r)/1000, 'DodgerBlue', linewidth=1.25)
ax.set_rmax(80)
ax.set_rticks(np.array([20, 80]))
ax.grid(True)
ax.set_xticklabels(['0$^{\circ}$', '45$^{\circ}$', '', '135$^{\circ}$', '180$^{\circ}$',
    '225$^{\circ}$', '270$^{\circ}$', '315$^{\circ}$'])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax.set_title( 'number of orbits = %.4g' % (phi[-1]/(2*pi)) )

# Save the figure.
fig.tight_layout()
plt.savefig('inspiral_diagram.pdf')


# Finally, write a gif illustrating the inspiral orbit.
try:
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    fig = plt.figure( figsize=(6, 6) )
    ax = fig.add_subplot(1, 1, 1)
    # start the visualization in the last 0.05 seconds
    # with a frame rate that samples at thrice fISCO
    start, step = -int(0.06/dt), int(1./(2*fISCO*dt))
    N = len(t[start::step])
    # plot the first frame
    x1, y1 = r[start] * np.cos(phi[start]), r[start] * np.sin(phi[start])
    star1 = ax.scatter(x1/2000, y1/2000, c='Tomato', s=3000, edgecolors='none')
    star2 = ax.scatter(-x1/2000, -y1/2000, c='Tomato', s=3000, edgecolors='none')
    ax.grid(False)
    ax.set_xlim([-50, 50])
    ax.set_xlabel('distance (km)')
    ax.set_ylim([-50, 50])
    ax.set_ylabel('distance (km)')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.set_title('%.2g--%.2g $M_{\odot}$ (slowed 100x)' % (m1/MSun, m2/MSun))
    # draw an inset showing the waveform
    inset_ax = inset_axes(ax, width="40%", height=1., loc=4)
    waveform, = inset_ax.plot([], [], 'DarkSlateGray')
    inset_ax.set_xlim([t[start], t[start] + 0.01])
    inset_ax.xaxis.get_major_formatter().set_useOffset(False)
    inset_ax.set_ylim([-1.05*h.max()/1e-20, 1.05*h.max()/1e-20])
    inset_ax.set_ylabel('strain (10$^{-20}$)')
    inset_ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))
    plt.setp(inset_ax.get_xticklabels(), visible=False)
    def update(i):
        """ Update the neutron star positions, then return a tuple of
            artists that have to be redrawn for this frame. """
        x1 = r[start+step*i] * np.cos(phi[start+step*i])
        y1 = r[start+step*i] * np.sin(phi[start+step*i])
        star1.set_offsets( (x1/2000, y1/2000) )
        star2.set_offsets( (-x1/2000, -y1/2000) )
        waveform.set_data(t[start:start+step*i], h[start:start+step*i]/1e-20)
        maximum = max( (t[start+step*i], t[start] + 0.01) )
        minimum = maximum - 0.01
        inset_ax.set_xlim([minimum, maximum])
        return star1, star2, ax, waveform, inset_ax
    anim = FuncAnimation(fig, update, frames=N, interval=30)
    anim.save('inspiral.gif', dpi=100, writer='imagemagick')
except:
    import sys; sys.exit(0)
