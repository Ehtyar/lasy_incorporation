import time
start = time.time()

from lasy.laser import Laser
from lasy.profiles.combined_profile import CombinedLongitudinalTransverseProfile
from lasy.profiles.longitudinal import GaussianLongitudinalProfile
from lasy.profiles.transverse import SuperGaussianTransverseProfile
from lasy.profiles import GaussianProfile
from lasy.propagators import AngularSpectrumPropagator
from lasy.optical_elements import Axiparabola
from lasy.utils.laser_utils import get_w0
from lasy.utils.grid import Grid

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
import radialGroupDelay as RGD
import full_field
from rectProfile import RectTransverseProfile
import axiparabola_theory as axi


dim = "rt"
do_rgd = True
l_w = 10.54e-7
f0 = 7e-2
delta = 2e-3
w = 1e-3
tau = 1.5e-14
E = 6.2
des_dt = 1.39e-16 # PIConGPU Standardwert
des_dt = 7.3443e-17 # spezieller Wert
w0 = f0 * l_w / w / np.pi
vf = 1.02 * c
print("w0 =", w0)
print("w/w0 =",w/w0)
if dim == "xyt":
    npoints = (int(2*w/w0), int(2*w/w0), 200)
    npoints_prop = (int(10*w/w0), int(10*w/w0), 200)
    hi = (1.1*w, 1.1*w, 4.5*tau)
    lo = (-1.1*w, -1.1*w, -5.*tau)
elif dim == "rt":
    p_per_r = 1.0/3
    picpoints_per_p = 2
    print("points in file:", int(1024/picpoints_per_p))
    spacing = 0.1772e-6 * p_per_r * 3 # PIConGPU Standardwert
    npoints = (int(1.1*w/spacing), 7500)
    cut_frac = 0.3
    hi = (1.1*w, 21*tau)
    lo = (0., -15*tau)
    offset_frac = hi[1]/4 / (hi[1]-lo[1])
    print(offset_frac)
print(npoints)

print(np.pi*w0**2/l_w)
print(100000*des_dt*c)

print("time:", (time.time()-start)/60, "min")
profile = CombinedLongitudinalTransverseProfile(l_w, (1,0),
    GaussianLongitudinalProfile(l_w, tau, 0),
    #SuperGaussianTransverseProfile(w, n_order=6),
    RectTransverseProfile(w),
    laser_energy=E)
#profile = GaussianProfile(l_w, (1,0), E, w, tau, 0.0)
#propagator = AngularSpectrumPropagator(profile.omega0, "xyt")

laser = Laser(dim, lo, hi, npoints, profile)
#laser.add_propagator(propagator)
#laser.show()
print("time:", (time.time()-start)/60, "min")
axiparabola = axi.Axiparabola_Ambat(f0, delta, w)
def tau_D(r):
    return RGD.tau_D_const_v(r, vf, axiparabola)

if do_rgd:
    radial_delay = RGD.RadialGroupDelay(tau_D, l_w)
    laser.apply_optics(radial_delay)
    #laser.show()
    def ztime(z):
        return (z-axiparabola.f0)/vf + axiparabola.f0/c
else:
    def ztime(z):
        r2 = (z-axiparabola.f0) / axiparabola.delta*axiparabola.R**2
        return 1/c*(z+r2/2/z-2*axiparabola.R**2/4/axiparabola.delta*np.log(1+axiparabola.delta/axiparabola.f0*r2/axiparabola.R**2))
print("time:", (time.time()-start)/60, "min")
laser.apply_optics(axiparabola)
#laser.show()
print("time:", (time.time()-start)/60, "min")
fig, ax = full_field.show_field(laser, linthresh_frac=0.01, ret_ax=True)
fig.savefig("flying_focus_img/axiparabola.png")
newGrid = Grid(laser.dim, (0., -1e-13), (0.0008, 2.5e-13), npoints, n_azimuthal_modes=1)
laser.propagate(f0, grid_out=newGrid)
#laser.show()
print("time:", (time.time()-start)/60, "min")
fig, ax = full_field.show_field(laser, linthresh_frac=0.01, title="At the focus of the axiparabola", ret_ax=True)
fig.savefig("flying_focus_img/focus.png")
print("w =", get_w0(laser.grid, laser.dim))
tps = full_field.get_tpeak(laser)
print(tps)
N=5
zs = np.zeros(N+1)
ts = np.zeros(N+1)
tes = np.zeros(N+1)
ws = np.zeros(N+1)
wes = np.zeros(N+1)

ws[0] = get_w0(laser.grid, laser.dim)
for n in range(N):
    laser.propagate(delta/N)
    fig, ax = full_field.show_field(laser, ret_ax=True)
    fig.savefig(f"flying_focus_img/focus_step{n+1}.png")
    #laser.show()
    ts[n+1] = full_field.get_tpeak(laser) - tps
    print("t:", ts[n+1])
    tes[n+1] = ztime(f0+(n+1)*delta/N) - (f0+(n+1)*delta/N) / c
    print("t expect", tes[n+1])
    ws[n+1] = get_w0(laser.grid, laser.dim)
    print("w:", ws[n+1])
    wes[n+1] = l_w*axiparabola.f0/np.pi/axiparabola.R*np.sqrt(axiparabola.delta/((n+1)*delta/N))
    print("w expect", wes[n+1])
    zs[n+1] = (n+1)*delta/N
    print("z:", zs[n+1]+f0)
    print("time:", (time.time()-start)/60, "min")

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(zs*1e3, tes*1e15, label="theoretical t-z/c")
ax.plot(zs*1e3, ts*1e15, ".", label="measured t-z/c")
ax.legend()
ax.set_xlabel("$(z-f_0)$/mm")
ax.set_ylabel("$(t-z/c)$/fs")

plt.savefig("flying_focus_img/ts.png")

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(zs[1:]*1e3, wes[1:]*1e6, label="theoretical w")
ax.plot(zs*1e3, ws*1e6, ".", label="measured w")
ax.legend()
ax.set_xlabel("$(z-f_0)$/mm")
ax.set_ylabel("$w/\\mu$m")

plt.savefig("flying_focus_img/ws.png")
