import time
start = time.time()

from lasy.laser import Laser
from lasy.profiles.combined_profile import CombinedLongitudinalTransverseProfile
from lasy.profiles.longitudinal import GaussianLongitudinalProfile
from lasy.profiles.transverse import SuperGaussianTransverseProfile
from lasy.profiles.gaussian_profile import GaussianProfile
from lasy.propagators import AngularSpectrumPropagator
from lasy.optical_elements import Axiparabola
from lasy.utils.laser_utils import get_w0

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
import radialGroupDelay as RGD
import full_field
import sys

dim = "rt"
if sys.argv[1] == "no":
    do_rgd = False
else:
    do_rgd = True
    vf = c/100*int(sys.argv[1])
    print("vf:", vf)


l_w = 10.54e-7
f0 = 7e-2
delta = 2e-3
w = 1e-3
tau = 1.5e-14
E = 6.2
des_dt = 1.39e-16 # PIConGPU Standardwert
w0 = f0 * l_w / w / np.pi
print("w0 =", w0)
print("w/w0 =",w/w0)
if dim == "xyt":
    npoints = (int(2*w/w0), int(2*w/w0), 200)
    npoints_prop = (int(10*w/w0), int(10*w/w0), 200)
    hi = (1.1*w, 1.1*w, 4.5*tau)
    lo = (-1.1*w, -1.1*w, -5.*tau)
elif dim == "rt":
    p_per_r = 1./3
    picpoints_per_p = 2
    print("points in file:", int(1024/picpoints_per_p))
    spacing = 0.1772e-6 * p_per_r * 3 # PIConGPU Standardwert
    npoints = (int(2*w/spacing), 10000)
    cut_frac = 0.3
    hi = (2*w, 9*tau)
    lo = (0., -15*tau)
    offset_frac = hi[1]/4 / (hi[1]-lo[1])
    print(offset_frac)
print(npoints)

print(np.pi*w0**2/l_w)
print(100000*des_dt*c)

print("time:", (time.time()-start)/60, "min")

tps = 0

N = int(sys.argv[2])
zs = np.zeros(N+1)
ts = np.zeros(N+1)
tes = np.zeros(N+1)
ws = np.zeros(N+1)
wes = np.zeros(N+1)
for n in range(N+1):
    profile = CombinedLongitudinalTransverseProfile(l_w, (1,0),
                    GaussianLongitudinalProfile(l_w, tau, 0),
                    SuperGaussianTransverseProfile(w, n_order=6),
                    laser_energy=E)
    #profile = GaussianProfile(l_w, (1,0), E, w, tau, 0.0)
    #propagator = AngularSpectrumPropagator(profile.omega0, "xyt")
    
    laser = Laser(dim, lo, hi, npoints, profile)
    #laser.add_propagator(propagator)
    #laser.show()
    print("time:", (time.time()-start)/60, "min")
    
    axiparabola = Axiparabola(f0, delta, 1.7*w)
    
    if do_rgd:
        def tau_D(r):
            return RGD.tau_D_const_v(r, vf, axiparabola)
        radial_delay = RGD.RadialGroupDelay(tau_D, l_w)
        laser.apply_optics(radial_delay)
        laser.show()
        def ztime(z):
            return (z-axiparabola.f0)/vf + axiparabola.f0/c
    else:
        def ztime(z):
            r2 = (z-axiparabola.f0) / axiparabola.delta*axiparabola.R**2
            return 1/c*(z+r2/2/z-2*axiparabola.R**2/4/axiparabola.delta*np.log(1+axiparabola.delta/axiparabola.f0*r2/axiparabola.R**2))
    print("time:", (time.time()-start)/60, "min")
    
    laser.apply_optics(axiparabola)
    print("time:", (time.time()-start)/60, "min")
    laser.propagate(f0+n*delta/N)
    if n == 0:
        tps = full_field.get_tpeak(laser)
    else:
        ts[n] = full_field.get_tpeak(laser) - tps
    print("t:", ts[n])
    if n > 0 or do_rgd:
        tes[n] = ztime(f0+(n)*delta/N) - (f0+(n)*delta/N) / c
    print("t expect", tes[n])
    ws[n] = get_w0(laser.grid, laser.dim)
    print("w:", ws[n])
    if n > 0:
        wes[n] = l_w*axiparabola.f0/np.pi/axiparabola.R*np.sqrt(axiparabola.delta/((n)*delta/N))
    print("w expect", wes[n])
    zs[n] = (n)*delta/N
    print("z:", zs[n])
    print("time:", (time.time()-start)/60, "min")

if do_rgd:
    name="flfoc_"+sys.argv[1]
else:
    name="axiparabola"

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(zs*1e3, tes*1e15, label="theoretical t-z/c")
ax.plot(zs*1e3, ts*1e15, ".", label="measured t-z/c")
ax.legend()
ax.set_xlabel("$(z-f_0)$/mm")
ax.set_ylabel("$(t-z/c)$/fs")

plt.savefig("flfoc_out/lasy_"+name+"_dir_ts.png")

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(zs[1:]*1e3, wes[1:]*1e6, label="theoretical w")
ax.plot(zs*1e3, ws*1e6, ".", label="measured w")
ax.legend()
ax.set_xlabel("$(z-f_0)$/mm")
ax.set_ylabel("$w/\\mu$m")

plt.savefig("flfoc_out/lasy_"+name+"_dir_ws.png")
