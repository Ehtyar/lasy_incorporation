import time
start = time.time()

from lasy.laser import Laser
from lasy.profiles.combined_profile import CombinedLongitudinalTransverseProfile
from lasy.profiles.longitudinal import GaussianLongitudinalProfile
from lasy.profiles.transverse import SuperGaussianTransverseProfile
from lasy.profiles.gaussian_profile import GaussianProfile
from lasy.propagators import AngularSpectrumPropagator
from lasy.optical_elements import Axiparabola

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c

dim = "rt"

l_w = 10.54e-7
f0 = 7e-2
delta = 5e-3
w = 1e-3
tau = 1.5e-14
E = 6.2
des_dt = 1.39e-16 # PIConGPU Standardwert
w0 = f0 * l_w / w / np.pi
vf = 300000000
print("w0 =", w0)
print("w/w0 =",w/w0)
if dim == "xyt":
    npoints = (int(2*w/w0), int(2*w/w0), 200)
    npoints_prop = (int(10*w/w0), int(10*w/w0), 200)
    hi = (1.1*w, 1.1*w, 4.5*tau)
    lo = (-1.1*w, -1.1*w, -5.*tau)
elif dim == "rt":
    p_per_r = 1
    picpoints_per_p = 2
    print("points in file:", int(1024/picpoints_per_p))
    spacing = 0.1772e-6 * p_per_r * 10  # PIConGPU Standardwert; grob für den Moment
    npoints = (int(5*w/spacing), 800)
    cut_frac = 0.3
    hi = (5*w, 9*tau)
    lo = (0., -15*tau)
    offset_frac = hi[1]/2 / (hi[1]-lo[1])
print(npoints)

print("time:", (time.time()-start)/60, "min")

profile = GaussianProfile(l_w, (1,0), E, w, tau, 0.0)

laser = Laser(dim, lo, hi, npoints, profile)
print("time:", (time.time()-start)/60, "min")

axiparabola = Axiparabola(f0, delta, 1.5*w)
import radialGroupDelay as RGD
def tau_D(r):
    return RGD.tau_D_const_v(r, vf, axiparabola)
radial_delay = RGD.RadialGroupDelay(tau_D, l_w)
laser.apply_optics(radial_delay)
print("time:", (time.time()-start)/60, "min")

laser.apply_optics(axiparabola)
print("time:", (time.time()-start)/60, "min")

laser.propagate(f0)
print("time:", (time.time()-start)/60, "min")

import full_field
full_field.laser_to_openPMD(laser, "fl_foc_rough10_comp", Nt=1024, Nx=int(1024/picpoints_per_p), Ny=int(1024/picpoints_per_p), conversion_safety=1.1,
                            points_between_r=p_per_r, forced_dt=des_dt, offset_frac=offset_frac, file_format="bp", data_step=picpoints_per_p)
print("time:", (time.time()-start)/60, "min")

print(laser.grid.npoints)
dist = des_dt * c * 6000
for n in range(16):
    laser.propagate(dist)
    full_field.laser_to_openPMD(laser, "fl_foc_rough10_comp", Nt=1024, iteration=n+1, Nx=int(1024/picpoints_per_p), Ny=int(1024/picpoints_per_p),
                                forced_dt=des_dt, offset_frac=offset_frac, file_format="bp", data_step=picpoints_per_p, append=True)
    print(f0+(n+1)*dist)
    print("time:", (time.time()-start)/60, "min")

