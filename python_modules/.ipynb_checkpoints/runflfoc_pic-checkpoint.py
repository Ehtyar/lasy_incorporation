import ptime

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
des_dt = 2.0564e-16 # spezieller Wert
w0 = f0 * l_w / w / np.pi
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
    #spacing = 0.1772e-6 * p_per_r * 3 # PIConGPU Standardwert
    spacing = 2e-6 * p_per_r
    npoints = (int(2*w/spacing), 7500)
    cut_frac = 0.3
    hi = (2*w, 21*tau)
    lo = (0., -15*tau)
    if do_rgd:
        if vf > c:
            offset_frac = 2/5
        else:
            offset_frac = 1/5
    else:
        offset_frac = 1/5
    print(offset_frac)
print(npoints)

print(np.pi*w0**2/l_w)
print(100000*des_dt*c)

print("steps to focus:",f0/c/des_dt)
print("steps through focus:", delta/c/des_dt)
print("total:", (f0+delta)/c/des_dt)

ptime.ptime()
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
ptime.ptime()
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
ptime.ptime()
laser.apply_optics(axiparabola)
#laser.show()
ptime.ptime()

full_field.laser_to_openPMD(laser, "flfoc_"+sys.argv[1]+"_PIC", Nt=1536, Nx=int(1024/picpoints_per_p), Ny=int(1024/picpoints_per_p), 
                            points_between_r=p_per_r, forced_dt=des_dt, offset_frac=1*offset_frac, file_format="bp", data_step=picpoints_per_p, show=True)
ptime.ptime()
