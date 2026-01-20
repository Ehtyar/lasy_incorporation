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
from lasy.utils.grid import Grid

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
import radialGroupDelay as RGD
import full_field
import sys

assert len(sys.argv) >= 4, "all arguments need to be set."
cluster = sys.argv[3]
if cluster == "rosi":
    nameplus = "rosi_"
else:
    nameplus = ""
    
lines = []
def printf(string, filename=nameplus+"flfoc_angspec_out/printout"):
    lines.append(string+"\n")
    file = open(filename, "w")
    file.writelines(lines)
    file.close()


printf("running on "+cluster)
printf("settings: v="+sys.argv[1]+"% c; N="+sys.argv[2])

dim = "xyt"
if sys.argv[1] == "no":
    do_rgd = False
else:
    do_rgd = True
    vf = c/100*int(sys.argv[1])
    printf(f"vf: {vf}")

l_w = 10.54e-7
f0 = 7e-2
delta = 2e-3
w = 1e-3
tau = 1.5e-14
E = 6.2
des_dt = 1.39e-16 # PIConGPU Standardwert
w0 = f0 * l_w / w / np.pi
printf(f"w0 = {w0}")
printf(f"w/w0 ={w/w0}")
if dim == "xyt":
    if cluster == "rosi":
        npoints = (330, 330, 300)
    elif cluster == "hemera":
        npoints = (3000, 3000, 2000)
    else:
        raise ValueError("cluster settings only defined for rosi and hemera")
    hi = (2*w, 2*w, 9*tau)
    lo = (-2*w, -2*w, -15*tau)
    offset_frac = hi[1]/4 / (hi[1]-lo[1])
elif dim == "rt":
    p_per_r = 1./3
    picpoints_per_p = 2
    printf(f"points in file:{int(1024/picpoints_per_p)}")
    spacing = 0.1772e-6 * p_per_r * 3 # PIConGPU Standardwert
    npoints = (int(2*w/spacing), 5000)
    cut_frac = 0.3
    hi = (2*w, 9*tau)
    lo = (0., -15*tau)
    offset_frac = hi[1]/4 / (hi[1]-lo[1])
    printf(str(offset_frac))
printf(str(npoints))

printf(str(np.pi*w0**2/l_w))
printf(str(100000*des_dt*c))

printf(f"time: {(time.time()-start)/60} min")

profile = CombinedLongitudinalTransverseProfile(l_w, (1,0),
    GaussianLongitudinalProfile(l_w, tau, 0),
    SuperGaussianTransverseProfile(w, n_order=6),
    laser_energy=E)
#profile = GaussianProfile(l_w, (1,0), E, w, tau, 0.0)
propagator = AngularSpectrumPropagator(profile.omega0, "xyt")

laser = Laser(dim, lo, hi, npoints, profile)
laser.add_propagator(propagator)
#laser.show()
printf(f"time: {(time.time()-start)/60} min")
printf(f"w = {get_w0(laser.grid, laser.dim)}")

axiparabola = Axiparabola(f0, delta, 1.7*w)

if do_rgd:
    def tau_D(r):
        return RGD.tau_D_const_v(r, vf, axiparabola)
    radial_delay = RGD.RadialGroupDelay(tau_D, l_w)
    laser.apply_optics(radial_delay)
    #laser.show()
    def ztime(z):
        return (z-axiparabola.f0)/vf + axiparabola.f0/c
else:
    def ztime(z):
        r2 = (z-axiparabola.f0) / axiparabola.delta*axiparabola.R**2
        return 1/c*(z+r2/2/z-2*axiparabola.R**2/4/axiparabola.delta*np.log(1+axiparabola.delta/axiparabola.f0*r2/axiparabola.R**2))
printf(f"time: {(time.time()-start)/60} min")
printf(f"w = {get_w0(laser.grid, laser.dim)}")

laser.apply_optics(axiparabola)
#laser.show()
printf(f"time: {(time.time()-start)/60} min")
printf(f"w = {get_w0(laser.grid, laser.dim)}")

newGrid = Grid(dim, (-0.5*w, -0.5*w, -5*tau), (0.5*w, 0.5*w, 5*tau), npoints, n_azimuthal_modes=1)
laser.propagate(f0, grid_out=newGrid)
#laser.show()
printf(f"time: {(time.time()-start)/60} min")

printf(f"w = {get_w0(laser.grid, laser.dim)}")

#full_field.laser_to_openPMD(laser, "fl_foc_"+sys.argv[1], Nt=1536, Nx=int(1024/picpoints_per_p), Ny=int(1024/picpoints_per_p), conversion_safety=1.1,
#                            points_between_r=p_per_r, forced_dt=des_dt, offset_frac=1*offset_frac, file_format="bp", data_step=picpoints_per_p)

#laser.show()
printf(f"time: {(time.time()-start)/60} min")

tps = full_field.get_tpeak(laser)
printf(f"{tps}")
N = int(sys.argv[2])
zs = np.zeros(N+1)
ts = np.zeros(N+1)
tes = np.zeros(N+1)
ws = np.zeros(N+1)
wes = np.zeros(N+1)

ws[0] = get_w0(laser.grid, laser.dim)

if do_rgd:
    name="flfoc_"+sys.argv[1]
else:
    name="axiparabola"

fig, ax = full_field.show_field(laser, linthresh_frac=1.,Nr=npoints[0]//2, ret_ax=True)
fig.savefig(nameplus+"flfoc_angspec_out/lasy_"+name+"_focus.png")

for n in range(N):
    laser.propagate(delta/N)
    fig, ax = full_field.show_field(laser, linthresh_frac=1.,Nr=npoints[0]//2, ret_ax=True)
    fig.savefig(nameplus+"flfoc_angspec_out/lasy_"+name+"_step"+str(n)+".png")
    ts[n+1] = full_field.get_tpeak(laser) - tps
    printf(f"t: {ts[n+1]}")
    tes[n+1] = ztime(f0+(n+1)*delta/N) - (f0+(n+1)*delta/N) / c
    printf(f"t expect {tes[n+1]}")
    ws[n+1] = get_w0(laser.grid, laser.dim)
    printf(f"w: {ws[n+1]}")
    wes[n+1] = l_w*axiparabola.f0/np.pi/axiparabola.R*np.sqrt(axiparabola.delta/((n+1)*delta/N))
    printf(f"w expect {wes[n+1]}")
    zs[n+1] = (n+1)*delta/N
    printf(f"z: {zs[n+1]+f0}")
    printf(f"time: {(time.time()-start)/60} min")



fig = plt.figure()
ax = fig.add_subplot()

ax.plot(zs*1e3, tes*1e15, label="theoretical t-z/c")
ax.plot(zs*1e3, ts*1e15, ".", label="measured t-z/c")
ax.legend()
ax.set_xlabel("$z-f_0$/mm")
ax.set_ylabel("$t-z/c$/fs")

plt.savefig(nameplus+"flfoc_angspec_out/lasy_"+name+"_ts.png")

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(zs[1:]*1e3, wes[1:]*1e6, label="theoretical w")
ax.plot(zs*1e3, ws*1e6, ".", label="measured w")
ax.legend()
ax.set_xlabel("$z-f_0$/mm")
ax.set_ylabel("$w/\\mu$m")

plt.savefig(nameplus+"flfoc_angspec_out/lasy_"+name+"_ws.png")
printf("done")
print("done")
