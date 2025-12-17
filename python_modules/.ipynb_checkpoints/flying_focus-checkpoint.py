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
import time

import radialGroupDelay as RGD
import full_field

def flying_focus(w_l, w, tau, E, f0, delta, vf, hi, lo, npoints, dim="rt", save=True,
                 show_after=False, show_after_steps=3, show_after_frac=1/3, npoints_file=(192,192,1024), 
                 p_per_r=1, des_dt=1.39e-16, print_time=True, show_laser=True, filename="fl_foc", 
                 offset_frac=0, save_wtest=False):
    """does the flying focus calculation..."""
    if print_time:
        start = time.time()

    # generate the laser
    print("generate the laser")
    profile = GaussianProfile(w_l, (1,0), E, w, tau, 0.0)
    laser = Laser(dim, lo, hi, npoints, profile)
    if print_time:
        print("time:", (time.time()-start)/60, "min")

    print("apply the RGD")
    # generate the axiparabola
    axiparabola = Axiparabola(f0, delta, w)
    # generate and apply the RGD
    def tau_D(r):
        return RGD.tau_D_const_v(r, vf, axiparabola)
    radial_delay = RGD.RadialGroupDelay(tau_D, l_w)
    laser.apply_optics(radial_delay)
    if show_laser:
        laser.show()
    if save_wtest:
        laser.write_to_file(write_dir="wtest_out/w"+str(w).replace(".","_"), file_prefix="RGD")
    if print_time:
        print("time:", (time.time()-start)/60, "min")

    print("apply the axiparabola")
    # apply the axiparabola
    laser.apply_optics(axiparabola)
    if show_laser:
        laser.show()
    if save_wtest:
        laser.write_to_file(write_dir="wtest_out/w"+str(w).replace(".","_"), file_prefix="Axiparabola")
    if print_time:
        print("time:", (time.time()-start)/60, "min")

    print("prpoagate the laser")
    # propagate the laser to the focus
    laser.propagate(f0)
    if show_laser:
        laser.show()
    if save_wtest:
        laser.write_to_file(write_dir="wtest_out/w"+str(w).replace(".","_"), file_prefix="Propagate")
    if print_time:
        print("time:", (time.time()-start)/60, "min")

    # save the laser field
    if save:
        full_field.laser_to_openPMD(laser, filename, Nt=npoints_file[-1], Nx=npoints_file[0], Ny=npoints_file[1],
                                    points_between_r=p_per_r, forced_dt=des_dt, offset_frac=offset_frac)
        if show_laser:
            laser.show()
        if save_wtest:
            laser.write_to_file(write_dir="wtest_out/w"+str(w).replace(".","_"), file_prefix="saved")
        if print_time:
            print("time:", (time.time()-start)/60, "min")

    # show what happens next
    if show_after:
        print("propagating further")
        for n in range(show_after_steps):
            laser.propagate(delta*show_after_frac)
            if show_laser:
                laser.show()
            if save_wtest:
                laser.write_to_file(write_dir="wtest_out/w"+str(w).replace(".","_"), file_prefix="after_step"+str(n))
            print(f0+(n+1)*delta*show_after_frac)
            if print_time:
                print("time:", (time.time()-start)/60, "min")

    return laser



if __name__ == "__main__":
    import sys
    
    
    l_w = 10.54e-7
    f0 = 70e-2
    delta = 2e-2
    w = 10e-3
    tau = 1.5e-14
    E = 6.2
    des_dt = 1.39e-16 # PIConGPU Standardwert
    w0 = f0 * l_w / w / np.pi
    vf = 300000000
    print("w0 =", w0)
    print("w/w0 =",w/w0)

    p_per_r = 2
    picpoints_per_p = 2
    print("points in file:", int(1024/picpoints_per_p))
    print("approximate file size: "+str(4*int(1024/picpoints_per_p)**2/1024/1024+0.2)+" GB")
    spacing = 0.1772e-6 * p_per_r * picpoints_per_p # PIConGPU Standardwert

    for n in sys.argv[1]:
        npoints = (int(1.1*w*(int(n)+1)/spacing), 1600)
        cut_frac = 0.3
        hi = (1.1*w*(int(n)+1), 9*tau)
        lo = (0., -39*tau)
        offset_frac = hi[1]/2 / (hi[1]-lo[1])
        print(offset_frac)
        print(npoints)
        print("w =", w*(int(n)+1))
        flying_focus(l_w, w*(int(n)+1), tau, E, f0, delta, vf, hi, lo, npoints, show_after=True, filename="wtest"+n,
                                  npoints_file=(int(1024/picpoints_per_p),int(1024/picpoints_per_p),1024), p_per_r=p_per_r, offset_frac=offset_frac,
                    save_wtest=True)
        plt.show()

        