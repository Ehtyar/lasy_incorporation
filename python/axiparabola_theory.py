"""Theoretical values of a pulse reflected by an axiparabola in the focal region.

Following mainly Ambat et al: "Programmable-trajectory ultrafast flying focus
pulses" in Optics express, 2023.

Supports numpy arrays.

Functions are named (with a few exceptions) as follows: 
[value you want to find][value you know]([value you know], axiparabola[, additional args])

Also provides a different axiparabola object that more directly follows Ambat et al.
"""

import numpy as np
from scipy.constants import c
from scipy.optimize import root_scalar
from lasy.optical_elements.optical_element import OpticalElement

class Axiparabola_Ambat(OpticalElement):
    """Axiparabola according to Ambat et al: "Programmable-trajectory ultrafast flying focus pulses" 
    in Vol. 31, No. 19 / 11 Sep 2023 / Optics Express

    Parameters:
    f0 : float
        Focal length of the axiparabola

    delta : float
        Length of the focal region

    R : float
        Radius of the axiparabola
    """
    def __init__(self, f0, delta, R):
        self.f0 = f0
        self.delta = delta
        self.R = R

    def amplitude_multiplier(self, x, y, omega):
        """
        Return the amplitude multiplier.

        Parameters
        ----------
        x, y, omega : ndarrays of floats
            Define points on which to evaluate the multiplier.
            These arrays need to all have the same shape.

        Returns
        -------
        multiplier : ndarray of complex numbers
            Contains the value of the multiplier at the specified points.
            This array has the same shape as the array omega.
        """
        r = np.sqrt(x**2+y**2)
        T = np.exp(-2j * (omega/c) * self.R**2/4/self.delta*np.log(1+self.delta/self.f0*(r/self.R)**2))
        T[x**2 + y**2 > self.R**2] = 0
        return T

def f(r, axiparabola):
    """The focal distance depends on the radius. Returns the focal distance.

    Parameters:
    r : float
        The radius to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    f : float
        The focal distance at that radius.
    """
    assert np.all(np.less_equal(r, axiparabola.R)), "r needs to be smaller or equal than the radius of the axiparabola R."
    assert np.all(np.greater_equal(r, 0.0)), "r needs to be at least 0"
    return axiparabola.f0 + axiparabola.delta * (r**2 / axiparabola.R**2)

zr = f

def rz(z, axiparabola):
    """The radius at which the pulse is focused to a length z
    
    Parameters:
    z : float
        The distance to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    r : float
        The radius on the axiparabola
    """
    assert np.all(np.less(z, axiparabola.f0 + axiparabola.delta)), "z needs to be smaller than f0 + delta."
    assert np.all(np.greater_equal(z, axiparabola.f0)), "z needs to be at least f0"
    return np.sqrt((z-axiparabola.f0) / axiparabola.delta * axiparabola.R**2)

def sf(r, axiparabola):
    """The sag function of the axiparabola.

    Parameters:
    r : float
        The radius to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    sf : float
        The value of the sag function at the radius.
    """
    assert np.all(np.less_equal(r, axiparabola.R)), "r needs to be smaller or equal than the radius of the axiparabola R."
    assert np.all(np.greater_equal(r, 0.0)), "r needs to be at least 0"
    return axiparabola.R**2/4/axiparabola.delta*np.log(1+axiparabola.delta/axiparabola.f0*(r/axiparabola.R)**2)

sfr = sf

def tfr(r, axiparabola):
    """The time at which the part of the pulse at the radius r is focused.
    
    Parameters:
    r : float
        The radius to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    tf : float
        The focus time.
    """
    assert np.all(np.less_equal(r, axiparabola.R)), "r needs to be smaller or equal than the radius of the axiparabola R."
    assert np.all(np.greater_equal(r, 0.0)), "r needs to be at least 0"
    return 1/c*(f(r, axiparabola)+r**2/2/f(r, axiparabola)-2*sf(r, axiparabola))

def vfr(r, axiparabola):
    """The velocity of the focus coming from the radius r.

    Parameters:
    r : float
        The radius to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    vf : float
        The focus velocity.
    """
    assert np.all(np.less_equal(r, axiparabola.R)), "r needs to be smaller or equal than the radius of the axiparabola R."
    assert np.all(np.greater_equal(r, 0.0)), "r needs to be at least 0"
    return c * (1 + r**2 / 2/f(r, axiparabola)**2)

def vftf(tf, axiparabola):
    """The velocity of the focus at the time tf.

    Parameters:
    tf : float
        The time to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    vf : float
        The focus velocity.
    """
    assert np.all(tf <= tfz((axiparabola.f0 + axiparabola.delta) * (1-1e-7), axiparabola)), "tf needs to be smaller or equal than the time to reach f0 + delta."
    assert np.all(tf >= tfz(axiparabola.f0, axiparabola)), "tf needs to be at least the time to reach f0"
    return vfr(rtf(tf, axiparabola), axiparabola)

def rtf(tf, axiparabola):
    """Numerical approximation of The radius at the axiparabola where the light at the focus at the time tf comes from.

    Parameters:
    tf : float
        The time to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    r : float
        The radius on the axiparabola
    """
    assert np.all(tf <= tfz((axiparabola.f0 + axiparabola.delta) * (1-1e-7), axiparabola)), "tf needs to be smaller or equal than the time to reach f0 + delta."
    assert np.all(tf >= tfz(axiparabola.f0, axiparabola)), "tf needs to be at least the time to reach f0"
    try:
        tf[0]
        mul = True
    except TypeError:
        mul=False
    if mul:
        r = np.zeros_like(tf)
        for n, t in enumerate(tf):
            def func(r):
                return tfr(r, axiparabola) - t
            sol = root_scalar(func, method="bisect", bracket=(0, axiparabola.R)) 
            if not sol.converged:
                raise ValueError("finding the radius did not converge:"+sol.flag)
            r[n] = sol.root
    else:
        def func(r):
            return tfr(r, axiparabola) - tf
        sol = root_scalar(func, method="bisect", bracket=(0, axiparabola.R)) 
        if not sol.converged:
            raise ValueError("finding the radius did not converge:"+sol.flag)
        r = sol.root
    return r

def tfz(z, axiparabola):
    """The time at which the focus will arrive at the distance z from the axiparabola.
    
    Parameters:
    z : float
        The distance to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    tf : float
        The focus time.
    """
    assert np.all(np.less_equal(z, axiparabola.f0 + axiparabola.delta)), "z needs to be smaller or equal than the f0 + delta."
    assert np.all(np.greater_equal(z, axiparabola.f0)), "z needs to be at least f0"
    return tfr(rz(z, axiparabola), axiparabola)

def ztf(tf, axiparabola):
    """The distance along the axis from the axiparabola to the focus point at the time tf.

    Parameters:
    tf : float
        The time to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    Returns:
    z : float
        The focus distance.
    """
    assert np.all(tf <= tfz((axiparabola.f0 + axiparabola.delta) * (1-1e-7), axiparabola)), "tf needs to be smaller or equal than the time to reach f0 + delta."
    assert np.all(tf >= tfz(axiparabola.f0, axiparabola)), "tf needs to be at least the time to reach f0"

    return f(rtf(tf, axiparabola), axiparabola)
    
def wz(z, axiparabola, lambda0):
    """The beam wais w at the distance z from the axiparabola.
    
    Parameters:
    z : float
        The distance to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    lambda0 : float
        The main wavelength of the pulse.

    Returns:
    w : float
        The beam waist.
    """
    assert np.all(np.less_equal(z, axiparabola.f0 + axiparabola.delta)), "z needs to be smaller or equal than the f0 + delta."
    assert np.all(np.greater(z, axiparabola.f0)), "z needs to be more than f0 for this approximation to hold"
    return lambda0*axiparabola.f0/np.pi/axiparabola.R*np.sqrt(axiparabola.delta/(z-axiparabola.f0))

def wtf(tf, axiparabola, lambda0):
    """The beam wais w at the time tf from the axiparabola.
    
    Parameters:
    tf : float
        The time to evaluate at.

    axiparabola : lasy axiparabola object
        The axiparabola for which to calculate this

    lambda0 : float
        The main wavelength of the pulse.

    Returns:
    w : float
        The beam waist.
    """
    assert np.all(tf <= tfz((axiparabola.f0 + axiparabola.delta) * (1-1e-7), axiparabola)), "tf needs to be smaller or equal than the time to reach f0 + delta."
    assert np.all(tf >= tfz(axiparabola.f0, axiparabola)), "tf needs to be at least the time to reach f0"
    return wz(ztf(tf, axiparabola), axiparabola, lambda0)
