from lasy.optical_elements.optical_element import OpticalElement
from scipy.constants import c, e, epsilon_0, m_e
from scipy.integrate import solve_ivp
import numpy as np

class RadialGroupDelay(OpticalElement):
    """class that represents a radial group delay echelon.

    Delays the laser puls depending on the distance from the center. 

    See "Programmable-trajectory ultrafast flying focus pulses" by Ambat et al in Vol. 31, No. 19 / 11 Sep 2023 / Optics Express

    
    Parameters
    ----------
    tau_func: function with on parameter r
        This is the function tau_D(r) that defines, how much the puls is supposed to be delayed.
        In seconds given r in meters.

    lambda_0: float (in meter)
        The wavelength for which the echelon should be constructed to preserve wavefronts.
    """
    def __init__(self, tau_func, lambda_0):
        self.tau = tau_func
        self.lambda_0 = lambda_0
        

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
        # calculate tau
        tau = self.tau(np.sqrt(x**2 + y**2))
        # apply phase according to Eq. 11 in Ambat et al (2023)
        return np.exp(2j * omega/c * self.lambda_0/4 * (np.ceil(c*tau/self.lambda_0) + np.floor(c*tau/self.lambda_0)))

def tau_D_const_v(r, v, axiparabola):
    """radial delay function that should (theoretically) result in a focus moving with 
    velocity v when given to the RadialGroupDelay class, if the laser is then passed through
    the delay echelon first and then through the Axiparabola axiparabola. 

    Formula taken from "Axiparabola: a new tool for high-intensity optics" by Kosta Oubrerie et al, 2022.

    use fuctools.partial or similar to define v and axiparabola before giving this to
    RadialGroupDelay."""
    # roughly following Eq. 19 in "Axiparabola: a new tool for high-intensity optics" by Kosta Oubrerie et al 2022 Journal of Optics 24 045503 DOI 10.1088/2040-8986/ac57d2
    A2 = - (v-c)/c**2 * axiparabola.delta/axiparabola.R**2
    A4 = 1/(2*axiparabola.f0**2*c) * ((v-c)/c + 0.5) * axiparabola.delta/axiparabola.R**2
    return A2 * r**2 + A4 * r**4 # + o(r**5)

def tau_D_integrated_palastro(r, v, axiparabola):
    """radial delay function that should (theoretically) result in a focus moving with 
    velocity v when given to the RadialGroupDelay class, if the laser is then passed through
    the delay echelon first and then through the Axiparabola axiparabola.
    
    Formula taken as differential equation from Palastro et al ("Dephasingless Laser Wakefield
    Acceleration" in Physical Review Letters 124)

    use fuctools.partial or similar to define v and axiparabola before giving this to
    RadialGroupDelay."""
    def dz(r):
        return ((1.0 / (4 * axiparabola.f0)) * r * 2
            - (axiparabola.delta / (8 * axiparabola.f0**2 * axiparabola.R**2)) * r**3 * 4
            + axiparabola.delta * (axiparabola.R**2 + 8 * axiparabola.f0 * axiparabola.delta)
            / (96 * axiparabola.f0**4 * axiparabola.R**4) * r**5 * 6)
    def dtau(r, tau):
        return [(c/v - (1-dz(r)**2) / (1+dz(r)**2)) * 2*axiparabola.delta*r/c/axiparabola.R**2,]
    sol = solve_ivp(dtau, (np.min(r), np.max(r)), [0,], t_eval=r)
    return sol.y[0]

def tau_D_integrated_ambat(r, v, axiparabola):
    """radial delay function that should (theoretically) result in a focus moving with 
    velocity v when given to the RadialGroupDelay class, if the laser is then passed through
    the delay echelon first and then through the Axiparabola axiparabola.
    
    Formula taken as differential equation from Ambat et al (see RadialGroupDelay)

    use fuctools.partial or similar to define v and axiparabola before giving this to
    RadialGroupDelay."""
    def f(r):
        return axiparabola.f0 + axiparabola.delta*(r/axiparabola.R)**2
    def dtau(r, tau):
        return [(1-v/c+r**2/2/f(r)**2) * 2*axiparabola.delta*r/c/axiparabola.R**2,]
    sol = solve_ivp(dtau, (np.min(r), np.max(r)), [0,], t_eval=r)
    return sol.y[0]

def v_plasma(lambda0, n_e):
    """returns the velocity required to still achieve vacuum c in the plasma.
    lambda0: wavelength of the laser
    n_e: electron density in the plasma"""
    return c * (1 + (n_e*e**2/epsilon_0/m_e)/2/(2*np.pi*c/lambda0)**2)


