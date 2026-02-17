import numpy as np

from lasy.profiles.transverse.transverse_profile import TransverseProfile
from lasy.profiles.longitudinal.longitudinal_profile import LongitudinalProfile
from lasy.profiles.profile import Profile

class RectTransverseProfile(TransverseProfile):
    """class for a transversal rectangular laser profile.

    Parameters:
    w0 : float
        The width of the rect.

    height : float (optional)
        The height of the rect.
    """
    def __init__(self, w0, height=1.):
        """initialize the RectProfile"""
        super().__init__()
        self.w0 = w0
        self.height = height

    def _evaluate(self, x, y):
        """
        Return the transverse envelope.

        Parameters
        ----------
        x, y : ndarrays of floats
            Define points on which to evaluate the envelope
            These arrays need to all have the same shape.

        Returns
        -------
        envelope: ndarray of complex numbers
            Contains the value of the envelope at the specified points
            This array has the same shape as the arrays x, y
        """
        return self.height * np.heaviside(self.w0**2 - x**2 - y**2, 1.).astype(np.complex128)

class RectLongitudinalProfile(LongitudinalProfile):
    """class for a longitudinal rectangular laser profile.

    Parameters:
    wavelength : float
        The wavelength of the laser.
        
    tau : float
        The length of the rect.

    height : float (optional)
        The height of the rect.
    """
    def __init__(self, wavelength, tau, height=1.):
        super().__init__(wavelength)
        self.tau = tau
        self.height = height

    def evaluate(self, t):
        """
        Return the longitudinal envelope.

        Parameters
        ----------
        t : ndarrays of floats
            Define points on which to evaluate the envelope

        Returns
        -------
        envelope: ndarray of complex numbers
            Contains the value of the longitudinal envelope at the
            specified points. This array has the same shape as the array t.
        """
        return self.height * np.heaviside(np.abs(self.tau) - np.abs(t), 1.).astype(np.complex128)

class RectProfile(Profile):
    """class for a complete cetangular laser profile

    Parameters:

    """
    def __init__(self, wavelength, pol, w0, tau, laser_energy, height=1.):
        super().__init__(wavelength, pol)
        self.w0 = w0
        self.tau = tau
        self.laser_energy = laser_energy
        self.height = height

    def evaluate(self, x, y, t):
        """
        Return the envelope field of the laser.

        Parameters
        ----------
        x, y, t : ndarrays of floats
            Define points on which to evaluate the envelope
            These arrays need to all have the same shape.

        Returns
        -------
        envelope : ndarray of complex numbers
            Contains the value of the envelope at the specified points
            This array has the same shape as the arrays x, y, t
        """
        return (self.height * np.heaviside(self.w0**2 - x**2 - y**2, 1.) * np.heaviside(np.abs(self.tau) - np.abs(t), 1.)).astype(np.complex128)