import numpy as np

from lasy.profiles.transverse.transverse_profile import TransverseProfile

class RectProfile(TransverseProfile):
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