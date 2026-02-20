"""
Complete drag polar construction.
Raymer Chapter 12, Section 12.3
"""

import numpy as np


class DragPolar:
    """
    Drag polar: CD = CD0 + K * CL^2   (Raymer Eq. 12.4, uncambered)
            or: CD = CDmin + K * (CL - CLmin_drag)^2  (Raymer Eq. 12.5, cambered)
    """

    def __init__(self, cd0, K, cl_min_drag=0.0):
        """
        Parameters
        ----------
        cd0 : float
            Zero-lift drag coefficient (or CDmin for cambered).
        K : float
            Drag-due-to-lift factor = 1/(pi*A*e).
        cl_min_drag : float
            CL at minimum drag. 0.0 for uncambered (Eq. 12.4),
            small positive value for cambered wings (Eq. 12.5).
        """
        self.cd0 = cd0
        self.K = K
        self.cl_min_drag = cl_min_drag

    def cd(self, cl):
        """
        Total drag coefficient at given CL.
        Raymer Eq. 12.4 or 12.5.
        """
        cl = np.asarray(cl, dtype=float)
        # Raymer Eq. 12.5 (reduces to Eq. 12.4 when cl_min_drag = 0)
        return self.cd0 + self.K * (cl - self.cl_min_drag) ** 2

    def ld(self, cl):
        """Lift-to-drag ratio at given CL."""
        return cl / self.cd(cl)

    def ld_max(self):
        """
        Maximum L/D and the CL at which it occurs.
        For Eq. 12.4: (L/D)max = 1 / (2*sqrt(CD0*K))
                       CL at (L/D)max = sqrt(CD0/K)

        Returns
        -------
        tuple : (ld_max, cl_at_ld_max)
        """
        cl_opt = np.sqrt(self.cd0 / self.K) + self.cl_min_drag
        ld = self.ld(cl_opt)
        return float(ld), float(cl_opt)

    def cl_for_cruise(self, weight, q, s_ref):
        """
        Required CL for level flight: CL = W / (q * S).

        Parameters
        ----------
        weight : float
            Aircraft weight (lbs).
        q : float
            Dynamic pressure (lb/ft^2).
        s_ref : float
            Wing reference area (ft^2).

        Returns
        -------
        float : Required CL
        """
        return weight / (q * s_ref)

    def summary(self):
        """Print a summary of the drag polar."""
        ld, cl = self.ld_max()
        print(f"Drag Polar Summary")
        print(f"  CD0          = {self.cd0:.5f}")
        print(f"  K            = {self.K:.5f}")
        print(f"  CL_min_drag  = {self.cl_min_drag:.4f}")
        print(f"  (L/D)_max    = {ld:.2f}")
        print(f"  CL at L/Dmax = {cl:.4f}")
        print(f"  CD at L/Dmax = {self.cd(cl):.5f}")
