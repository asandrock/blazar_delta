"""
Optical depth for γγ absorption in blazar jets

Currently, only the optical depth on broad line region photons is calculated.
The BLR is modelled according to Finke, ApJ 830 (2016) 94.
"""
from __future__ import division, print_function
from numpy cimport ndarray
from numpy import empty
import numpy as np
cimport numpy as np
cdef extern:
    void optical_depth_blr(double *tau_gg_blr, double *eps, double *z_jet,
                           double *redshift, double *L_Hbeta, double *R_Hbeta)

def tau_blr(
    double eps,
    double z_jet,
    double redshift,
    double L_Hbeta,
    double R_Hbeta,
):
    """
    Calculate the absorption optical depth on BLR photons

    Parameters
    ----------
    eps
        Energy of VHE-photon in units of electron rest energy
    z_jet
        Distance from black hole to photon production zone in cm
    redshift
        Redshift of emitting blazar. Set to zero, if you want the energies in
        the reference frame of the host galaxy.
    L_Hbeta
        Luminosity of the Hβ line in erg/sec.
    R_Hbeta
        Radius of the Hβ line emitting shell.

    Output
    ------
    tau_blr
        Optical depth for absorption on BLR photons.
    """

    cdef double tau_blr
    optical_depth_blr(
        &tau_blr,
        &eps,
        &z_jet,
        &redshift,
        &L_Hbeta,
        &R_Hbeta
    )
    return (tau_blr)
