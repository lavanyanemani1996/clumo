#Imports
import numpy as np
import math

from astropy.cosmology import LambdaCDM, w0waCDM
import astropy.units as u
from astropy.cosmology import z_at_value
from astropy import constants as const

import hmf
from hmf import MassFunction
from hmf import transfer

import scipy
from scipy.integrate import simps, trapz

#Functions

#Default Cosmology (WMAP-9)
H0, Ob0, Om0, Ode0, Tcmb0, m_nu = 69.7, 0.0464, 0.1138 / (69.7/100)**2, 0.721, 2.73, 0
default_cosmo = LambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0, Ode0 = Ode0, Tcmb0 = Tcmb0, \
                            m_nu = m_nu*u.eV)

def load_cosmology(H0 = 69.7, Ob0 = 0.0464, Om0 = 0.1138 / (69.7/100)**2, Ode0 = 0.721, \
                    Tcmb0 = 2.73, m_nu = 0):

    new_model = LambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0, Ode0 = Ode0, Tcmb0 = Tcmb0, \
                            m_nu = m_nu*u.eV)

    return new_model

def load_transferfunction(transfer_model_sel = 3, cosmo_model = default_cosmo, sigma_8 = 0.821, \
                            n = 0.972, lnk_min = -4, lnk_max = 2, dlnk = 0.0026, \
                            transfer_params = None, takahashi = True):

    transfer_model_array = ['BBKS', 'BondEfs', 'CAMB', 'EH']
    transfer_model = transfer_model_array[transfer_model_sel]

    t = transfer.Transfer(transfer_model = transfer_model, cosmo_model = cosmo_model, \
                            sigma_8 = sigma_8, n = n, lnk_min = lnk_min, lnk_max = lnk_max, \
                            dlnk = dlnk, transfer_params = transfer_params, takahashi = takahashi)
    return t

def load_massfunction(cosmo_model = default_cosmo, hmf_model_sel = 4, sigma_8 = 0.821, n = 0.972, \
                        delta_h = 500, delta_wrt_sel = 1, Mmin = 13, Mmax = 16, dlog10m = 0.001, \
                        transfer_model_sel = 3, lnk_min = -4, lnk_max = 2, dlnk = 0.0026, \
                        transfer_params = None, takahashi = True):

    hmf_model_array = ['SMT', 'Jenkins', 'Warren', 'Tinker08', 'Tinker10']
    hmf_model = hmf_model_array[hmf_model_sel]

    delta_wrt_array = ['mean', 'crit']
    delta_wrt = delta_wrt_array[delta_wrt_sel]

    transfer_model_array = ['BBKS', 'BondEfs', 'CAMB', 'EH']
    transfer_model = transfer_model_array[transfer_model_sel]

    h = MassFunction(cosmo_model = cosmo_model, hmf_model = hmf_model, sigma_8 = sigma_8, n = n, \
                    delta_h = delta_h, delta_wrt = delta_wrt, Mmin = Mmin, Mmax = Mmax, \
                    dlog10m = dlog10m, transfer_model = transfer_model, lnk_min = lnk_min, \
                    lnk_max = lnk_max, dlnk = dlnk, transfer_params = transfer_params, \
                    takahashi = takahashi)

    return h


def dnbydVdm_z(z, mass, h):
    """
    Calculates for every redshift the mass function.
    """

    s = (np.shape(z)[0], np.shape(mass)[0])
    dndm = np.zeros(s)

    i = 0
    while i < np.shape(z)[0]:
        z_hmf = z[i]
        h.update(z = z_hmf)
        dndm[i] = h.dndm

        i = i+1

    return dndm


def dVbydz_z(z, cosmo):
    """
    Calculates dV/dz [units of Mpc^3 h^(-3)] for width of dz redshift comoving
    volume element for given cosmology cosmo at each given redshift z.
    """

    n_z = np.shape(z)[0]

    D_m = np.zeros(n_z)
    c_km_s = const.c.to('km/s')
    D_h = (c_km_s/cosmo.H0)
    Ok0 = cosmo.Ok0
    E_z = cosmo.efunc(z)

    if np.absolute(Ok0) < 1e-5:
        D_m = [ cosmo.comoving_distance(z[i]) / u.Mpc for i in range(n_z) ]

    elif Ok0 > 0:
        D_m = [ (D_h / u.Mpc) * (1 / math.sqrt(Ok0) ) * \
               math.sinh(math.sqrt(Ok0) * (cosmo.comoving_distance(z[i]) / u.Mpc) / (D_h / u.Mpc)) \
               for i in range(n_z) ]

    elif Ok0 < 0:
        Ok0 = np.absolute(Ok0)
        D_m = [ (D_h / u.Mpc) * (1 / math.sqrt(Ok0)) * \
               math.sin(math.sqrt(Ok0) * (cosmo.comoving_distance(z[i]) / u.Mpc) / (D_h / u.Mpc)) \
               for i in range(n_z) ]

    dVdz_z = (4 * np.pi * (D_h / u.Mpc) * (np.array(D_m)**2) * cosmo.h**3 ) / E_z

    return dVdz_z


def dNbydV_z(zmf, mass):
    """
    Integration wrt mass using trapezium rule at each redshift.
    Returns dN/dV at each redshiftself.
    """

    dNdV_z = np.zeros(np.shape(zmf)[0])

    i = 0
    while i < np.shape(zmf)[0]:
        dNdV_z[i] = trapz(zmf[i], mass)
        i = i+1

    return dNdV_z

def dN_dz_bins(dNdz_z, z_, l):

    J = int(np.shape(dNdz_z)[0] / l)
    dN, dz, z = np.zeros(J), np.zeros(J), np.zeros(J)

    i = 0
    j = 0
    while i < np.shape(dNdz_z)[0]:

        if j == (J-1):
            dN[j] = trapz(dNdz_z[i:],  z_[i:])
            dz[j] = z_[-1] - z_[i]
            z[j] = 0.5 * (z_[-1] + z_[i])

        else:
            dN[j] = trapz(dNdz_z[i:i+l],  z_[i:i+l])
            dz[j] = z_[i+l-1] - z_[i]
            z[j] = 0.5 * (z_[i+l-1] + z_[i])

        j = j+1
        i = i+l

    return dN, dz, z

def create_dN_dV_dz_z(z_min, z_max, n_z, type, h, t, cosmo):

    l = 10
    n_ = n_z * l

    if type == 0:
        z_ = np.linspace(z_min, z_max, n_, endpoint = True)
    else:
        z_ = np.logspace(np.log10(z_min), np.log10(z_max), n_, endpoint = True, base = 10.0)

    mass = h.m
    s = (np.shape(z_)[0], np.shape(mass)[0])
    zmf = np.zeros(s)

    #Creating mass function at each redshift
    zmf = dnbydVdm_z(z_, mass, h)

    #Integrating wrt to mass
    dNdV_z = dNbydV_z(zmf, mass)

    #Volume dV/dz at each redshift
    dVdz_z = dVbydz_z(z_, cosmo)

    #Number of clusters per redshift
    dNdz_z = dNdV_z * dVdz_z

    #Integrating for dN for dz and also dV for the respective bins
    dN, dz, z = dN_dz_bins(dNdz_z, z_, l)
    dV = dVbydz_z(z, cosmo) * dz

    return dN, dV, dz, z
