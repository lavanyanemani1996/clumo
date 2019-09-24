import numpy as np
from astropy.cosmology import LambdaCDM, w0waCDM
import astropy.units as u

def LM_Schellenberger_2017(cosmo, M, z, A = 0.84, B = 1.35, e = 0.81):
    """
    0	Schellenberger (2017)
    Given masses in M_solar h^{-1} returns [X-ray L in ergs/s] at given redshift z.

	Lx    --> X-ray Luminosity
	h     --> 0.71 used
	E_z   --> E(z) Cosmology function
	e_L   --> redshift evolution of Normalization --> −0.81 +/- 0.12 in M-L (A. Reichert(2018))
	M     --> Mass in SolarM / h

	Norm  --> 0.84 +/- 0.06
	Slope --> 1.35 +/- 0.07

	log10(Lx (erg/s) / (h^(-2) * 1e44)) = Norm*(E_z**e_L) + Slope*log10(M (SolarM) / (h^(-1) * 1e15))
    """

    E_z = cosmo.efunc(z)
    h_relation = 0.71
    M = M * h_relation / cosmo.h
    L = h_relation**(-2) * 1e44 * 10**( A*(E_z**e) + B*np.log10(M/1e15) )

    return L

def YM_Plank_2013(cosmo, M, z, A = -0.175, B = 1.77, e = -2/3):
    """
    0 Plank (2013)
    Given masses in M_solar h^-1 returns Y compton at given z.


	D_A --> Angular diameter distance in Mpc
	E_z --> E(z) Cosmology function
	e_Y --> redshift evolution of Normalization [e_Y = -2/3]
	M   --> Mass in SolarM / h

	Norm --> -0.175 +/- 0.011
    Slope  --> 1.77 +/- 0.06

	E_z^(e_Y) * [ (D_A)^2 (Mpc^2) * Y_500 / 1e-4 ] = 10^(Norm) * [ M (SolarM) / 6e14 ]^(Slope)

    """

    M = M / cosmo.h #Units [M_solar]
    D_A = cosmo.angular_diameter_distance(z).value #Units [Mpc]
    E_z = cosmo.efunc(z)

    Y = 10**(A) * (M / 6e14)**(B) * 1e-4 / ( E_z**e * D_A**2)

    return Y

def MT_Lovisari_2014(cosmo, M, z, A = 0.20, B = 1.71, e = -1.04):
    """
    0	Lovisari (2014)
    Given masses in M_solar h^-1 returns T in keV in z.

	Norm  --> 0.20 +/- 0.02
	Slope --> 1.71 +/- 0.04

	e_T   --> -1.04 +/- 0.07 M-T (A. Reichert, H. Bohringer, R. Fassbender, M. Muhlegger (2018))

    np.log10(Y/C1) = B*np.log10(X/C2) + (A)
    Y = M (M_solar h70^-1)
    X = T (keV)
    C1 = 5 · 10^13 h_70 ^−1 M⊙  C2 = 2 keV

	T = 10 ** ( (1/B) * ( np.log10(M/C1) - ( A * (E_z**e) ) ) )

    """

    h_70 = cosmo.H0.value / 70
    M = M * h_70 / cosmo.h  #M now in units of M_solar h70^-1
    E_z = cosmo.efunc(z)
    C1, C2 = 5e13, 2

    T = 10 ** ( (1/B) * ( np.log10(M/C1) - ( A * (E_z**e) ) ) )

    return T
