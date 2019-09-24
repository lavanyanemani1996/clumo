import numpy as np
from astropy.cosmology import LambdaCDM, w0waCDM
import astropy.units as u

def LM_0(cosmo, M, z, A = 0.84, B = 1.35, e = 0.81):
    """
    0	Format & Default : Schellenberger (2017)
    Given masses in M_solar h^{-1} returns [X-ray L in ergs/s] at given redshift z.

	Lx    --> X-ray Luminosity
	h     --> 0.71 used
	E_z   --> E(z) Cosmology function
	e_L   --> redshift evolution of Normalization --> −0.81 +/- 0.12 in M-L (A. Reichert(2011))
	M     --> Mass in SolarM / h

	Norm  --> 0.84 +/- 0.06
	Slope --> 1.35 +/- 0.07

	log10(Lx (erg/s) / (h^(-2) * 1e44)) = Norm*(E_z**e_L) + Slope*log10(M (SolarM) / (h^(-1) * 1e15))
    """

    E_z = cosmo.efunc(z)
    h_relation = 0.71
    M = M * h_relation / cosmo.h
    A = A*(E_z**e)
    L = h_relation**(-2) * 1e44 * 10**( A + B*np.log10(M/1e15) )

    return L

def LM_1(cosmo, M, z, A = 1.91, B = 0.66, e = -0.81):
    """
    1	Format & Default : Reichert (2011)
        Format also used in Pratt (2009) scaling relations.

    Given masses in M_solar h^{-1} returns [X-ray L in ergs/s] at given redshift z.

	Lx    --> X-ray Luminosity
	E_z   --> E(z) Cosmology function
	e_L   --> redshift evolution of Normalization --> −0.81 +/- 0.12 in M-L (A. Reichert(2011))
	M     --> Mass in SolarM / h

	Norm  --> 1.91 +/- 0.104
	Slope --> 0.66 +/- 0.04

	M = (1e14 Solar mass) * (Norm) * (Lx[1e44 ergs/s])**(Slope)
    """

    M = M / cosmo.h
    E_z = cosmo.efunc(z)
    A = A*(E_z**e)

    L = 1e44 * ((M / (A*1e14)) ** (1/B))

    return L

def YM_0(cosmo, M, z, A = -0.19, B = 1.79, e = -2/3):
    """
    0 Format & Default : Plank (2014)
    Given masses in M_solar h^-1 returns Y * E(z)**e * D_A**2 [Mpc^2] compton at given z.


	D_A --> Angular diameter distance in Mpc
	E_z --> E(z) Cosmology function
	e_Y --> redshift evolution of Normalization [e_Y = -2/3]
	M   --> Mass in SolarM / h

	Norm --> -0.19 +/- 0.011
    Slope  --> 1.79 +/- 0.06

	E_z^(e_Y) * [ (D_A)^2 * Y_500 / 1e-4 ] = 10^(Norm) * [ M (SolarM) / 6e14 ]^(Slope)

    """

    M = M / cosmo.h #Units [M_solar]
    D_A = cosmo.angular_diameter_distance(z).value #Units [Mpc]
    E_z = cosmo.Om0 * (1+z)**3 + cosmo.Ode0

    Y_pivot = 1e-4

    M_pivot = 6e14

    A_SZ = 10**A
    B_SZ = 1.79

    Y = Y_pivot * A_SZ * (M / M_pivot)**B_SZ

    return Y

def YM_1(cosmo, M, z, A = 0.86, B = 1.51, e = -2/3):
    """
    1 Format & Default : Nagarajan (2018)
    Given masses in M_solar h^-1 returns Y * E(z)**e * D_A**2 [Mpc^2] compton at given z.


	E_z --> E(z) Cosmology function
	e_Y --> redshift evolution of Normalization [e_Y = -2/3]
	M   --> Mass in SolarM / h

	Norm --> 0.86 + 0.18 - 0.21
    Slope --> 1.51 + 0.28 - 0.24

	(Y_500 / Y_pivot) * (E_z ** e_Y) = Norm * (M[Solarmass]/M_pivot)**(Slope)

    Y_pivot --> 7.93e-5 [Mpc^2]
    M_pivot --> 7e14 [Solarmass]

    """

    M = M / cosmo.h #Units [M_solar]
    E_z = cosmo.efunc(z)

    D_A = cosmo.angular_diameter_distance(z).value

    Y_pivot = 7.93e-5
    M_pivot = 7e14

    A_SZ = A
    B_SZ = B

    Y = Y_pivot * A_SZ * (M / M_pivot)**B_SZ

    return Y


def MT_0(cosmo, M, z, A = 0.20, B = 1.71, e = -1.04):
    """
    0 Format & Default : Lovisari (2015)
    Given masses in M_solar h^-1 returns T in keV in z.

	Norm  --> 0.20 +/- 0.02
	Slope --> 1.71 +/- 0.04

	e_T   --> -1.04 +/- 0.07 M-T (A. Reichert, H. Bohringer, R. Fassbender, M. Muhlegger (2011))

    Y = M (M_solar h70^-1)
    X = T (keV)
    C1 = 5 · 10^13 h_70 ^−1 M⊙  C2 = 2 keV

    np.log10(Y/C1) = B*np.log10(X/C2) + (A)

    """

    h_70 = cosmo.H0.value / 70
    M = M * h_70 / cosmo.h  #M now in units of M_solar h70^-1
    E_z = cosmo.efunc(z)
    C1, C2 = 5e13, 2

    slope = B
    norm = A * E_z**e

    T = C2 * np.exp( (np.log(M/C1) - norm) / B )

    return T

def MT_1(cosmo, M, z, A = 0.236, B = 1.76, e = -1.04):
    """
    1 Format & Default : Reichert (2011)
    Given masses in M_solar h^-1 returns T in keV in z.

	Norm  --> 0.236 +/- 0.031
	Slope --> 1.76 +/- 0.08

	e_T   --> -1.04 +/- 0.07 M-T (A. Reichert, H. Bohringer, R. Fassbender, M. Muhlegger (2011))

    M[Solarmass] = 1e14[Solarmass] * (Norm) * (T)**(Slope)
    """

    M = M /cosmo.h  #M now in units of M_solar
    E_z = cosmo.efunc(z)
    C1, C2 = 1e14, 1

    slope = B
    norm = np.log( A * E_z**e )

    T = C2 * np.exp( (np.log(M/C1) - norm) / B )

    return T
