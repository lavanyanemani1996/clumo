#Imports
import numpy as np
import math
import random
import time
import multiprocessing as mp

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from itertools import cycle
import itertools

from astropy.cosmology import LambdaCDM, w0waCDM
import astropy.units as u
from astropy.cosmology import z_at_value
from astropy import constants as const

import hmf
from hmf import MassFunction
from hmf import transfer

import scipy
from scipy import interpolate, optimize
from scipy.optimize import curve_fit
import scipy.integrate as integrate
from scipy.stats import norm
from scipy.interpolate import interp1d
from scipy.integrate import simps, trapz
from scipy.interpolate import InterpolatedUnivariateSpline as _spline

from functools import partial
import camb
import pynverse

import scaling_relations_module as sm
#Functions


def create_rand_m(dN_bin, z_bin, h):

    mass_sampled_array = np.zeros(shape = (np.size(z_bin), int(np.max(dN_bin))))

    start = time.time()

    n = np.size(z_bin)
    i = 0
    while i < n:

        h.update(z=z_bin[i])
        icdf = _spline((h.ngtm / h.ngtm[0])[::-1], np.log10(h.m[::-1]), k=3)

        N = int(dN_bin[i])
        x = np.random.random(N)
        ind = int(dN_bin[i])

        mass_sampled_array[i, 0:ind] = 10 ** icdf(x)

        i = i+1

    end = time.time()

    print(end-start)

    return mass_sampled_array


def scalingrelations(cosmo, M, sr_rel, sr_rel_params, z_i):

    s = (np.shape(M)[0], np.shape(sr_rel)[0])
    array = np.zeros(s)

    i = 0
    while i < np.shape(sr_rel)[0]:

        #sr_file = 0
        if sr_rel[i] == 0:

            A, B, e = sr_rel_params[:, i]
            if i == 0:
                array[:, i] = sm.LM_0(cosmo, M, z_i, A, B, e)
            elif i == 1:
                array[:, i] = sm.YM_0(cosmo, M, z_i, A, B)
            elif i == 2:
                array[:, i] = sm.MT_0(cosmo, M, z_i, A, B, e)

        #sr_file = 1
        if sr_rel[i] == 1:

            A, B, e = sr_rel_params[:, i]
            if i == 0:
                array[:, i] = sm.LM_1(cosmo, M, z_i, A, B, e)
            elif i == 1:
                array[:, i] = sm.YM_1(cosmo, M, z_i, A, B)
            elif i == 2:
                array[:, i] = sm.MT_1(cosmo, M, z_i, A, B, e)

        #add your own scaling relations

        i = i+1

    return array[:,0], array[:,1], array[:,2]

def cov(R, SD):

    D = np.diag(SD)

    return np.matmul(D, np.matmul(R,D))

def multivariate_new(x):

    cov, r = x
    mean = np.zeros(3)

    return r.multivariate_normal(mean, cov)

def sigma_evol(cosmo, M, z, sig):

    """
    Given M in M_solar h^-1 units, redshift and sigma0 in log10, return sigma(M, z)
    """

    alpha_L = 0.53
    beta_L = -0.5
    M_pivot_L = 1e13 / cosmo.h

    a = sig * ( ( 1 + np.log10(M/M_pivot_L) )**beta_L ) * ( (1+z)**alpha_L )

    return a

def sd_array(cosmo, M, z_i, sr_sd, sr_sd_evol):

    s = (np.shape(M)[0], np.shape(sr_sd)[0])
    srsd_arr = np.zeros(s)

    j = 0
    while j < np.shape(sr_sd_evol)[0]:

        if sr_sd_evol[j] == 1:
            srsd_arr[:, j] = sigma_evol(cosmo, M, z_i, sr_sd[j])
        else:
            srsd_arr[:, j] = np.zeros(np.shape(M)[0]) + sr_sd[j]

        j = j+1

    return srsd_arr

def create_rand_srprop(sr, mass_sampled_array, dN, z, h, cosmo, k = 20):

    pool = mp.Pool(processes = k)

    sr_rel, sr_sd, sr_sd_evol, sr_rel_params, sr_corr = sr[0], sr[1], sr[2], sr[3:6, :], sr[6:, :]

    s = np.shape(mass_sampled_array)
    Lssarr, Yssarr, Tssarr = np.zeros(s), np.zeros(s), np.zeros(s)

    n = np.shape(z)[0]
    i = 0
    while i < n:

        print(i)

        ind = int(dN[i])
        z_i = z[i]
        M = mass_sampled_array[i, 0:ind]

        L, Y, T = np.log10(scalingrelations(cosmo, M, sr_rel, sr_rel_params, z_i))

        srsd_arr = sd_array(cosmo, M, z[i], sr_sd, sr_sd_evol)

        cov_l = [ cov(sr_corr, srsd_arr[i]) for i in range(np.shape(M)[0]) ]
        rngs_l = [ np.random.RandomState(o) for o in range(np.shape(M)[0]) ]

        result = pool.map(multivariate_new, list(zip(cov_l, rngs_l)))
        result_array = np.array(result)

        Lssarr[i, 0:ind] = result_array [:,0] + L
        Yssarr[i, 0:ind] = result_array [:,1] + Y
        Tssarr[i, 0:ind] = result_array [:,2] + T

        i = i+1

    return Lssarr, Yssarr, Tssarr
