import cluster_module as clumo_s
import sample_module as clumo_p
import xspec_module as clumo_x

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import os, sys
import shutil

def copytree(src, dst, symlinks=False, ignore=None):

    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def create_dir(directory):

    parent_dir = os.getcwd()

    path = os.path.join(parent_dir, directory)

    os.mkdir(path)
    print("Directory '%s' created" %directory)

    path_1 = os.path.join(parent_dir, 'cluster_files')
    path_2 = os.path.join(parent_dir, 'parameter_files')

    path_a =  os.path.join(path, 'cluster_files')
    path_b =  os.path.join(path, 'parameter_files')

    os.mkdir(path_a)
    os.mkdir(path_b)

    copytree(path_1, path_a)
    copytree(path_2, path_b)

    print("Copy of paramter files loaded in directory '%s'. Please change parameters only inside this directory." %directory)

    print("Your machine has {} CPUs / cores".format(mp.cpu_count()))
    
    print("Now in your created directory change parameters as required and then use function clumorun(directoryname, numberofcores) to create cluster sample")

    return None


def clumorun(direc, k):

    directory = str(direc)

    cosmoparams = np.loadtxt(directory+'/parameter_files/cosmoparams.dat')
    H0, Ob0, Om0, Ode0, Tcmb0, m_nu = cosmoparams
    cosmo = clumo_s.load_cosmology(H0, Ob0, Om0, Ode0, Tcmb0, m_nu)

    tfparams = np.loadtxt(directory+'/parameter_files/tfparams.dat')
    n_s, sigma_8, lnk_min, lnk_max, dlnk, trans_model_sel = tfparams
    t = clumo_s.load_transferfunction(int(trans_model_sel), cosmo, sigma_8, n_s, \
                                lnk_min, lnk_max, dlnk)

    mfparams = np.loadtxt(directory+'/parameter_files/mfparams.dat')
    hmf_model_sel, delta_h, delta_wrt_sel, Mmin, Mmax, dlog10m = mfparams
    h = clumo_s.load_massfunction(cosmo, int(hmf_model_sel), sigma_8, n_s, \
                            delta_h, int(delta_wrt_sel), Mmin, Mmax, dlog10m)

    sampleparams = np.loadtxt(directory+'/parameter_files/sampleparams.dat')
    z_min, z_max, n_z, type = sampleparams
    dN, dV, dz, z = clumo_s.create_dN_dV_dz_z(z_min, z_max, n_z, type, h, t, cosmo)

    np.save(directory+'/cluster_files/dN.npy', dN)
    np.save(directory+'/cluster_files/dz.npy', dz)
    np.save(directory+'/cluster_files/z.npy', z)

    sample_m = clumo_p.create_rand_m(dN, z, h)
    np.save(directory+'/cluster_files/sample_m.npy', sample_m)


    sr = np.loadtxt(directory+'/parameter_files/srparams.dat', usecols = (1,2,3))
    Lssarr, Yssarr, Tssarr = clumo_p.create_rand_srprop(sr, sample_m, dN, z, h, \
                                                    cosmo, k)

    np.save(directory+'/cluster_files/Lssarr.npy', Lssarr)
    np.save(directory+'/cluster_files/Yssarr.npy', Yssarr)
    np.save(directory+'/cluster_files/Tssarr.npy', Tssarr)

    xp = np.genfromtxt(directory+'/parameter_files/xspecparams.dat', dtype = 'str')
    rsp = directory+'/parameter_files/eROSITA/rmf01_sdtq.fits'
    arf = directory+'/parameter_files/eROSITA/erosita_7arf.fits'
    flux, cts = clumo_x.create_flux_cts(cosmo, rsp, arf, z, dN, Tssarr, Lssarr, \
                                    sample_m, xp, k, directory)

    np.save(directory+'/cluster_files/flux.npy', flux)
    np.save(directory+'/cluster_files/cts.npy', cts)

    print("Cluster sample created. Check created directory for results")

    return None
