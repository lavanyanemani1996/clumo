import xspec
import numpy as np
import multiprocessing as mp
import time
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline as _spline

def xspec_flux(x):

    Temperature, z_, ind_file, El_, Eh_, Z_ = x

    L_factor, norm, abundance = 1e44, 1.0, Z_
    E_norm = [El_, Eh_]

    xspec.Model("apec")
    m1 = xspec.AllModels(1)
    m1.setPars(Temperature, abundance, z_, norm)

    filename = mydirectory+"/parameter_files/xspec/" + str(ind_file) + ".fak"
    fs1 = xspec.FakeitSettings(response, ancillary, fileName = filename,  \
                                exposure = 100000.0)
    xspec.AllData.fakeit(1,fs1,applyStats=False, noWrite = True)

    spec1 = xspec.AllData(1)
    xspec.AllModels.setEnergies(".01 100. 1000 log")

    xspec.AllModels.calcLumin("%f %f %f" %(E_norm[0],E_norm[1],z_))
    L = spec1.lumin[0]
    xspec.AllModels.calcFlux("%f %f" %(E_norm[0],E_norm[1]))
    F = spec1.flux[0]

    new_norm = L_factor / (L * 1e44)
    m1.setPars({4:new_norm})

    xspec.AllModels.calcLumin("%f %f %f" %(E_norm[0],E_norm[1],z_))
    L_new = spec1.lumin[0]
    xspec.AllModels.calcFlux("%f %f" %(E_norm[0],E_norm[1]))
    F_new = spec1.flux[0]

    spec1.ignore("**-%f %f-**" %(E_norm[0],E_norm[1]))
    rate = spec1.rate[3]

    xspec.AllData.clear()

    return [F_new, rate]

def create_it_flux_ctr(cosmo, z, dN, Tssarr, xp, k):


    abund, El, Eh, Z, ni = xp[0], float(xp[1]), float(xp[2]), float(xp[3]), float(xp[4])
    strH0, strOde0 = str(cosmo.H0.value), str(cosmo.Ode0)

    print(strH0 + " 0.0 " + strOde0)
    xspec.Xset.cosmo = strH0 + " 0.0 " + strOde0
    xspec.Xset.abund = abund
    xspec.Xset.chatter = 0

    s = ( np.shape(z)[0], int(np.max(dN)/ni) )
    Tip, fluxip, ctrip = np.zeros(s), np.zeros(s), np.zeros(s)

    ni_copy = ni
    n = np.shape(z)[0]
    i = 0
    pool = mp.Pool(processes = k)
    while i < n:

        print(i)

        ind = int(dN[i])
        Ts = Tssarr[i, 0:ind]
        Tmin, Tmax = np.amin(Ts), np.amax(Ts)

        if dN[i]/ni < ni_copy:
            ni = int(dN[i]/ni_copy)
        else:
            ni = ni_copy

        Tip[i, 0:int(ind/ni)] = np.logspace(Tmin - 0.2, Tmax + 0.1, \
                                            num = int(ind/ni), endpoint = True)
        B = Tip[i, 0:int(ind/ni)]
        A = np.size(B)

        z_i = np.zeros(A) + z[i]
        ind_file = np.linspace(i*A, (i+1)*A, num = A, endpoint = False)

        Ellist = np.zeros(A) + El
        Ehlist = np.zeros(A) + Eh
        Zlist = np.zeros(A) + Z

        result = pool.map(xspec_flux, list(zip(B, z_i, ind_file, \
                                                Ellist, Ehlist, Zlist)))
        result_array = np.array(result)

        fluxip[i, 0:int(ind/ni)] = result_array[:,0]
        ctrip[i, 0:int(ind/ni)] = result_array[:,1]

        ni = ni_copy
        i = i+1

    return Tip, fluxip, ctrip

def create_flux_cts(cosmo, rsp, arf, z, dN, Tssarr, Lssarr, sample_m, xp, k, directory):

    global response
    response = rsp

    global ancillary
    ancillary = arf

    global mydirectory
    mydirectory = directory

    s = np.shape(Tssarr)
    flux_ipt, cts_ipt = np.zeros(s), np.zeros(s)

    Tip, fluxip, ctrip = create_it_flux_ctr(cosmo, z, \
                                            dN, Tssarr, xp, k)
    expt = float(xp[5])

    n = np.size(z)
    ni = float(xp[4])
    ni_copy = ni

    i = 0
    while i < n:

        ind = int(dN[i])
        M = np.log10(sample_m[i, 0:ind])
        T_scatter = 10 ** Tssarr[i, 0:ind]
        L_scatter = 10 ** Lssarr[i, 0:ind]

        if dN[i]/ni < ni_copy:
            ni = int(dN[i]/ni_copy)
        else:
            ni = ni_copy

        T_ip_z = Tip[i, 0:int(ind/ni)]
        fluxes_ip_z = fluxip[i, 0:int(ind/ni)]
        ctr_ip_z = ctrip[i, 0:int(ind/ni)]

        f1 = _spline(T_ip_z, fluxes_ip_z)
        f2 = _spline(T_ip_z, ctr_ip_z)

        flux_ipt[i, 0:ind] = f1(T_scatter) * L_scatter / 1e44
        ctr_ipt = f2(T_scatter) * L_scatter / 1e44
        cts_ipt[i, 0:ind] = ctr_ipt * expt

        ni = ni_copy
        i = i+1

    return flux_ipt, cts_ipt
