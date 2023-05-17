#!/usr/bin/python

import numpy as np
import scipy as sp
from . import rates_functions as fn_
from . import jet_functions as fn
from . import constants as cs
from scipy.integrate import simps, trapz

def p_z(zs):
    return np.exp(0.3-(zs/8.9))*(1.-(0.41*np.exp(-((zs-1.6)**2)/0.11)))

def p_F(Fs):
    return np.min([(0.27+(Fs/2.e-6)),np.ones(Fs.shape)], axis=0)

def p_Fz(Fs, zs):
    return p_F(Fs)*p_z(zs)

def model_interp(JS_model, theta_c, theta_j=np.pi/2., theta_o=np.pi/2., s=1., a=1., LFc=1000., eps_BNS=1., N_beam=1001, N_z=50, N_theta=100, N_L=1000, beaming=False, interpfunc=None, z_min=1.e-9, z_max=cs.z_max, Ntheta_c=40, L_min=np.power(10.,cs.l10_min), L_max=np.power(10.,cs.l10_max), gamma=cs.gamma, L0=np.power(10.,cs.l10_L_0), detector='BAT'):

    if interpfunc==None:
        return rate_integral_SF(JS_model, theta_c, theta_j=theta_j, theta_o=theta_o, s=s, a=a, LFc=LFc, eps_BNS=eps_BNS, N_beam=N_beam, N_z=N_z, N_theta=N_theta, N_L=N_L, beaming=beaming, interpfunc=interpfunc, z_min=z_min, z_max=z_max, L_min=L_min, L_max=L_max, gamma=gamma, L0=L0, detector=detector)
    if JS_model=='GJ' or JS_model=='GB':
        if beaming:
            #return np.array(interpfunc(np.array([theta_c, LFc, R_BNS, R_NSBH, eps_NSBH]).T))
            return np.array(interpfunc(np.array([theta_c, LFc]).T))
        else:
            raise ValueError('Only beaming implemented so far.')
    elif JS_model=='TH':
        if beaming:
            return np.array(interpfunc(np.array([theta_j, np.log10(LFc)]).T))
    elif JS_model=='PZ':
        if beaming:
            return np.array(interpfunc(np.array([theta_j, np.log10(LFc)]).T))
    elif JS_model=='DG':
        if beaming:
            return np.array(interpfunc(np.array([Ntheta_c*np.power(theta_c/(np.pi/2.),1./3), theta_j, np.log10(s), np.log10(LFc), np.log10(a)]).T))

    else:
        raise ValueError(JS_model+' is not yet implemented.')
        exit()

def rate_integral_SF(JS_model, theta_c, theta_j=np.pi/2., theta_o=np.pi/2., s=1., a=1, LFc=1000., eps_BNS=1., N_beam=6001, N_z=11, N_theta=51, N_L=1001, beaming=False, interpfunc=None, z_min=1.e-9, z_max=1.4, L_min=np.power(10.,cs.l10_min), L_max=np.power(10.,cs.l10_max), gamma=cs.gamma, L0=np.power(10.,cs.l10_L_0), detector='BAT'):

    if JS_model=='PL' and theta_j<theta_c:
        return np.array([0])

    try:
        Ncoord = len(theta_j)
    except:
        Ncoord = 1

    thetas = (np.pi/2.)*np.power(np.linspace(1.,np.ones(Ncoord)*N_theta,N_theta)/N_theta,5)
    zs = np.linspace(z_min, z_max, N_z)

    THETAS, ZS = np.meshgrid(thetas, zs, indexing='ij', sparse=True)

    THETAS = THETAS.reshape((N_theta,1,Ncoord))
    ZS = np.broadcast_to(ZS.T, (Ncoord,N_z,1)).T
    Ls = np.logspace(np.log10(L_min), np.log10(L_max), N_L, axis=0)

    BEAM = fn.relativistic_beaming(THETAS, theta_c, LFc=LFc, JSfunc=fn.beam_function(JS_model), theta_j=theta_j, theta_o=theta_o, s=s, a=a, LFfunc=fn.LF_function(JS_model), N_theta=N_beam)
   
    Fs = BEAM/(4.*np.pi*(fn_.d_L(ZS)**2)*fn_.k_cor(ZS))
    Fs = Ls*np.expand_dims(Fs,axis=0).transpose(1,2,3,0)
    Fs = Fs.transpose(3,0,1,2)

    integrand = (fn_.RF(ZS)/fn_.RF(0.))*(1./(1.+ZS))*fn_.CMV(ZS)*fn_.LF_gamma(np.broadcast_to(Ls,(Ncoord,1,1,N_L)).T, L0=L0, k=gamma, Lmin=L_min, Lmax=L_max)*fn_.selection_func(Fs, detector=detector)*np.sin(THETAS)

    integral__ = simps(integrand, Ls, axis=0)
    integral_ = simps(integral__, THETAS, axis=0)
    integral = simps(integral_, zs, axis=0)

    return integral

def L_iso(JS_model, theta_c, theta_j=np.pi/2., s=1., a=1., LFc=1000., eps_BNS=1., N_beam=1001, N_z=50, N_theta=100, N_L=1000, beaming=False, interpfunc=None, z_min=1.e-9, z_max=cs.z_max, Ntheta_c=40, L_min=np.power(10.,cs.l10_min), L_max=np.power(10.,cs.l10_max), gamma=cs.gamma, L0=np.power(10.,cs.l10_L_0), detector='BAT'):

    if interpfunc==None:
        return L_iso_(JS_model, theta_c, theta_j=theta_j, s=s, a=a, LFc=LFc, eps_BNS=eps_BNS, N_beam=N_beam, N_z=N_z, N_theta=N_theta, N_L=N_L, beaming=beaming, interpfunc=interpfunc, z_min=z_min, z_max=z_max, L_min=L_min, L_max=L_max, gamma=gamma, L0=L0, detector=detector)
    if JS_model=='GJ' or JS_model=='GB':
        if beaming:
            #return np.array(interpfunc(np.array([theta_c, LFc, R_BNS, R_NSBH, eps_NSBH]).T))
            return np.array(interpfunc(np.array([theta_c, LFc]).T))
        else:
            raise ValueError('Only beaming implemented so far.')
    elif JS_model=='TH':
        if beaming:
            return np.array(interpfunc(np.array([theta_j, np.log10(LFc)]).T))
    elif JS_model=='PZ':
        if beaming:
            return np.array(interpfunc(np.array([theta_j, np.log10(LFc)]).T))
    elif JS_model=='DG':
        if beaming:
            return np.array(interpfunc(np.array([Ntheta_c*np.power(theta_c/(np.pi/2.),1./3), theta_j, np.log10(s), np.log10(LFc), np.log10(a)]).T))

    else:
        raise ValueError(JS_model+' is not yet implemented.')
        exit()

def L_iso_(JS_model, theta_c, theta_j=np.pi/2., s=1., a=1, LFc=1000., eps_BNS=1., N_beam=6001, N_z=11, N_theta=51, N_L=1001, beaming=False, interpfunc=None, z_min=1.e-9, z_max=1.4, L_min=np.power(10.,cs.l10_min), L_max=np.power(10.,cs.l10_max), gamma=cs.gamma, L0=np.power(10.,cs.l10_L_0), detector='BAT'):

    if JS_model=='PL' and theta_j<theta_c:
        return np.array([0])

    try:
        Ncoord = len(theta_j)
    except:
        Ncoord = 1

    thetas = (np.pi/2.)*np.power(np.linspace(1.,np.ones(Ncoord)*N_theta,N_theta)/N_theta,5)
    zs = np.linspace(z_min, z_max, N_z)

    THETAS, ZS = np.meshgrid(thetas, zs, indexing='ij', sparse=True)

    THETAS = THETAS.reshape((N_theta,1,Ncoord))
    ZS = np.broadcast_to(ZS.T, (Ncoord,N_z,1)).T
    Ls = np.logspace(np.log10(L_min), np.log10(L_max), N_L, axis=0)

    BEAM = fn.relativistic_beaming(THETAS, theta_c, LFc=LFc, JSfunc=fn.beam_function(JS_model), theta_j=theta_j, s=s, a=a, LFfunc=fn.LF_function(JS_model), N_theta=N_beam)

    Fs = BEAM/(4.*np.pi*(fn_.d_L(ZS)**2)*fn_.k_cor(ZS))
    Fs = Ls*np.expand_dims(Fs,axis=0).transpose(1,2,3,0)
    Fs = Fs.transpose(3,0,1,2)

    integrand = p_Fz(Fs,ZS)*(fn_.RF(ZS)/fn_.RF(0.))*(1./(1.+ZS))*fn_.CMV(ZS)*fn_.LF_gamma(np.broadcast_to(Ls,(Ncoord,1,1,N_L)).T, L0=L0, k=gamma, Lmin=L_min, Lmax=L_max)*fn_.selection_func(Fs, detector=detector)*np.sin(THETAS)

    integrand = integrand.transpose(1,2,0,3)

    integral_ = simps(integrand, thetas[:,0], axis=0)
    integral = simps(integral_, zs, axis=0)

    return Ls,integral

