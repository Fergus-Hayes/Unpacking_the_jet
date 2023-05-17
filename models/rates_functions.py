import numpy as np
import scipy as sp
from scipy.integrate import simps, trapz
from scipy.misc import derivative
from scipy.special import gamma as gamma_fn
from scipy.special import gammainc as gammainc_fn
from scipy.special import gammaincc as gammaincc_fn
from scipy.interpolate import InterpolatedUnivariateSpline
from . import constants as cs
from . import jet_functions as jsfn

def interp_CMV(N=1000):
    z_max=cs.z_max
    zs = np.linspace(1.e-9,z_max,N)
    CMVs = []
    for z in zs:
        CMVs.append(CMV(z))
    CMVs=np.array(CMVs)
    return InterpolatedUnivariateSpline(zs, CMVs, k=1, ext="zeros")

def interp_RF(N=1000):
    z_max=cs.z_max
    zs = np.linspace(1.e-9,z_max,N)
    RFs = []
    for z in zs:
        RFs.append(RF(z))
    RFs=np.array(RFs)
    return InterpolatedUnivariateSpline(zs, RFs, k=1, ext="zeros")

def H_z(z):
    c,H_0,Omega_M,Omega_D = cs.cosmo
    return H_0*np.sqrt((Omega_M*(1.+z)**3)+Omega_D)

def CHW(z):
    fr_a, fr_b, fr_c, fr_d = cs.fr_paras
    c,H_0,Omega_M,Omega_D = cs.cosmo
    return ((fr_a+(fr_b*z))/(1.+np.power(z/fr_c,fr_d))*H_z(z))

def T_z(z, N=1000):
    c,H_0,Omega_M,Omega_D = cs.cosmo
    zz=np.linspace(0.,z,N)
    return np.where(z==0,0.,simps(1./((1.+zz)*H_z(zz)),x=zz, axis=0))

def E_t(t):
    c,H_0,Omega_M,Omega_D = cs.cosmo
    return np.exp(np.log((1.+np.sqrt(Omega_D))/(1.-np.sqrt(Omega_D)))-3*H_0*np.sqrt(Omega_D)*t)

def Z_t(t):
    global cosmo
    c,H_0,Omega_M,Omega_D = cs.cosmo
    return (np.power(Omega_D/Omega_M,1./3)*np.power(np.power((1.+E_t(t))/(1.-E_t(t)),2)-1,1./3))-1

def P_t(t):
    return 1./t

def RF(z, N=1000):
    c,H_0,Omega_M,Omega_D = cs.cosmo
    ts = np.linspace(1.e-3,1.e-1,10000)
    T_inf = ts[np.argmax(Z_t(ts))]
    ts = np.linspace(cs.t_m,T_inf-T_z(z),N)
    zs = Z_t(T_z(z)+ts)
    return simps(CHW(zs)*P_t(ts)/(1+zs), x=ts, axis=0)

def RF_norm():
    return 1./RF(0.)

def LF_BPL(L):
    L_0, alpha, beta, L_min, L_max = cs.lf_paras
    return np.where(L<L_min,0.,np.where(L>L_max,0.,np.where(L<L_0,np.power(L/L_0,-alpha),
               np.power(L/L_0,-beta))))

def LF_Sch(L, L0=np.power(10.,cs.l10_L_0), gamma=cs.gamma, L_min=0.):
    return np.where(L>L_min, np.power(L/L0,-gamma)*np.exp(-L/L0),0.)

def LF_norm(N_L=1000, func=LF_Sch, L_min=np.power(10.,cs.l10_min), L_max=np.power(10.,cs.l10_max), L0=np.power(10.,cs.l10_L_0), gamma=cs.gamma):
    Ls=np.logspace(np.log10(L_min),np.log10(L_max),N_L)
    return 1./(2*simps(np.power(Ls,3./2.)*func(Ls, L0=L0, gamma=gamma),x=Ls))

def LF_SPLEC(Ls, L0=np.power(10.,cs.l10_L_0), gamma=cs.gamma):
    return np.power(L0,gamma)*np.power(Ls,-gamma-1.)*np.exp(-L0/Ls)/gamma_fn(gamma)

def LF_SPLEC_inv(Ls, L0=np.power(10.,cs.l10_L_0), gamma=cs.gamma):
    return np.power(L0,-gamma)*np.power(Ls,gamma-1.)*np.exp(-Ls/L0)/gamma_fn(gamma)

def LF_gamma(Ls, L0, k, Lmin=0., Lmax=np.inf):
    x = Ls/L0
    alpha = k-1
    norm = 1./(gamma_fn(k)*(gammaincc_fn(k,Lmin/L0)-gammaincc_fn(k,Lmax/L0)))
    return norm*np.power(x,alpha)*np.exp(-x)/L0

def spec_band(E):
    A, E_p, E_min, E_max, alpha, beta = cs.spec_paras
    return A*np.where(E<((alpha-beta)*E_p)/(alpha+2.),
               np.power(E/100.,alpha)*np.exp(-((alpha+2.)*E)/(E_p)),
               np.power(E/100.,beta)*np.exp(beta-alpha)
               *np.power(((alpha-beta)*E_p)/(100*(alpha+2.)),alpha-beta))

def k_cor(z, N=1000, E_min=cs.E_min, E_max=cs.E_max, e_min=cs.e_min, e_max=cs.e_max):
    Es = np.linspace(E_min/(1.+z),E_max/(1.+z),N)
    Es_ = np.linspace(e_min,e_max,N)
    return simps(Es*spec_band(Es),x=Es, axis=0)/simps(Es_*spec_band(Es_),x=Es_)

def z_d(rs):
    return (rs/3.086e24)*cs.H_0/cs.c

def I_z(z, N=1000):
    c,H_0,Omega_M,Omega_D = cs.cosmo
    zz=np.linspace(0,z,N)
    return H_0*simps(1./(H_z(zz)),x=zz,axis=0)

def d_L(z):
    c,H_0,Omega_M,Omega_D = cs.cosmo
    return (1+z)*(c/H_0)*I_z(z)*3.086e24

def selection_func(F, detector='BAT'):
    #if detector=='BAT':
    #    Fth = 2.8e-8
    #    return np.where(F<Fth,0., 1.)
    if detector=='BAT':
        a, b, c, d, F0, Fth = cs.SF_BAT
        return np.where(F<Fth,0.,a*(b+(c*F/F0))/(1+(F/(d*F0))))
    #    Fth = 5.e-9
    #    return np.where(F<Fth,0., 1.)
    if detector=='FERMI':
        Fth = cs.F_min
        return np.where(F<Fth,0., 1.)

def CMV(zs):
    c,H_0,Omega_M,Omega_D = cs.cosmo
    return 4*np.pi*(c/H_z(zs))*(d_L(zs)/(3.086e24*(1.+zs)))**2
