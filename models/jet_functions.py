import numpy as np
import scipy.integrate as integrate
from . import rates_model as rtmd
from scipy.special import erf

beam_labels = {'GB': '"Gaussian Beam"','SG': '"Uniform Beam"', 'PL': '"Power-Law Beam"', 'DG': '"Double-Gaussian Beam"'}
para_labels = {'GB':{'cosc':r'$\theta_{w}$'}, 'SG':{'cosj':r'$cos\theta_{out}$'}, 'PL':{'cosc':r'$\theta_{in}$','cosj':r'$\theta_{out}$'}, 'DG':{'cosc':r'$\theta_{in}$','cosj':r'$\theta_{out}$','C':r'\mathcal{C}'}}

def top_hat(theta, theta_c, theta_j=np.pi/2., s=1, a=1):
    theta = np.abs(theta)
    theta_j_ = np.where(theta_c<theta_j, theta_c, theta_j)
    return (erf(-(theta-theta_j_)*1000.)/2.)+.5

def Gaussian_beam(theta, theta_c, theta_j=np.pi/2., s=1, a=1):
    theta = np.abs(theta)
    #return np.where(theta<theta_j,np.exp(-0.5*(theta/theta_c)**2),0.)
    return np.exp(-0.5*(theta/theta_c)**2)

def Powerlaw_beam(theta, theta_c, theta_j=np.pi/2., s=2, a=1):
    theta = np.abs(theta)
    #return np.where(theta_c >= theta, 1.0,
    #                 np.where(theta_j >= theta,
    #                 np.power(theta/theta_c,-s),
    #                 0.0))
    return np.where(theta_c >= theta, 1.0, np.power(theta/theta_c,-s))

def Powerlaw_beam_LF(theta, theta_c, theta_j=np.pi/2., s=2, a=1):
    theta = np.abs(theta)
    #return np.where(theta_c >= theta, 1.0,
    #                 np.where(theta_j >= theta,
    #                 np.power(theta/theta_c,-a),
    #                 0.0))
    return np.where(theta_c >= theta, 1.0, np.power(theta/theta_c,-a))

def Double_Gaussian_beam(theta, theta_c, theta_j=np.pi/2., s=0, a=1):
    theta = np.abs(theta)
    return (1.-s)*Gaussian_beam(theta, theta_c) + s*Gaussian_beam(theta, theta_j)

def Double_Gaussian_LF(theta, theta_c, theta_j=np.pi/2., s=1., a=1):
    theta = np.abs(theta)
    out = Double_Gaussian_beam(theta, theta_c, theta_j, s)/Double_Gaussian_beam(theta, theta_c, theta_j, s/a)
    out[np.isnan(out)]=a
    return out

def Paz_func(theta, theta_c):
    theta = np.abs(theta)
    return np.sqrt(1.+(theta/theta_c)**2)

def Paz_PL_beam(theta, theta_c, theta_j=np.pi/2., s=1, a=1):
    theta = np.abs(theta)
    return np.where(theta<theta_j,Paz_func(theta, theta_c)**(-a),0.)

def Paz_PL_LF(theta, theta_c, theta_j=np.pi/2., s=1, a=1):
    theta = np.abs(theta)
    return np.where(theta<theta_j,Paz_func(theta, theta_c)**(-s),0.)

def beam_function(struct):
    if struct=='GJ' or struct=='GB':
        return Gaussian_beam
    elif struct=='PL':
        return Powerlaw_beam
    elif struct=='DG':
        return Double_Gaussian_beam
    elif struct=='PZ':
        return Paz_PL_beam
    elif struct=='TH':
        return top_hat
    else:
        raise ValueError(struct+' is not an implemented jet structure.')

def LF_function(struct):
    if struct=='GJ' or struct=='GB':
        return Gaussian_beam
    elif struct=='PL':
        return Powerlaw_beam_LF
    elif struct=='DG':
        return Double_Gaussian_LF
    elif struct=='PZ':
        return Paz_PL_LF
    elif struct=='TH':
        return top_hat
    else:
        raise ValueError(struct+' is not an implemented jet structure.')


def f_mu(rs, theta_vs, E, theta_c, k_cor, zs, JSfunc, theta_j=np.pi/2., theta_o=np.pi/2., s=1., a=1., LFc=1000., LFfunc=None):
    if not LFfunc:
        LFfunc=JSfunc
    beam = relativistic_beaming(theta_vs, theta_c, LFc, JSfunc, LFfunc=LFfunc, theta_j=theta_j, theta_o=theta_o, s=s, a=a, N_theta=6000)
    return beam*E/(4*np.pi*k_cor*rs**2.)

def beta_JS(theta, LFc, LFfunc, JSparas):
    return (1-((LFc-1)*LFfunc(theta, *JSparas)+1)**-2)**0.5

def BF(theta, theta_v, JSparas, LFfunc, LFc):
    beta = beta_JS(theta, LFc, LFfunc, JSparas)
    A = 1-beta*np.cos(theta)*np.cos(theta_v)
    B = -beta*np.sin(theta)*np.sin(theta_v)
    return (1./(4*np.pi*((LFc-1)*LFfunc(theta, *JSparas)+1)**6))*((np.pi/(A**2 - B**2)**2)*(5*(A/np.sqrt(A**2 - B**2))**3 - 3*(A/np.sqrt(A**2 - B**2))))

def relativistic_beaming_old(theta_v, theta_c, LFc, JSfunc, LFfunc=None, theta_j=np.pi/2., s=1., a=1., N_theta=6000):
    if not LFfunc:
        LFfunc=JSfunc
    theta_max = theta_j
    thetas = np.tile(theta_max*np.power(np.arange(1.,N_theta+1)/N_theta,3), np.expand_dims(theta_v.T,axis=-1).shape).T
    JSparas = (theta_c, theta_j, s, a)
    beam = integrate.simps(np.sin(thetas)*JSfunc(thetas, *JSparas)*BF(thetas, theta_v, JSparas, LFfunc, LFc)/2., x=thetas, axis=0)
    return beam

def relativistic_beaming(theta_v, theta_c, LFc, JSfunc, LFfunc=None, theta_j=np.pi/2., theta_o=np.pi/2., s=1., a=1., N_theta=6000):
    if not LFfunc:
        LFfunc=JSfunc
    theta_max = theta_o
    thetas = np.tile(theta_max*np.power(np.arange(1.,N_theta+1)/N_theta,3), np.expand_dims(theta_v.T,axis=-1).shape).T
    JSparas = (theta_c, theta_j, s, a)
    #beam = integrate.quad(lambda thetas: np.sin(thetas)*JSfunc(thetas, *JSparas)*BF(thetas, theta_v, JSparas, LFfunc, LFc)/2., 0., theta_max)
    beam = integrate.simps(np.sin(thetas)*JSfunc(thetas, *JSparas)*BF(thetas, theta_v, JSparas, LFfunc, LFc)/2., x=thetas, axis=0)
    return beam


def interp_function(func):
    if func=='rates':
        return rtmd.model_vector_fixed
    if func=='rates_int':
        return rtmd.rate_integral
    elif func=='beaming':
        return relativistic_beaming
    else:
        raise ValueError(func+' is not implemented. Must either be "rates" or "beaming".')
