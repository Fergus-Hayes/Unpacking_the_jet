import numpy as np
import pandas as pd
from . import constants as cs
from . import rates_functions as rtfn
from .jet_functions import relativistic_beaming, f_mu, beam_function, LF_function
from .rates_model import model_interp
from scipy.special import gammaln, gamma, gammainc, gammaincc, gammaincinv, gammainccinv, logsumexp, erf
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import os, sys, time, bilby, inspect, pickle
np.random.seed(int(time.time()))

background = 0.#5900.*6.1e-10

def logpow(x, m):
    """
    Calculates log(x**m) since m*log(x) will fail when m, x = 0.
    """
    # return m * log(x)
    return np.where(np.equal(x, 0), np.where(np.equal(m, 0), 0.0, -np.inf), m * np.log(x))

def factln(n):
    return gammaln(n + 1)

def log_sum_exp(x, axis=None):
    '''
    Returns the logorithm of the sum of the exponents of x.
    '''
    x_max = x.max(axis=axis, keepdims=True)
    return x_max + np.log(np.exp(x-x_max).sum(axis=axis))

def logp(y, mu, sig):
    return (np.log(1./(mu.shape[0]*np.sqrt(2.*np.pi*sig**2.)))+np.log(2./(1.-erf(-y/(sig*np.sqrt(2.)))))+logsumexp(-0.5*((y-mu)/(sig))**2.,axis=0))
    
def logq1(y, mu, sig):
    return np.sum(logp(y, mu, sig))

def logq2(Nz_obs, Nz):
    return (np.where(np.logical_and(Nz >= 0, Nz_obs >= 0), logpow(Nz, Nz_obs) - factln(Nz_obs) - Nz, -np.inf))

def trunc_gamma(shape, scale=1.0, a=0., b=np.inf, size=None):
   return scale*gammainccinv(shape,(gammaincc(shape,a/scale)-(np.random.uniform(0,1,size=size)*(gammaincc(shape,a/scale)-gammaincc(shape,b/scale)))))

def theta_j_cond(theta, LFc, L0, JSfunc, LFfunc, JSparas, T=1):
    return (19.1*((L0*JSfunc(theta, *JSparas))/(T*1.e51))**(1./6)) - ((LFc-1)*LFfunc(theta, *JSparas)+1)

def joint_func_interp(x, theta_c, theta_j, s, a, LFc, R_BNS, R_NSBH, eps_NSBH, L0, alpha):
    rs, theta_vs, JS, rate_paras, b, k_cor, zs, JSfunc, LFfunc, interpfunc, rates_, GWtrig_  = x
    Lmin = L0*1.e-3
    Lmax = L0*1.e2
    JSparas = (theta_c, theta_j, s, a)
    theta_o = np.pi/2.
    if JS!='TH':
        L0_mu = np.mean(trunc_gamma(alpha, scale=L0, a=Lmin, b=Lmax, size=10000))
        theta_o = optimize.fsolve(theta_j_cond, x0=0.44, args=(LFc, L0_mu, JSfunc, LFfunc, JSparas))
        if theta_o<theta_c or theta_o>np.pi/2.:
            theta_o = np.pi/2.

    if rates_:
        Nz = rate_paras*(R_BNS+(eps_NSBH*R_NSBH))*model_interp(JS, theta_c, theta_j=theta_j, theta_o=theta_o, s=s, a=a, LFc=LFc, N_beam=6000, beaming=True, interpfunc=interpfunc, L0=L0, gamma=alpha, L_min=Lmin, L_max=Lmax)
        Nz = np.array([Nz])
    else:
        Nz = False
    if GWtrig_:
        N,S = rs.shape
        E = trunc_gamma(alpha, scale=L0, a=Lmin, b=Lmax, size=(N,S))
        mu = f_mu(rs, theta_vs, E, theta_c, k_cor, zs, JSfunc, LFfunc=LFfunc, theta_j=theta_j, theta_o=theta_o, s=s, a=a, LFc=LFc).T + b
    else:
        mu = False

    return Nz, mu

def convert_x_y_to_z(parameters):
    """
    Function to convert between sampled parameters and constraint parameter.

    Parameters
    ----------
    parameters: dict
        Dictionary containing sampled parameter values, 'x', 'y'.

    Returns
    -------
    dict: Dictionary with constraint parameter 'z' added.
    """
    converted_parameters = parameters.copy()
    converted_parameters['z'] = parameters['theta_j'] - parameters['theta_c']
    return converted_parameters

def base_representation(target, base, ndims):

    div = base**(ndims-1)
    res = target

    if target > base**ndims:
        raise ValueError('Input is greater than allowed output.')

    output = []
    for i in np.arange(ndims):
        mod = np.floor(res/div)
        res = res - (mod*div)
        div = div/base
        output.append(mod)

    return np.array(output).astype(int)

def normal_cdf(x):
    return 0.5*(1+erf(x/np.sqrt(2.)))

class Schechter(bilby.prior.Prior):
    def __init__(self, alpha, L0, minimum=0., maximum=np.inf, name=None,
                 latex_label=None, unit=None, boundary=None):        
        
        """

        https://en.wikipedia.org/wiki/Luminosity_function_(astronomy)

        Parameters
        ==========
        alpha: float
            Gradient parameter
        L0: float
            Break parameter
        minimum: float
            Minimum parameter
        maximum: float
            Maximum parameter
        name: str
            See superclass
        latex_label: str
            See superclass
        unit: str
            See superclass
        boundary: str
            See superclass
        """
        
        super(Schechter, self).__init__(name=name, latex_label=latex_label, unit=unit,
                                        minimum=minimum, maximum=maximum, boundary=boundary)
        
        self.alpha = alpha
        self.L0 = L0
        self.minimum = minimum
        self.maximum = maximum

    @property
    def normalisation(self):
        """ Calculates the proper normalisation of the Schechter function
        Returns
        =======
        float: Proper normalisation of the Schechter function
        """
        return self.L0 * gamma(self.alpha+1.) * (gammaincc(self.alpha+1.,self.minimum/self.L0) - gammaincc(self.alpha+1.,self.maximum/self.L0))

    def rescale(self, val):
        """
        'Rescale' a sample from the unit line element to the appropriate Schechter prior.
        This maps to the inverse CDF. This has been analytically solved for this case.
        """
        return self.L0*gammainccinv(self.alpha+1.,gammaincc(self.alpha+1.,self.minimum/self.L0)-(val*self.normalisation/(self.L0*gamma(self.alpha+1.))))

    def prob(self, val):
        """Return the prior probability of val.
        Parameters
        ==========
        val: Union[float, int, array_like]
        Returns
        =======
        float: Prior probability of val
        """
        return np.exp(self.ln_prob(val))#np.power(val/self.L0,self.alpha) * np.exp(-val/self.L0) * self.is_in_prior_range(val) / self.normalisation 

    def ln_prob(self, val):
        """
        Return the logarithmic prior probability of val
        Parameters
        ==========
        val: Union[float, int, array_like]
        Returns
        =======
        float:

        """
        with np.errstate(divide='ignore', invalid='ignore'):
            ln_in_range = np.log(1. * self.is_in_prior_range(val))
        return - np.log(self.normalisation) + (self.alpha*np.log(val/self.L0)) - (val/self.L0) + ln_in_range

class skewlognormal(bilby.prior.Prior):
    def __init__(self, mu, sigma, minimum=-np.inf, maximum=np.inf, name=None, latex_label=None, unit=None, boundary=None):
        """

        https://en.wikipedia.org/wiki/Luminosity_function_(astronomy)

        Parameters
        ==========
        mu: float
            Location parameter
        sigma: float
            Scale parameter
        minimum: float
            Minimum parameter
        maximum: float
            Maximum parameter
        name: str
            See superclass
        latex_label: str
            See superclass
        unit: str
            See superclass
        boundary: str
            See superclass
        """
        super(skewlognormal, self).__init__(name=name, latex_label=latex_label, minimum=minimum, maximum=maximum, unit=unit, boundary=boundary)

        self.mu = mu
        self.sigma = sigma

    def rescale(self, val):
        try: 
            L = len(val)
        except TypeError:
            L = 1 
        return np.zeros(L)

    def prob(self, val):
        """
        """
        eps = (val-self.mu)/self.sigma
        return np.exp(-0.5*eps**2)*np.exp(0.5*(eps+self.sigma)**2)*(1-normal_cdf(eps+self.sigma))

    def ln_prob(self, val):
        """
        """
        eps = (val-self.mu)/self.sigma
        return (-0.5*eps**2) + (0.5*(eps+self.sigma)**2) + np.log(1-normal_cdf(eps+self.sigma))

class GRBGWLikelihood(bilby.Likelihood):
    def __init__(self, x, y, function, beta=1., sigma=None):

        self.x = x
        self.JS = self.x[2]
        self.Nz_obs, self.f_obs = y
        self.rates = x[-2]
        self.GWtrig = x[-1]

        self.function = function
        self.sigma = sigma
        self.beta = beta

        parameters = ['theta_c', 'theta_j', 's', 'a', 'invLFc', 'logR_BNS', 'logR_NSBH', 'eps_NSBH', 'L0', 'alpha']
        self.parameters = dict.fromkeys(parameters)
        self.function_keys = self.parameters.keys()
        # These lines of code infer the parameters from the provided function
        #inspect.getargspec(function).args
        #self.parameters = dict.fromkeys(parameters[1:])
        #self.function_keys = self.parameters.keys()
        if self.sigma is None:
            self.parameters['sigma'] = None
        self._marginalized_parameters = []

    @property
    def marginalized_parameters(self):
        return self._marginalized_parameters

    def log_likelihood(self):
        sigma = self.parameters.get('sigma', self.sigma)
        L0 = np.power(10.,49+self.parameters['L0'])
        #1.e52*self.parameters['L0']
        LFc = 1./self.parameters['invLFc']
        R_BNS = np.power(10.,self.parameters['logR_BNS'])
        R_NSBH = np.power(10.,self.parameters['logR_NSBH'])

        if self.JS=='DG':
            a = np.exp(self.parameters['a'])
            s = np.exp(self.parameters['s'])
        else:
            a = self.parameters['a']
            s = self.parameters['s']

        Nz, mu = self.function(self.x, self.parameters['theta_c'], theta_j=self.parameters['theta_j'], s=s, a=a, LFc=LFc, R_BNS=R_BNS, R_NSBH=R_NSBH, eps_NSBH=self.parameters['eps_NSBH'], L0=L0, alpha=self.parameters['alpha'])
        
        logL = 0.
        if self.GWtrig:
            logL += logq1(self.f_obs, mu, self.sigma)
        if self.rates:
            logL += logq2(self.Nz_obs, Nz)
        logL = np.array(logL)

        return self.beta*logL

def condition_func_gamma(reference_paras, alpha):
    return dict(k=alpha)

def lognorm_rate(Rmin, Rmax, std_N=3./2):
    lnRmin = np.log(Rmin)
    lnRmax = np.log(Rmax)
    lnRmu = lnRmin + ((lnRmax-lnRmin)/2.)
    sig = (lnRmu - lnRmin)/std_N
    return np.exp(lnRmu), sig

def get_priors(JS, rate_width=1, eps=sys.float_info.epsilon, Lfixed=False, default=True):
    reparameterisations = {}
    priors = {}
    if JS=='GJ' or JS=='GB':
        priors['s'] = 0.
        priors['a'] = 0.
        priors['theta_c'] = bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_c')
        #bilby.core.prior.Gamma(2., np.pi/8., 'theta_c')
        priors['theta_j'] = np.pi/2.#bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_j')
        reparameterisations['theta_c'] = {'reparameterisation': 'angle'}
        #{'reparameterisation': 'scale', 'scale': np.pi/4.}
    elif JS=='PL':
        #priors = bilby.core.prior.PriorDict(conversion_function=convert_x_y_to_z)
        priors['s'] = bilby.core.prior.Gamma(2., 4., 's')
        priors['a'] = bilby.core.prior.Gamma(2., 1., 'a')
        priors['theta_c'] = bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_c')
        priors['theta_j'] = np.pi/2.#bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_j')
        #priors['z'] = bilby.core.prior.Constraint(minimum=eps, maximum=(np.pi/2.)-eps)
        #reparameterisations['theta_j'] = {'reparameterisation': 'angle'}
        reparameterisations['s'] = {'reparameterisation': 'scale', 'scale': 2.}
        reparameterisations['a'] = {'reparameterisation': 'scale', 'scale': 2.}
        reparameterisations['theta_c'] = {'reparameterisation': 'angle'}
    elif JS=='TH':
        priors['s'] = 0.
        priors['a'] = 0.
        priors['theta_j'] = bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_j')#bilby.core.prior.Gamma(2., np.pi/8., 'theta_j')
        priors['theta_c'] = np.pi/2.#bilby.core.prior.Uniform(1.e-9, np.pi/2., 'theta_j')
        reparameterisations['theta_j'] = {'reparameterisation': 'angle'}
    elif JS=='PZ':
        #priors = bilby.core.prior.PriorDict(conversion_function=convert_x_y_to_z)
        priors['s'] = bilby.core.prior.Gamma(2., 1., 's')#bilby.core.prior.Uniform(1.e-3, 20., 's')
        priors['a'] = bilby.core.prior.Gamma(2., 1., 'a')#bilby.core.prior.Uniform(1.e-3, 20., 'a')
        priors['theta_c'] = bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_c')
        priors['theta_j'] = np.pi/2.#bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_j')
        #priors['z'] = bilby.core.prior.Constraint(minimum=eps, maximum=(np.pi/2.)-eps)
        reparameterisations['theta_c'] = {'reparameterisation': 'angle'}
        #reparameterisations['theta_j'] = {'reparameterisation': 'angle'}
        reparameterisations['s'] = {'reparameterisation': 'scale', 'scale': 2.}
        reparameterisations['a'] = {'reparameterisation': 'scale', 'scale': 2.}
    elif JS=='DG':
        priors = bilby.core.prior.PriorDict(conversion_function=convert_x_y_to_z)
        priors['s'] = bilby.core.prior.Uniform(np.log(1.e-6), np.log(1.), 's')
        priors['a'] = bilby.core.prior.Uniform(np.log(1.e-6), np.log(1.), 'a')
        priors['theta_c'] = bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_c')
        priors['theta_j'] = bilby.core.prior.Uniform(eps, (np.pi/2.)-eps, 'theta_j')
        priors['z'] = bilby.core.prior.Constraint(minimum=eps, maximum=(np.pi/2.)-eps)
        reparameterisations['theta_c'] = {'reparameterisation': 'angle'}
        reparameterisations['theta_j'] = {'reparameterisation': 'angle'}
        reparameterisations['s'] = {'reparameterisation': 'to-cartesian'}
        reparameterisations['a'] = {'reparameterisation': 'to-cartesian'}
    if rate_width>0.:
        BNS_mu, BNS_sig = lognorm_rate((320.-240.)*1.e-9, (320.+490.)*1.e-9) #O3 rate
        NSBH_mu, NSBH_sig = lognorm_rate((45.-33.)*1.e-9, (130.+112.)*1.e-9) #O3 rate
        
        priors['logR_BNS'] = bilby.core.prior.Normal(np.log10(BNS_mu), BNS_sig, 'logR_BNS')
        priors['logR_NSBH'] = np.log10(NSBH_mu)#bilby.core.prior.Normal(np.log(NSBH_mu), NSBH_sig, 'logR_NSBH')
        
        reparameterisations['logR_BNS'] = {'reparameterisation': 'scale', 'scale': -np.log(BNS_mu)}
        #reparameterisations['logR_NSBH'] = {'reparameterisation': 'scale', 'scale': -np.log(NSBH_mu)}
    else:
        priors['R_BNS'] =  BNS_mu#np.log(1.e-6)
        priors['R_NSBH'] = NSBH_mu#np.log(1.e-7)

    priors['eps_NSBH'] = 0.#bilby.core.prior.Uniform(0., 1., 'eps_NSBH')
    #reparameterisations['eps_NSBH']={'reparameterisation': 'to-cartesian'}

    if default:
        priors['invLFc'] = np.log(400.)
    else:    
        priors['invLFc'] = bilby.core.prior.Gamma(2., 1./400, 'invLFc')#bilby.core.prior.Uniform(np.log(5.), np.log(10000.), 'LFc')
        reparameterisations['invLFc']={'reparameterisation': 'scale', 'scale': 0.25/100}#{'reparameterisation': 'to-cartesian'}

    if Lfixed:
        priors['alpha'] = 0.55
        priors['L0'] = 2.6

    else:
        priors['alpha'] = bilby.core.prior.Uniform(1.e-9, 1., 'alpha')
        priors['L0'] = bilby.core.prior.Normal(2.6, 1., 'L0')

        reparameterisations['alpha']={'reparameterisation': 'to-cartesian'}
        reparameterisations['L0']={'reparameterisation': 'scale', 'scale': 2.}#{'reparameterisation': 'to-cartesian'}

    return priors, reparameterisations

def JS_model_interp_nessai(JS, thetavs, Fs, rs, sigmas, outdir, label, Nz_obs=False, eps=sys.float_info.epsilon, beta=1., draws=2000, burn=1000, progressbar=False, random_seed=np.random.randint(1.e8), debug=False, load=False, sampler='nessai', interp=True, cores=None, n_pool=None, Lfixed=False, default=True):

    label=label.replace("/","_")
    bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)

    T_0, f_r, f_FOV, eps_BNS, R_GW = cs.T_0, cs.f_r, cs.f_FOV, cs.eps_BNS, cs.R_GW
    rate_paras = T_0*f_FOV
    background = 0.

    JSfunc = beam_function(JS)
    LFfunc = LF_function(JS)

    if interp:
        ranges, Z = fetch_interp(JS)
        rate_interp = interpolate.RegularGridInterpolator(ranges, Z, bounds_error=False, fill_value=0.)
    else:
        rate_interp = None

    if Nz_obs:
        rates_ = True
        Nz_obs = np.array([int(Nz_obs)])
    else:
        rates_ = False
        Nz_obs = False

    if np.any(Fs):
        GWtrigs_ = True
        zs = (rs/3.086e24)*cs.H_0/cs.c
        k_cor = rtfn.k_cor(zs, e_min=50., e_max=300.)
        f_obs = Fs
    else:
        GWtrigs_ = False
        zs = False
        k_cor = False
        f_obs = False

    Y = [Nz_obs, f_obs]

    likelihood = GRBGWLikelihood((rs, thetavs, JS, rate_paras, background, k_cor, zs, JSfunc, LFfunc, rate_interp, rates_, GWtrigs_), Y, joint_func_interp, beta=1, sigma=sigmas)

    rate_width = 1.
    priors, reparameterisations = get_priors(JS, rate_width=rate_width, Lfixed=Lfixed, default=default)

    nlive=2000
    maximum_uninformed=4000
    update_poolsize=True
    poolsize=2000
    flow_config={"batch_size": 1000, "max_epochs": 500, "patience": 20,
                    "model_config": {"n_neurons": 8, "n_blocks": 4, "n_layers": 2,
                                    "kwargs": {"batch_norm_between_layers": True,
                                                "linear_transform": 'lu'}
                                    }
                }

    if sampler =='nessai':
        result = bilby.run_sampler(likelihood=likelihood, priors=priors, sampler='nessai', nlive=nlive, poolsize=poolsize,
                               outdir=outdir, label=label, plot=True, maximum_uninformed=maximum_uninformed,
                               flow_config=flow_config, analytic_priors=False, resume=True, update_poolsize=update_poolsize,
                               reparameterisations=reparameterisations, max_threads=cores, n_pool=n_pool, allow_multi_valued_likelihood=True)
    elif sampler=='dynesty':
        result = bilby.run_sampler(likelihood=likelihood, priors=priors, sampler='dynesty', nlive=draws,
                                   outdir=outdir, label=label, plot=True)

    return result

def fetch_interp(JS):
    if JS=='GJ' or JS=='GB':
        interppath = '/home/fergus.hayes/public_html/LVC/AutoDoc/2021_04_15/SIMULATION_theta_c_LF_GJ_largest/JS_GB_func_rates_int_Nthv_0_Nthj_0_Nthc_2_Ns_0_Na_0_NLF_2/13_45_11/'
    elif JS=='TH':
        interppath = '/home/fergus.hayes/public_html/LVC/AutoDoc/2021_05_14/SIMULATION_interpolate_TH_1000_300x300/JS_TH_func_rates_int_Nthv_0_Nthj_300_Nthc_0_Ns_0_Na_0_NLF_300/12_39_34/'
    elif JS=='DG':
        interppath = '/home/fergus.hayes/public_html/LVC/AutoDoc/2021_05_19/SIMULATION_interpolate_DG_1000ish_10x10x10x10x10/JS_DG_func_rates_int_Nthv_0_Nthj_10_Nthc_10_Ns_10_Na_10_NLF_10/09_59_26/'

    rundirs = np.sort(np.array(next(os.walk(interppath))[1]).astype(int)).astype(str)

    if JS=='GJ' or JS=='GB' or JS=='TH':
        Z = []
        for rundir in rundirs:
            interpfiles = next(os.walk(interppath+rundir+'/'))[2]
            for interpfile in interpfiles:
                if interpfile.endswith('ranges.npz'):
                    keys = np.array(np.load(interppath+rundir+'/'+interpfile, allow_pickle=True)['arr_0']).astype(str)
                    ranges = np.load(interppath+rundir+'/'+interpfile, allow_pickle=True)['arr_1']
                elif interpfile.endswith('Z.npz'):
                    Z.append(np.load(interppath+rundir+'/'+interpfile, allow_pickle=True)['arr_0'])

        Z = np.concatenate(Z).astype(float)
        Z = Z.reshape([range_.shape[0] for range_ in ranges])

    elif JS=='DG':
        Ns = np.array(interppath.split('/')[-3].split('_')[-11:][::2]).astype(int)
        varlist = np.array(interppath.split('/')[-3].split('_')[-12:][::2])[Ns!=0]
        Ns = Ns[Ns!=0]

        Z = np.zeros((int(np.power(len(rundirs),1./len(Ns)))*Ns))

        for rundir in rundirs:
            interpfiles = next(os.walk(interppath+rundir+'/'))[2]
            for interpfile in interpfiles:
                if interpfile.endswith('ranges.npz'):
                    keys = np.array(np.load(interppath+rundir+'/'+interpfile, allow_pickle=True)['arr_0']).astype(str)
                    ranges = np.load(interppath+rundir+'/'+interpfile, allow_pickle=True)['arr_1']
                    foundZ = True
                elif interpfile.endswith('Z.npz'):
                    G=np.load(interppath+rundir+'/'+interpfile, allow_pickle=True)['arr_0']
                    #print(G.shape)
                    coord = base_representation(int(rundir), int(np.power(len(rundirs),1./len(Z.shape))), len(Z.shape))
                    #print(coord)
                    if len(Ns)==1:
                        Z[coord[0]*Ns[0]:(coord[0]+1)*Ns[0]] = G.reshape((Ns))
                    elif len(Ns)==2:
                        Z[coord[0]*Ns[0]:(coord[0]+1)*Ns[0],coord[1]*Ns[1]:(coord[1]+1)*Ns[1]] = G.reshape((Ns))
                    elif len(Ns)==3:
                        Z[coord[0]*Ns[0]:(coord[0]+1)*Ns[0],coord[1]*Ns[1]:(coord[1]+1)*Ns[1],coord[2]*Ns[2]:(coord[2]+1)*Ns[2]] = G.reshape((Ns))
                    elif len(Ns)==4:
                        Z[coord[0]*Ns[0]:(coord[0]+1)*Ns[0],coord[1]*Ns[1]:(coord[1]+1)*Ns[1],coord[2]*Ns[2]:(coord[2]+1)*Ns[2],coord[3]*Ns[3]:(coord[3]+1)*Ns[3]] = G.reshape((Ns))
                    elif len(Ns)==5:
                        Z[coord[0]*Ns[0]:(coord[0]+1)*Ns[0],coord[1]*Ns[1]:(coord[1]+1)*Ns[1],coord[2]*Ns[2]:(coord[2]+1)*Ns[2],coord[3]*Ns[3]:(coord[3]+1)*Ns[3],coord[4]*Ns[4]:(coord[4]+1)*Ns[4]] = G.reshape((Ns))
                    elif len(Ns)==6:
                        Z[coord[0]*Ns[0]:(coord[0]+1)*Ns[0],coord[1]*Ns[1]:(coord[1]+1)*Ns[1],coord[2]*Ns[2]:(coord[2]+1)*Ns[2],coord[3]*Ns[3]:(coord[3]+1)*Ns[3],coord[4]*Ns[4]:(coord[4]+1)*Ns[4],coord[5]*Ns[5]:(coord[5]+1)*Ns[5]] = G.reshape((Ns))

        Z = Z.flatten()
    
        dims = []
        for para_range in ranges:
            dims.append(len(para_range))
        dims = np.array(dims)
        Z = Z.reshape(dims)

    return ranges, Z
