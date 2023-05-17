"""
Model Comparison
=============================

This file initialises the sampling of the different models.

Copyright 2018 Fergus Hayes

"""

import matplotlib
matplotlib.use("agg")

import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from models.models import JS_model_interp_nessai
import argparse, os, time

np.random.seed(int(time.time()))

def main(path, structure, GW170817, GW190425, rates, beta=1, burn=2000, draws=1000, progressbar=False, interp=False, cores=None, Lfixed=False, default=False, random_seed=None, datapath='./inputdata/'):

    ######### Sample model #########

    Y = []
    F_sig = []
    d = []
    cosi = []
    Nz = False

    label=''

    if rates:
        label+='rates_'
        Nz = np.load(datapath+'rates.npy')
    if GW170817:
        label+='GW170817_'
        Y.append(np.load(datapath+'GW170817_Flux.npy').squeeze())
        F_sig.append(np.load(datapath+'GW170817_Flux_uncertainty.npy').squeeze())
        d_, cosi_ = np.load(datapath+'GW170817_dcosi_samples.npy')
        d.append(d_)
        cosi.append(cosi_)
    if GW190425:
        label+='GW190425_'
        Y.append(np.load(datapath+'GW190425_Flux.npy').squeeze())
        F_sig.append(np.load(datapath+'GW190425_Flux_uncertainty.npy').squeeze())
        d_, cosi_ = np.load(datapath+'GW190425_dcosi_samples.npy')
        d.append(d_)
        cosi.append(cosi_)

    if label=='':
        raise ValueError('Please specify which data to analyze with at least one of the flags: --GW170817 --GW190425 --rates')

    label+=structure

    if len(Y)==0:
        Y, F_sig, d, thetav = False, False, False, False
    else:
        Y, F_sig, d, thetav = np.array(Y), np.array(F_sig), np.array(d), np.arccos(cosi)

    if isinstance(cores,str):
        cores = int(cores)

    if cores==None:
        n_pool=None
    elif cores==1:
        n_pool=None
    elif cores>1:
        n_pool=cores-1

    if not structure:
        structure = str(path.split('/')[-3].split('_')[0])

    if Lfixed:
        Lstr='schechter_luminosity'
    else:
        Lstr='hierarchical_luminosity'

    path = path+'/'+Lstr+'/'+label+'/temp/'+str(structure)+'/'+str(int(np.round(beta)))+'/'

    if not os.path.exists(path):
        os.makedirs(path)

    # load or create random seed
    if os.path.exists(path+'random_seed.npy'):
        random_seed = int(np.load(path+'random_seed.npy'))
    elif type(random_seed)!=int:
        random_seed = np.random.randint(1.e8)
        np.save(path+'random_seed',random_seed)

    # Run sampler on given jet structure models. The models are specified in "./models/models.py"
    
    result = JS_model_interp_nessai(structure, thetav, Y, d, F_sig, Nz_obs=Nz, beta=beta, outdir=path+'trace/', label=label, burn=burn, draws=draws, progressbar=False, random_seed=random_seed, interp=interp, sampler='nessai', cores=cores, n_pool=n_pool, Lfixed=Lfixed, default=default)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Runs MCMC given structure model.")
    parser.add_argument('--path', help="Path to run.")
    parser.add_argument('--label', help="Run label.")
    parser.add_argument('--datapath', help="Path to data.", default=False)

    parser.add_argument('--GW170817', help="Add GW170817 data.", action='store_true', default=False)
    parser.add_argument('--GW190425', help="Add GW190425 data.", action='store_true', default=False)
    parser.add_argument('--rates', help="Add rates data.", action='store_true', default=False)

    parser.add_argument('--beta', help="Temperature of sampling.", default='1')
    parser.add_argument('--struct', help="Structure to analyse.", default=False)
    parser.add_argument('--burn', help="Number of samples to discard.", default=500000)
    parser.add_argument('--draws', help="Number of samples to keep.", default=100000)
    parser.add_argument('--progressbar', help="See progressbar.", default=False)
    parser.add_argument('--interp', help="Interpolation.", action='store_true', default=False)
    parser.add_argument('--Lfixed', help="Fix luminosity function.", action='store_true', default=False)
    parser.add_argument('--cores', help="Number of cores to use.", default=None)
    parser.add_argument('--default', help="Default options (DEPRICATED).", action='store_true', default=False)

    parser.add_argument('--randomseed', help="Random seed.", default=None)

    opt = parser.parse_args()

    main(opt.path, opt.struct, opt.GW170817, opt.GW190425, opt.rates, beta=float(opt.beta), burn=int(opt.burn), draws=int(opt.draws), progressbar=bool(opt.progressbar), interp=opt.interp, cores=opt.cores, Lfixed=opt.Lfixed, default=opt.default, random_seed=opt.randomseed, datapath='./inputdata/')
