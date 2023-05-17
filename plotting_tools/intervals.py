#!/usr/bin/python

import argparse
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt

def credible_interval(samples, alpha=0.9, options='HPD', N=1000):
    """
    Find the Bayesian credible interval from histograming given a sample from a posterior.
    
    Parameters:
    samples: 1D array of samples from the posterior.
    alpha (default=0.9): Confidence interval.
    options (default='HPD'): Options for the credible interval calculation.
                             - 'HPD': Narrowest range to give CI alpha
    N (default=1000): Histogram binning to calculate the CI over. 
    """
    if options=='HPD':
        x_min=np.min(samples)
        x_max=np.max(samples)
        counts,bins=np.histogram(samples,bins=N)
        counts_T=len(samples)
        for i,ibin in enumerate(bins):
            x_range=x_max-x_min
            for j_,jbin in enumerate(bins[i:]):
                j=i+j_
                range=bins[j]-bins[i]
                if range>=x_range:
                    break
                if (np.sum(counts[i:j])/counts_T)>=alpha:
                    x_min=bins[i]
                    x_max=bins[j]
                    break
    return x_min, x_max                


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', required=True, type=str, help='numpy file containing samples')
    parser.add_argument('--para', required=True, type=str, help='Name of parameter considered')
    parser.add_argument('--alpha', type=float, default=0.9, help='Credible interval')
    parser.add_argument('--options', type=str, default='HPD',
                        help='Type of credible interval')
    parser.add_argument('--N', type=int, default=1000,
                        help='Number of intervals to consider')
    opt = parser.parse_args()

    #samples=np.load(opt.samples)[opt.para]
    samples = np.random.beta(a=2,b=5,size=10000)

    print(credible_interval(samples,opt.alpha,opt.options,opt.N))
