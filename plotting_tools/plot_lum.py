import numpy as np
import matplotlib.pyplot as plt
import datetime, argparse, json, matplotlib
from scipy.special import gammaln, gamma, gammainc, gammaincc, gammaincinv, gammainccinv, logsumexp, erf
from intervals import credible_interval

def trunc_gamma(shape, scale=1.0, a=0., b=np.inf, size=None):
   return scale*gammainccinv(shape,(gammaincc(shape,a/scale)-(np.random.uniform(0,1,size=size)*(gammaincc(shape,a/scale)-gammaincc(shape,b/scale)))))


def main(locs, labels, outpath, N=2000):

    if len(locs)!=len(labels):
        raise(ValueError)

    fontsize = 28
    ticksize = 22
    figsize = (8,12)

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    L0s = []
    cmins = []
    cmaxs = []

    CI = 0.9#5

    for loc in locs:
        with open(loc) as f:
            sample_file = json.load(f)
        L0 = np.random.choice(np.power(10.,np.array(sample_file['posterior']['content']['L0']).flatten()+49.), N, replace=False)
        Lmin = L0*1.e-3
        Lmax = L0*1.e2
        alpha = np.random.choice(np.array(sample_file['posterior']['content']['alpha']).flatten(), N, replace=False)
        
        L0_mu = []
        for i in np.arange(len(L0)):
            L0_mu.append(np.log10(np.median(trunc_gamma(alpha[i], L0[i], a=Lmin[i], b=Lmax[i], size=1000))))
        L0s.append(np.array(L0_mu))
        
        ci_min,ci_max=credible_interval(L0s[-1], alpha=CI, options='HPD', N=1000)
        cmins.append(ci_min)
        cmaxs.append(ci_max)

        print(np.round(np.median(L0_mu),1), np.round(np.median(L0_mu)-ci_min,1), np.round(ci_max-np.median(L0_mu),1), loc)

    labels.append('Prior')

    print('Prior')

    L0 = np.power(10.,np.random.normal(49.+2.6,1,N))#np.power(10.,np.random.normal(51,1,N))#np.random.uniform(49,53,N))
    alpha = np.random.uniform(1.e-9,1.,N)#np.random.uniform(2,3,N)
    Lmin = L0*1.e-3
    Lmax = L0*1.e2
    L0_mu = []
    for i in np.arange(len(L0)):
        L0_mu.append(np.log10(np.median(trunc_gamma(alpha[i], L0[i], a=Lmin[i], b=Lmax[i], size=1000))))
    L0s.append(np.array(L0_mu))
    
    #L0_mu = np.log10(((1+alpha)*L0))#np.log10(L0/(alpha-2))

    ci_min,ci_max=credible_interval(L0s[-1], alpha=CI, options='HPD', N=1000)
    cmins.append(ci_min)
    cmaxs.append(ci_max)

    print(np.round(np.median(L0_mu),1), np.round(np.median(L0_mu)-ci_min,1), np.round(ci_max-np.median(L0_mu),1))

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    pos = np.arange(len(L0s))*2
    parts = ax.violinplot(L0s, positions=pos, showmeans=True, widths=1., vert=False)

    colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', 'black', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']

    color='black'

    for i,pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor(colors[i])
    for partname in ('cbars','cmins','cmaxes','cmeans'):
        vp = parts[partname]
        vp.set_edgecolor(color)
        vp.set_linewidth(1)

    dp = 0.2
    for cmin, cmax, p in zip(cmins,cmaxs,pos):
        ax.vlines(cmin,(p-dp),(p+dp),color=color,ls='--',lw=2)
        ax.vlines(cmax,(p-dp),(p+dp),color=color,ls='--',lw=2)


    ax.tick_params(axis='both', labelsize=ticksize)
    ax.set_xlabel(r'$\log_{10}(\hat{L}_{0}/\regular{erg}\,\regular{s}^{-1})$', fontsize=fontsize)
    ax.set_ylabel(r'', fontsize=fontsize)
    
    ax.set_yticks(pos)
    ax.set_yticklabels(labels)
    
    ax.set_xticks(np.round(np.linspace(np.min(L0s),np.max(L0s),7),1))

    ax.legend(fontsize=ticksize, framealpha=0.)

    fig.savefig(outpath+'L0s.png', bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Use Bibly to estimate jet structure parameters from GRB and GW NSNS rates.")
    parser.add_argument('--loc', help="Result files", type=str, nargs='+', required=True)
    parser.add_argument('--labels', help="Labels", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="Labels", type=str, required=False, default='/home/fergus/Documents/Projects/GRB/Rates/new_results/plots/')

    opt = parser.parse_args()

    main(opt.loc, opt.labels, opt.outpath)

