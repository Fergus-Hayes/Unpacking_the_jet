import numpy as np
import matplotlib.pyplot as plt
import datetime, argparse, json, matplotlib
from intervals import credible_interval

def main(locs, labels, outpath, N=2000):

    if len(locs)!=len(labels):
        raise(ValueError)

    fontsize = 28
    ticksize = 22
    figsize = (8,12)

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    BNSrates = []
    cmins = []
    cmaxs = []

    CI = 0.9#5

    colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', 'black', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']


    for i,loc in enumerate(locs):
        with open(loc) as f:
            sample_file = json.load(f)
        BNSrates.append(np.random.choice(np.array(sample_file['posterior']['content']['logR_BNS']).flatten(), N, replace=False))
        ci_min,ci_max=credible_interval(BNSrates[-1], alpha=CI, options='HPD', N=1000)
        cmins.append(ci_min)
        cmaxs.append(ci_max)

        print(np.round(np.median(BNSrates[-1]),2), np.round(cmaxs[-1]-np.median(BNSrates[-1]),2), np.round(np.median(BNSrates[-1])-cmins[-1],2), loc)
        #print(np.round(1.e9*np.power(10.,cmins[-1]),2), np.round(1.e9*np.power(10.,cmaxs[-1]),2))
        print(np.median(np.power(10,BNSrates[-1])*1.e9),credible_interval(np.power(10,BNSrates[-1])*1.e9, alpha=CI, options='HPD', N=1000),labels[i])

    BNSrates.append(np.random.normal(-6.6,0.7,2*N))
    
    ci_min,ci_max=credible_interval(BNSrates[-1], alpha=CI, options='HPD', N=1000)
    cmins.append(ci_min)
    cmaxs.append(ci_max)

    labels.append('Prior')

    print(np.round(np.median(BNSrates[-1]),2), np.round(cmaxs[-1]-np.median(BNSrates[-1]),2), np.round(np.median(BNSrates[-1])-cmins[-1],2), loc)
    print(np.median(np.power(10,BNSrates[-1])*1.e9),credible_interval(np.power(10,BNSrates[-1])*1.e9, alpha=CI, options='HPD', N=1000),labels[-1])

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    pos = np.arange(len(BNSrates))*2
    parts = ax.violinplot(BNSrates, positions=pos, showmeans=True, widths=1., vert=False)

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
    ax.set_xlabel(r'$\log_{10}(R_{\regular{BNS}}/\regular{Mpc}^{-3}\,\regular{yr}^{-1})$', fontsize=fontsize)
    ax.set_ylabel(r'', fontsize=fontsize)
    
    ax.set_yticks(pos)
    ax.set_yticklabels(labels)
    
    ax.set_xticks(np.round(np.linspace(np.min(BNSrates[-1])-0.2,np.max(BNSrates[-1])+0.2,7),1))

    ax.legend(fontsize=ticksize, framealpha=0.)

    fig.savefig(outpath+'BNS_rate.png', bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Use Bibly to estimate jet structure parameters from GRB and GW NSNS rates.")
    parser.add_argument('--loc', help="Result files", type=str, nargs='+', required=True)
    parser.add_argument('--labels', help="Labels", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="Labels", type=str, required=False, default='/home/fergus/Documents/Projects/GRB/Rates/new_results/plots/')

    opt = parser.parse_args()

    main(opt.loc, opt.labels, opt.outpath)

