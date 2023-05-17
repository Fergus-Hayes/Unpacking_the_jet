import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import corner, json
import datetime, argparse, os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as mpimg

def overlaid_corner(samples_list_, sample_labels_, labels, truths=None, label_kwargs={"fontsize": 16}, show_titles=True, JS=False):
    """Plots multiple corners on top of each other"""

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    # get some constants
    cmap = plt.cm.gist_rainbow
    colors_ = np.array(['black','blue','violet', 'orange', 'red'])[::-1]
    
    samples_list = []
    sample_labels = []
    colors = []

    for i,samples_ in enumerate(samples_list_):
        if np.all(np.logical_not(samples_==None)):
            samples_list.append(samples_)
            sample_labels.append(sample_labels_[i])
            colors.append(colors_[i])

    colors = np.array(colors)

    _, ndim = samples_list[0].shape
    max_len = max([len(s) for s in samples_list])
    n = len(samples_list)

    plot_range = []
    for dim in range(ndim):
        plot_range.append(
            [
                min([min(samples_list[i].T[dim]) for i in range(n)]),
                max([max(samples_list[i].T[dim]) for i in range(n)]),
            ]
        )

    fig = corner.corner(
        samples_list[0],
        color=colors[0],
        weights=get_normalisation_weight(len(samples_list[0]), max_len),
        labels=labels,
        truths=truths,
        show_titles=show_titles,
        smooth=0.9,
        range=plot_range,
        label_kwargs=dict(fontsize=5+int(ndim*3)),
        title_kwargs=dict(fontsize=3+int(ndim*2))
        )

    for idx in range(1, n):
        fig = corner.corner(
            samples_list[idx],
            fig=fig,
            weights=get_normalisation_weight(len(samples_list[idx]), max_len),
            color=colors[idx],
            labels=labels,
            truths=truths,
            show_titles=show_titles,
            smooth=0.9,
            range=plot_range,
            label_kwargs=dict(fontsize=5+int(ndim*3)),
            title_kwargs=dict(fontsize=3+int(ndim*2))
            )

    plt.legend(
        handles=[
            mlines.Line2D([], [], color=colors[i], label=sample_labels[i])
            for i in range(n)
        ],
        fontsize=5+int(ndim*3), frameon=False,
        bbox_to_anchor=(1, ndim),
    )
    
    for ax in fig.get_axes():
        ax.tick_params(axis='both', labelsize=3+int(ndim*2))

    return fig


def get_normalisation_weight(len_current_samples, len_of_longest_samples):
    return np.ones(len_current_samples) * (len_of_longest_samples / len_current_samples)

def get_samples(loc, JS=False):

    #L0 = 10.**51.5

    injlabels = np.array(['theta_c', 'log R_BNS', 'log R_NSBH', 'eps_NSBH', 'LFc', 'E0'])

    try:
        with open(loc) as f:
            sample_file = json.load(f)
    except OSError:
        return None, None

    samples_ = sample_file['posterior']['content']
    labels = np.array(sample_file['parameter_labels'])
    labels_ = np.array(list(samples_.keys()))

    theta_labels = []
    ntheta_labels = []
    for label in labels_:
        if label[:5]=='theta':
            theta_labels.append(label)
        else:
            ntheta_labels.append(label)
    labels = np.concatenate([theta_labels,ntheta_labels])

    print(loc,np.round(sample_file['log_evidence'],2), np.round(sample_file['log_evidence_err'],2))

    if labels[-1]==None:
        labels[-1] = 'E0'

    plot_labels = []
    samples = []
    
    for i,label in enumerate(labels):
        if label=='log_prior' or label=='log_likelihood':
            continue
        if label[:5]=='theta':
            samples.append(180.*np.array(samples_[label])/np.pi)
            plot_labels.append(r'$\theta_{'+label.split('_')[1]+'}$ '+u'(\u00b0)')
        if label[:2]=='R_':
            samples.append(np.log10(np.array(samples_[label])))
            plot_labels.append('$R_{'+label.split('_')[1]+'}$')# (Gpc$^{-3}\,$y$^{-1}$)')
        if label[:4]=='eps_':
            plot_labels.append('$\epsilon_{'+label.split('_')[1]+'}$')
            samples.append(np.array(samples_[label]))
        if label[:3]=='log':
            samples.append(np.array(samples_[label]))#np.log10(np.exp(np.array(samples_[label]))))
            plot_labels.append(r'log$_{10}$ $R_{'+label.split('_')[1]+'}$')# (Gpc$^{-3}\,$y$^{-1}$)')
        if label=='LFc':
            plot_labels.append(r'$\log_{10}$ $\Gamma_{0}$')
            samples.append(np.log10(np.array(samples_[label])))
        if label=='invLFc':
            plot_labels.append(r'$\log_{10}$ $\Gamma_{0}$')
            samples.append(np.log10(1./np.array(samples_[label])))
        if label=='a':
            if np.any(np.array(labels)=='s')*np.any(np.array(labels)=='a'):
                if JS=='DG':
                    samples_[label] = np.exp(samples_[label])
                samples.append(np.log10(samples_[label]))
                plot_labels.append(r'log$_{10}$ $\mathcal{A}$')
            else:
                samples.append(samples_[label])
                plot_labels.append(r'$a$')
        if label=='s':
            if np.any(np.array(labels)=='s')*np.any(np.array(labels)=='a'):
                if JS=='DG':
                    samples_[label] = np.exp(samples_[label])
                samples.append(np.log10(samples_[label]))
                plot_labels.append(r'log$_{10}$ $\mathcal{C}$')
            else:
                samples.append(samples_[label])
                plot_labels.append(r'$s$')
        if label=='L0':
            samples.append(np.array(samples_[label])+49.)
            plot_labels.append(r'$\log_{10}$ $L_{0}^{\ast}$')

        if label=='alpha':
            samples.append(np.array(samples_[label]))
            plot_labels.append(r'$\gamma$')

    samples = np.array(samples).T

    plot_labels = np.array(plot_labels)

    if np.any(plot_labels==r'$\theta_{c}$ '+u'(\u00b0)'):
        samples = samples[np.array(samples[:,plot_labels==r'$\theta_{c}$ '+u'(\u00b0)']<30.).squeeze()]
    if np.any(plot_labels==r'$\theta_{j}$ '+u'(\u00b0)'):
        samples = samples[np.array(samples[:,plot_labels==r'$\theta_{j}$ '+u'(\u00b0)']<90.).squeeze()]
    if np.any(plot_labels==r'log$_{10}$ $\Gamma_{0}$'):
        samples = samples[np.array(samples[:,plot_labels==r'log$_{10}$ $\Gamma_{0}$']<np.log10(1000.)).squeeze()]

    if JS:
        if JS=='GJ':
            plot_labels[plot_labels==r'$\theta_{c}$ '+u'(\u00b0)'] = r'$\theta_{\sigma}$ '+u'(\u00b0)'
        if JS=='DG':
            plot_labels[plot_labels==r'$\theta_{c}$ '+u'(\u00b0)'] = r'$\theta_{in}$ '+u'(\u00b0)'
            plot_labels[plot_labels==r'$\theta_{j}$ '+u'(\u00b0)'] = r'$\theta_{out}$ '+u'(\u00b0)'
        if JS=='PL':
            plot_labels[plot_labels==r'log$_{10}$ $\mathcal{A}$'] = r'$a$'
            samples[:,plot_labels==r'$a$'] = np.power(10.,samples[:,plot_labels==r'$a$'])
            plot_labels[plot_labels==r'log$_{10}$ $\mathcal{C}$'] = r'$s$'
            samples[:,plot_labels==r'$s$'] = np.power(10.,samples[:,plot_labels==r'$s$'])
        if JS=='PZ':
            Ss = np.power(10.,samples[:,plot_labels==r'log$_{10}$ $\mathcal{A}$'])
            As = np.power(10.,samples[:,plot_labels==r'log$_{10}$ $\mathcal{C}$'])
            plot_labels[plot_labels==r'log$_{10}$ $\mathcal{A}$'] = r'$a$'
            plot_labels[plot_labels==r'log$_{10}$ $\mathcal{C}$'] = r'$s$'
            samples[:,plot_labels==r'$s$'] = Ss
            samples[:,plot_labels==r'$a$'] = As

    return samples, plot_labels

def main(locs, data_labels_, outpath, injectfile=False, JS=False):

    data_labels_ = [r'$\mathcal{D}_{170817}$', r'$\mathcal{D}_{190425}$', r'$\mathcal{D}_{R}$', r'$\mathcal{D}_{170817+R}$', r'$\mathcal{D}_{\regular{all}}$']

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    data_labels = data_labels_

    injections_ = None
    injections=[]

    if injectfile:
        injections_ = np.load(injectfile)
        print(injections_)
        injections_[0] = np.rad2deg(injections_[0])
        injections_[1:3] = np.log10(np.exp(injections_[1:3]))
        injections_[-1] = np.log10(injections_[-1])
    injections = injections_

    data_samples = []
    for loc in locs:
        samples_, plot_labels_ = get_samples(loc, JS=JS)
        data_samples.append(samples_)
        if not np.any(plot_labels_==None):
            plot_labels = plot_labels_    


    fig = overlaid_corner(data_samples, data_labels, labels=plot_labels, truths=injections, label_kwargs={"fontsize": 16}, show_titles=True, JS=JS)
    fig.savefig(outpath+'joint_corner.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Use Bibly to estimate jet structure parameters from GRB and GW NSNS rates.")
    parser.add_argument('--loc', help="Path where files are saved", type=str, nargs='+', required=True)
    parser.add_argument('--labels', help="List of dataset names", type=str, nargs='+', required=True)
    parser.add_argument('--parafile', help="Injection files", default=False, type=str)
    parser.add_argument('--JS', help="Jet structure", default=False, type=str)
    parser.add_argument('--outpath', help="Location to save plot", required=True, type=str)

    opt = parser.parse_args()

    main(opt.loc, opt.labels, opt.outpath, opt.parafile, opt.JS)

