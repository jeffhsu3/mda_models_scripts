import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors
import numpy as np

def plot_allele_frequency(SFS):
    """Plots site frequency spectrum and joint SFS

    Parameters
    ---------
    SFS : obj
        from figs_popgen.site_freqspec_fx
    jSFS : obj
        from figs_popgen.site_freqspec_fx
    Returns
    -------
    SFS : fig
        plot of SFS with each a different color on same plot
    jSFS : fig
        joint SFS between 2 villages
    """

    return(None)


def plot_allele_trace(allele_trace):
    """Plots trace of fitness influencing alleles

    Parameters
    ----------
    allele_trace : obj
        obj created from figs_popgen.sel_trace_fx

    Returns
    -------
    allele_trace : fig

    """
    #dfworm.sel['1S'] > 0
    #dfworm.sel['2S'] > 0
    #dfworm.sel['1F'] > 0
    #np.where(dfworm.sel['2F'] > 0)
    #fig.savefig(outfile)


def plot_coordinates_host(dfHost, thetaHost, height=7, width=7):
    """Plots host coordinates with size by theta

    Parameters
    ----------
    dfHost : df
        df with host data, needs to use coordinates
    thetaHost : dict w/ numpy array
        ['hostidx'], ['theta'] theta values by hostidx for MF stage
    Returns
    -------
    host_plot : fig

    """
    #TO DO: Add bubble size for theta value
    cmap = colors.ListedColormap(['k','g','y','r','b'])
    bounds=[0,1,2,3,4]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots()
    ax.scatter([i[0] for i in dfHost.coordinates],
            [i[1] for i in dfHost.coordinates],
            c = dfHost.village, cmap=cmap)
    ax.set_ylabel('y position')
    ax.set_xlabel('x position')
    fig.savefig('test.png')
    plt.tight_layout()
    return(fig)


def plot_hapnetwork(hapnetobj):
    '''Plots haplotype network

    Parameters
    ----------
    hapnetwork : df
        hapnetwork values from MF from figs_popgen.haplotype_net_fx,

    Returns
    -------
    hap_network : fig
    '''
    return(None)


def plot_pairwise(pairwise):
    '''Plots pairwise matrix as surface

    Parameters
    ----------
    pairwise : matrix
        matrix of values to plot for each host from figs_popgen

    Returns
    -------
    pairwise_stat : fig
    '''
    return(None)
