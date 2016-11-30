import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.constants import c
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.isgm.abscomponent import AbsComponent
from linetools.analysis import voigt as lav
from linetools import utils as ltu
from pyntejos import utils as pntu

# color style for spectra
COLOR_FLUX = 'k'
COLOR_SIG = 'g'
COLOR_BLEND = '#CCCCCC'
COLOR_RESIDUAL = 'grey'
COLOR_MODEL = 'r'


def common_labels(fig, xlabel='',ylabel='',title='',fontsize=20,xlabelpad=None,ylabelpad=None):
    """Plots only labels in a parent figure. Special for subplots with
    shared axes.
    
    Mode of use:

    fig  = plt.figure()
    common_labels(fig, xlabel='x label',ylabel='y label',fontsize=20)
    
    ax1 = fig.add_subplot(1,2,1)
    plt.plot(bla bla)
    ax2 = fig.add_subplot(1,2,2)
    plt.plot(blu blu)
"""
    ax = fig.add_subplot(1,1,1,frameon=False)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='none', top='off', bottom='off',
    left='off', right='off') 
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if xlabelpad:
        ax.xaxis.labelpad = xlabelpad
    if ylabelpad:
        ax.yaxis.labelpad = ylabelpad
    ax.set_xlabel(xlabel,fontsize=fontsize)
    ax.set_ylabel(ylabel,fontsize=fontsize)
    #ax.xaxis.set_label_coords(0.5,-0.06)
    ax.set_title(title)


def plot_sq(x0,y0,side,**kwargs):
    """Plots a square of side centred at (x0,y0)."""
    
    x1 = x0-side/2.
    y1 = y0-side/2.
    x2 = x0+side/2.
    y2 = y1
    x3 = x2
    y3 = y0+side/2.
    x4 = x1
    y4 = y3
    
    plt.fill([x1,x2,x3,x4],[y1,y2,y3,y4],**kwargs)
    

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1."""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)


def plot_spectrum(ax, spec, complist=None, plot_res=True, fwhm=3):
    """ Plots a spectrum in a given axis

    Parameters
    ----------
    ax : axis
        matplotlib axis
    spec: XSpectrum1D
        Spectrum to plot
    complist : list of AbsComponents, optional
        If given, a model is computed and plotted
    plot_res : bool, optional
        Whether to plot residuals, only works if components
        is not None.
    fwhm : int, optional
        FWHM in pixels

    Returns
    -------
    ax : axis
        Axis with the plot.
    """
    # checks
    if not isinstance(spec, XSpectrum1D):
        raise IOError('Input spec must be XSpectrum1D object.')
    if complist is not None:
        if not isinstance(complist[0], AbsComponent):
            raise IOError('components must be a list of AbsComponent objects.')
        plot_model = True
    else:
        plot_model = False

    if plot_model:
        # create a model
        model = lav.voigt_from_components(spec.wavelength, complist, fwhm=fwhm)

    if spec.co_is_set:
        spec.normalize(co=spec.co)
    ax.plot(model.wavelength, spec.flux, '-', color=COLOR_FLUX, lw=1)
    if spec.sig_is_set:
        ax.plot(model.wavelength, spec.sig, '-', color=COLOR_SIG, lw=0.5)
    if plot_model:
        ax.plot(model.wavelength, model.flux, '-', color=COLOR_MODEL, lw=0.5)
        if plot_res:
            residual = spec.flux - model.flux
            ax.plot(spec.wavelength, residual, '.', color=COLOR_RESIDUAL, ms=2)
            ax.plot(spec.wavelength, -1*spec.sig, '-', drawstyle='steps-mid', color=COLOR_SIG, lw=0.5)


def plot_vel(ax, spec, iline, z, dvlims, complist=None, fwhm=3):

    # first normalize
    if not spec.normed:
        spec.normalize(co=spec.co)

    # first identify good components
    good_comps = []
    good_comps_aux = pntu.get_components_at_z(complist, z, dvlims)
    for comp in good_comps_aux:
        for aline in comp._abslines:
            if aline.name == iline['name']:
                good_comps += [comp]
                break
    # bad comps will be those unrelated to the given redshift/transition
    bad_comps = []
    for comp in complist:
        if comp not in good_comps:
            bad_comps += [comp]
    # import pdb; pdb.set_trace()
    # only good comps will have a model
    if len(good_comps) > 0:
        model_spec = lav.voigt_from_components(spec.wavelength, good_comps, fwhm=fwhm)
    else:
        model_spec = XSpectrum1D.from_tuple((spec.wavelength, np.ones(spec.npix)))

    # main plot
    velo = ltu.give_dv(spec.wavelength/iline['wrest'] - 1. , z)
    ax.plot(velo, spec.flux, drawstyle='steps-mid', color=COLOR_FLUX)
    plot_spec_complist(ax, velo, spec, bad_comps, min_ew = 0.3*u.AA, color=COLOR_BLEND, lw=2, drawstyle='steps-mid')
    plot_spec_complist(ax, velo, model_spec, good_comps, min_ew = None, color=COLOR_MODEL, lw=1, ls='-')

    if spec.sig_is_set:
        ax.plot(velo, spec.sig, drawstyle='steps-mid', color=COLOR_SIG, lw=0.5)

    # Zero velocity line
    ax.plot([0., 0.], [-1e9, 1e9], ':', color='gray')
    # Unity flux level line
    ax.plot([-1e9, 1e9], [1, 1], ':', color='b', lw=0.5)
    # Zero flux level line
    ax.plot([-1e9, 1e9], [0, 0], '--', color='k', lw=1)


def plot_spec_comp(ax, x, spec, comp, min_ew=None, label=False, **kwargs):
    """Plots absorption lines within a component over a given spectrum model"""

    for aline in comp._abslines:
        # check min ew if given
        if min_ew is not None:
            try:
                b = aline.attrib['b']
            except:
                b = 10 * u.km/u.s
            Wr = aline.get_Wr_from_N_b(aline.attrib['N'], b)
            if Wr < min_ew:
                pass
        wvlim = aline.limits.wvlim
        cond = (spec.wavelength > wvlim[0]) & (spec.wavelength < wvlim[1])
        ax.plot(x[cond], spec.flux[cond], **kwargs)
        if label:
            # import pdb; pdb.set_trace()
            s = '{}z={:.3f}'.format(aline.name, aline.z)
            x_s = np.mean(x[cond].value)
            y_s = 0.5
            ax.annotate(s, (x_s, y_s), rotation=90, fontsize=6)

def plot_spec_complist(ax, x, spec, complist, min_ew=None, **kwargs):
    for comp in complist:
        plot_spec_comp(ax, x, spec, comp, min_ew=min_ew, **kwargs)

