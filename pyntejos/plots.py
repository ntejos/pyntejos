import matplotlib.pyplot as plt
from astropy import units as u
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.isgm.abscomponent import AbsComponent
from linetools.analysis import voigt as lav

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
        gdlin = []
        for comp in complist:
            for ii, line in enumerate(comp._abslines):
                wvobs = (1 + line.attrib['z']) * line.wrest
                if (wvobs > spec.wvmin) & (wvobs <spec.wvmax):
                    line.attrib['N'] = 10.**line.attrib['logN'] / u.cm**2
                    gdlin.append(line)
        model = lav.voigt_from_abslines(spec.wavelength, gdlin, fwhm=fwhm)


    if spec.co_is_set:
        spec.normalize(co=spec.co)
    ax.plot(spec.wavelength, spec.flux, '-', drawstyle='steps-mid', color='k')
    if spec.sig_is_set:
        ax.plot(spec.wavelength, spec.sig, '-', drawstyle='steps-mid', color='g', lw=0.5)
    if plot_model:
        ax.plot(model.wavelength, model.flux, '-', color='r', lw=0.5)
        if plot_res:
            residual = spec.flux - model.flux
            ax.plot(spec.wavelength, residual, '.', color='grey', ms=2)
            ax.plot(spec.wavelength, -1*spec.sig, '-', drawstyle='steps-mid', color='g', lw=0.5)
            





