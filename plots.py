import matplotlib.pyplot as pl

def common_labels(fig,xlabel='',ylabel='',title='',fontsize=20,xlabelpad=None,ylabelpad=None):
    """Plots only labels in a parent figure. Special for subplots with
    shared axes.
    
    Mode of use:

    fig  = pl.figure()
    common_labels(fig, xlabel='x label',ylabel='y label',fontsize=20)
    
    ax1 = fig.add_subplot(1,2,1)
    pl.plot(bla bla)
    ax2 = fig.add_subplot(1,2,2)
    pl.plot(blu blu)
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
    
    pl.fill([x1,x2,x3,x4],[y1,y2,y3,y4],**kwargs)
    

