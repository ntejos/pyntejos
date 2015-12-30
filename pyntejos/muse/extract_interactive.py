# Extract spectra and 2d images of objects in a MUSE data cube.
# Usage: extract.py [filename] [sexcat] [MmmYYYY]
# SExtractor cataogue required and must contain x, y coordinates (Cols 6 & 7), flags (Col10) and 90% flux radius (Col15).

import numpy as np

class PointBrowser:
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """
    def __init__(self):
        self.lastind = 0

        self.text = ax2.text(0.05, 0.95, 'selected: none',transform=ax2.transAxes, va='top')
        xs = wl[0]
        ys = spec_smooth[0]
        if redshift[i] > 2.96:
            xs = (redshift[i]+1.)*1215.67
            distances = np.abs(wl-xs)
            indmin = distances.argmin()
            ys = spec_smooth[indmin]
        print xs,ys
        self.selected, = ax2.plot([xs,xs], [0,ys], '-',color='red')

    def onpress(self, event):
        global id
        if self.lastind is None: return
        if event.key not in ('1', '2','3','9','a','b','C','d','m','n','p'): return
        if event.key=='1':
            var = raw_input("Wavelength for Lya line: ")
            redshift[i] = float(var)/1215.67-1.
            self.update()
        if event.key=='2':
            var = raw_input("Wavelength for Oii line: ")
            redshift[i] = float(var)/3727.-1.
            self.update()
        if event.key=='3':
            var = raw_input("Wavelength for Oiii line: ")
            redshift[i] = float(var)/5007.-1.
            self.update()
        if event.key=='9':
            var = raw_input("Wavelength for HK break: ")
            redshift[i] = float(var)/4000.-1.
            self.update()
        if event.key=='a':
            var = raw_input("Wavelength for Halpha: ")
            redshift[i] = float(var)/6562.8-1.
            self.update()
        if event.key=='b':
            var = raw_input("Wavelength for HBeta: ")
            redshift[i] = float(var)/4861.-1.
            self.update()
        if event.key=='d':
            var = raw_input("Wavelength for HDelta: ")
            redshift[i] = float(var)/4101.7-1.
            self.update()
        if event.key=='C':
            var = raw_input("Wavelength for CIII (1907AA): ")
            redshift[i] = float(var)/1907.-1.
            self.update()
        if event.key=='m':
            var = raw_input("Wavelength for MgII (2796AA): ")
            redshift[i] = float(var)/2796.-1.
            self.update()
        if event.key=='n':
            id=id+1
            self.replot()
        if event.key=='p':
            id=id-1
            self.replot()

    def onpick(self, event):
        print "You clicked!"
        if event.artist!=line: return True
        N = len(event.ind)
        if not N: return True
            # the click locations
        xclick = event.mouseevent.xdata
        yclick = event.mouseevent.ydata
        print xclick,yclick
        distances = np.hypot(xclick-wl[event.ind], yclick-smooth_spec[event.ind])
        indmin = distances.argmin()
        dataind = event.ind[indmin]
        self.lastind = dataind
    #        self.update()

    def update(self):
        global id,wl,spec_smooth,perc
        if self.lastind is None: return
        dataind = self.lastind

        ax2.cla()
        plt.xlabel(r'$\lambda$ ($\AA$)',fontsize=18)
        plt.ylabel(r'Flux ($10^{-18}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)',fontsize=18)
        plt.axis([4680,9280,-0.05*perc,1.05*perc])
        self.skylines()
        ax2.plot([4680,9280],[0.,0.],linestyle=':',zorder=2,color='k')
        ax2.plot(wl,spec_smooth,label='{0:03d},{1:6.4f}'.format(id,redshift[i]))
        lindex = 0
        for line in lines:
            if (line*(1.+redshift[i]) > 4680) & (line*(1.+redshift[i]) < 9280):
                sep = np.abs(wl-line*(1.+redshift[i]))
                lind = sep.argmin()
                ypos = np.min([spec_smooth[lind]+0.02*perc,0.86*perc])
                if ltype[lindex] == 1:
                    ax2.plot(np.array([1,1])*line*(1.+redshift[i]),np.array([ypos,ypos+0.12*perc]),color='red')
                    ax2.text(line*(1.+redshift[i]),ypos+0.14*perc, names[lindex], ha="center", size=10,color='red')
                if ltype[lindex] == 2:
                    ax2.plot(np.array([1,1])*line*(1.+redshift[i]),np.array([ypos,ypos+0.08*perc]),color='blue')
                    ax2.text(line*(1.+redshift[i]),ypos+0.10*perc, names[lindex], ha="center", size=10,color='blue')
            lindex = lindex+1
        ax2.plot([4680,9280],[0,0],linestyle='..',color='k')
        ax2.legend()

        fig.canvas.draw()

    def replot(self):
        global id,wl,spec_smooth,perc,halfwidth
        ax1.cla()
        ax2.cla()
        off = np.abs(id-ids)
        i   = off.argmin()
        print "{0:5d} {1:9.3f} {2:9.3f} {3:4d} {4:12.8f}".format(id,x[i],y[i],flag[i],frad90[i])
        print id,frad90[i],y[i],x[i]
        if (2.*frad90[i]-x[i] > 0.) | (2.*frad90[i]-y[i] > 0): frad90[i] = np.min([x[i]/2.-1,y[i]/2.-1])
        print id,frad90[i],y[i],x[i]
        if (frad90[i] > 0.5) & (x[i] > 12) & (y[i] > 12):
            halfwidth = int(np.max([2.*frad90[i],12]))
            print halfwidth
            minicube = cube[:,y[i]-halfwidth:y[i]+halfwidth,x[i]-halfwidth:x[i]+halfwidth]
            shape = minicube.shape
            image    = np.nansum(minicube,axis=0)
            spectrum = np.zeros(naxis3)
            bgd = np.zeros(naxis3)
            radial_mask = minicube[2000,:,:]*0.
            bgd_mask    = minicube[2000,:,:]*0.
            print radial_mask.shape,minicube.shape
            for j in xrange(shape[1]):
                for k in xrange(shape[2]):
                    if ((j-halfwidth)**2+(k-halfwidth)**2)**0.5 < frad50[i]:
                        radial_mask[j,k] = 1.
                    if ((j-halfwidth)**2+(k-halfwidth)**2)**0.5 > frad90[i]:
                        bgd_mask[j,k] = 1.
            for j in xrange(naxis3):
                masked_square = minicube[j,:,:]*radial_mask
                spectrum[j] = np.nansum(minicube[j,:,:]*radial_mask)
                bgd[j]      = np.nanmedian(minicube[j,:,:]*bgd_mask)*np.nansum(radial_mask)
            spectrum = spectrum-bgd
            bgd_smooth = bgd
            spec_smooth = spectrum
            for k in xrange(len(spectrum)-smsc):
                spec_smooth[k+smsc/2] = np.median(spectrum[k:k+smsc])
                bgd_smooth[k+smsc/2] = np.median(bgd[k:k+smsc])
            mincut = np.percentile(image,5.)
            maxcut = np.percentile(image,95.)
            image[image < mincut] = mincut
            image[image > maxcut] = maxcut
            ax1.imshow(image,interpolation="bessel")
            circle1 = plt.Circle((halfwidth,halfwidth),frad90[i],color='r',fill=False)
            ax1.add_artist(circle1)
            if np.max(spec_smooth) < 400:
                perc = np.max([np.percentile(spec_smooth,99.5),np.median(spec_smooth)*4.])
            else:
                perc = np.max(spec_smooth)
            ax2.cla()
            plt.xlabel(r'$\lambda$ ($\AA$)',fontsize=18)
            plt.ylabel(r'Flux ($10^{-18}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)',fontsize=18)
            plt.axis([4680,9280,-0.05*perc,1.05*perc])
            self.skylines()
            ax2.plot([4680,9280],[0.,0.],linestyle=':',zorder=2,color='k')
            ax2.plot(wl,spec_smooth,label='{0:03d},{1:6.4f}'.format(id,redshift[i]),zorder=4)
            ax2.plot(wl,bgd_smooth,zorder=1,color='yellow')
            lindex = 0
            for line in lines:
                if (line*(1.+redshift[i]) > 4680) & (line*(1.+redshift[i]) < 9280):
                    sep = np.abs(wl-line*(1.+redshift[i]))
                    lind = sep.argmin()
                    ypos = np.min([spec_smooth[lind]+0.02*perc,0.86*perc])
                    if ltype[lindex] == 1:
                        ax2.plot(np.array([1,1])*line*(1.+redshift[i]),np.array([ypos,ypos+0.12*perc]),color='red')
                        ax2.text(line*(1.+redshift[i]),ypos+0.14*perc, names[lindex], ha="center", size=10,color='red')
                    if ltype[lindex] == 2:
                        ax2.plot(np.array([1,1])*line*(1.+redshift[i]),np.array([ypos,ypos+0.08*perc]),color='blue')
                        ax2.text(line*(1.+redshift[i]),ypos+0.10*perc, names[lindex], ha="center", size=10,color='blue')
                lindex = lindex+1
            ax2.plot([4680,9280],[0,0],linestyle='..',color='k')
            ax2.legend()

            fig.canvas.draw()
        else: print "Warning: Object too close to image edge, not plotted."


    def skylines(self):
        skylines = [5577.,5890.,6310.,6828.,7243.,7526.,7720.,7780.,7920.,8300.,8636.,8675.,8770.,8850.]
        for sline in skylines:
            ax2.fill_between([sline-5.,sline+5.],[0,0],[0.64*perc,0.64*perc],color='#A692AE',zorder=1)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import argparse
    import pyfits
    import sys
    import matplotlib.gridspec as gridspec

#    def ontype(event):
#        """Deal with keyboard events"""
#        print("You pressed key {:s}".format(event.key))
#        if event.key == 'K':
#            sys.exit()
#        if event.key == 'l':


    parser = argparse.ArgumentParser(description='Process some strings.')
    parser.add_argument('filename', metavar='N', type=str, nargs='+',
                        help='Input MUSE datacube')
    parser.add_argument('sexcat', metavar='N', type=str, nargs='+',
                        help='SExtractor catalogue of sources from the data cube.')
    parser.add_argument('date', metavar='N', type=str, nargs='+',
                        help='Date in format MmmYYYY.')
    args = parser.parse_args()

    print 'Reading from ',args.sexcat[0]
    ids,x,y,flag,frad50,frad90,redshift = np.genfromtxt(args.sexcat[0],unpack=True,usecols=(0,5,6,9,13,14,15),comments='#')
    x = x
    y = y
    ids = ids.astype(int)
    flag = flag.astype(int)
    nobj = len(ids)

    # Read in atomic line catalogue.
    lines,ltype = np.genfromtxt('galaxylines_py.dat',unpack=True,usecols=(0,1),comments='#')
    names = np.genfromtxt('galaxylines_py.dat',unpack=True,usecols=(2),dtype="S25",comments='#')

    # Brace and read in MUSE datacube:
    print "Reading in datacube: ",args.filename[0]
    hdulist = pyfits.open(args.filename[0])
    cube    = hdulist[1].data
    print cube.shape
    crval3  = hdulist[1].header['crval3']
    crpix3  = hdulist[1].header['crpix3']
    cd3_3   = hdulist[1].header['cd3_3']
    naxis3  = hdulist[1].header['naxis3']
    wl      = (np.arange(naxis3)+1-crpix3)*cd3_3+crval3
    wlplt   = (np.arange(naxis3/2)+1-crpix3/2.)*cd3_3*2.+crval3
    #fullimage = np.nansum(cube,axis=0)
    #mincut = np.percentile(fullimage,5.)
    #maxcut = np.percentile(fullimage,95.)
    #fullimage[fullimage < mincut] = mincut
    #fullimage[fullimage > maxcut] = maxcut
    #plt.figure(1)
    #plt.imshow(fullimage,interpolation='bessel')
    #fig = plt.gcf()
    #fig.set_size_inches(6.4,6.4)
    #plt.savefig('plots/Oct2014/image_Q0302_full.pdf',format='pdf',bbox_inches='tight')

    pixel = np.where(np.abs(wl-(3.29+1.)*1215.67) == np.min(np.abs(wl-(3.29+1.)*1215.67)))[0]

    fullimage = np.nansum(cube[pixel-5:pixel+5,:,:],axis=0)
    mincut = np.percentile(fullimage,5.)
    maxcut = np.percentile(fullimage,95.)
    fullimage[fullimage < mincut] = mincut
    fullimage[fullimage > maxcut] = maxcut

    smsc = 4

    xsize = len(cube[0,0,:])
    ysize = len(cube[0,:,0])
    nproc = 0
    id = 1
    while id < nobj:
        off = np.abs(id-ids)
        i = off.argmin()
        print "{0:5d} {1:9.3f} {2:9.3f} {3:4d} {4:12.8f}".format(id,x[i],y[i],flag[i],frad90[i])
        nproc = nproc + 1
        if (2.*frad90[i]-x[i] > 0.) | (2.*frad90[i]-y[i] > 0): frad90[i] = np.min([x[i]/2.-1,y[i]/2.-1])
        print id,frad90[i],y[i],x[i]
        if (frad90[i] > 0.5) & (x[i] > 12) & (y[i] > 12):
            halfwidth = int(np.max([2.*frad90[i],12]))
            minicube = cube[:,int(y[i])-halfwidth:int(y[i])+halfwidth,int(x[i])-halfwidth:int(x[i])+halfwidth]
            shape = minicube.shape
            image    = np.nansum(minicube,axis=0)
            spectrum = np.zeros(naxis3)
            bgd = np.zeros(naxis3)
            radial_mask = minicube[2000,:,:]*0.
            bgd_mask    = minicube[2000,:,:]*0.
            print radial_mask.shape,minicube.shape
            for j in xrange(shape[1]):
                for k in xrange(shape[2]):
                    if ((j-halfwidth)**2+(k-halfwidth)**2)**0.5 < frad50[i]:
                        radial_mask[j,k] = 1.
                    if ((j-halfwidth)**2+(k-halfwidth)**2)**0.5 > frad90[i]:
                        bgd_mask[j,k] = 1.
            for j in xrange(naxis3):
                masked_square = minicube[j,:,:]*radial_mask
                spectrum[j] = np.nansum(minicube[j,:,:]*radial_mask)
                bgd[j]      = np.nanmedian(minicube[j,:,:]*bgd_mask)*np.nansum(radial_mask)
            spectrum = spectrum-bgd
            bgd_smooth = bgd
            spec_smooth = spectrum
            for k in xrange(len(spectrum)-smsc):
                spec_smooth[k+smsc/2] = np.median(spectrum[k:k+smsc])
                bgd_smooth[k+smsc/2] = np.median(bgd[k:k+smsc])
            fig = plt.figure(1,figsize=(16,6))
            gs     = gridspec.GridSpec(1,4)
            plt.subplots_adjust(left=0.04, right=0.98, wspace=0.47)
            ax1    = plt.subplot(gs[:,0])
            mincut = np.percentile(image,5.)
            maxcut = np.percentile(image,95.)
            image[image < mincut] = mincut
            image[image > maxcut] = maxcut
            ax1.imshow(image,interpolation="bessel")
            ax2 = plt.subplot(gs[:,1:])
            if np.max(spec_smooth) < 400:
                perc = np.percentile(spec_smooth,98.)
            else:
                perc = np.max(spec_smooth)
            plt.xlabel(r'$\lambda$ ($\AA$)',fontsize=18)
            plt.ylabel(r'Flux ($10^{-18}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)',fontsize=18)
            plt.axis([4680,9280,-0.05*perc,1.05*perc])
            skylines = [5577.,5890.,6310.,6828.,7243.,7526.,7720.,7780.,7920.,8300.,8636.,8675.,8770.,8850.]
            for sline in skylines:
                ax2.fill_between([sline-5.,sline+5.],[0,0],[0.64*perc,0.64*perc],color='#A692AE',zorder=1)
            ax2.plot([4680,9280],[0.,0.],linestyle=':',zorder=2,color='k')
            ax2.plot(wl,spec_smooth,label='{0:03d},{1:6.4f}'.format(id,redshift[i]),zorder=4)
            ax2.plot(wl,bgd_smooth,zorder=1,color='yellow')
            lindex = 0
            for line in lines:
                sep = np.abs(wl-line*(1.+redshift[i]))
                lind = sep.argmin()
                ypos = np.min([spec_smooth[lind]+0.02*perc,0.86*perc])
                if ltype[lindex] == 1:
                    ax2.plot(np.array([1,1])*line*(1.+redshift[i]),np.array([ypos,ypos+0.12*perc]),color='red')
                    ax2.text(line*(1.+redshift[i]),ypos+0.14*perc, names[lindex], ha="center", size=10,color='red')
                if ltype[lindex] == 2:
                    ax2.plot(np.array([1,1])*line*(1.+redshift[i]),np.array([ypos,ypos+0.08*perc]),color='blue')
                    ax2.text(line*(1.+redshift[i]),ypos+0.10*perc, names[lindex], ha="center", size=10,color='blue')
                lindex = lindex+1
            ax2.plot([4680,9280],[0,0],linestyle='..',color='k')
            ax2.legend()

            browser = PointBrowser()
            fig.canvas.mpl_connect('pick_event', browser.onpick)
            fig.canvas.mpl_connect('key_press_event', browser.onpress)
            plt.show()
        else: print "Warning! Skipping object as too near to image boundary!"
        var = raw_input("enter a to stop: ")
        id = int(var)

    print nobj,i
    print 'Done ',nproc
