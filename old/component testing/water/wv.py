import os
import numpy
import pylab
from scipy import interpolate
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import plot_dmags as pd

def read_wvfiles():
    tt = os.listdir('.')
    files = []
    for t in tt:
        if t.startswith('am') & t.endswith('.plt'):
            files.append(t)
    print files
    atms = {}
    wvs = []
    for f in files:
        wv = f[9:12]
        print wv, f
        wvs.append(wv)
        atms[wv] = Bandpass()
        atms[wv].readThroughput(f, wavelen_step=0.1)
    std = Bandpass()
    std.readThroughput(os.path.join(os.getenv('LSST_THROUGHPUTS_ATMOS'), 'atmos_12.dat'))
    #std.readThroughput('atmos_85143574.dat', wavelen_step=.1)

    # calculate effect of absorption ..
    condition = ((atms[wvs[0]].wavelen >= 850) & (atms[wvs[0]].wavelen <= 1100))
    abs = {}
    for w in wvs:
        abs[w] = (1-atms[w].sb[condition]).sum() * 0.1
        print w, abs[w]

    # check with a plot
    max = len(wvs)-1
    pylab.figure()
    for w in [wvs[0], wvs[max]]:
        pylab.plot(atms[w].wavelen, atms[w].sb, label=w)
    pylab.plot(std.wavelen, std.sb, 'k-', label='std')
    pylab.legend(loc='lower right')
    pylab.xlabel('Wavelength')
    pylab.ylabel('Transmission')
    pylab.xlim(900, 1000)
    pylab.figure()
    for w in [wvs[0], wvs[max]]:
        pylab.plot(atms[w].wavelen, atms[w].sb, label=w)
    pylab.plot(std.wavelen, std.sb, 'k-', label='std')
    pylab.legend(loc='lower right')
    pylab.xlabel('Wavelength')
    pylab.ylabel('Transmission')
    #pylab.show()

    return atms, wvs, std



if __name__ == '__main__':

    filterlist = ('u', 'g', 'r', 'i', 'z', 'y3', 'y4')
    sys = pd.read_hardware(flist=filterlist)    
    atms, wvs, std = read_wvfiles()
    
    w = wvs[0]
    lsst = pd.combine_throughputs(atms[w], sys, flist=filterlist)

    stars, starlist, temperature, metallicity, logg = pd.read_kurucz()
    refmags = pd.calc_mags(stars, starlist, lsst, filterlist)

    gi = refmags['g'] - refmags['i']

    mags = {}
    dmags = {}
    for w in wvs:
        lsst = pd.combine_throughputs(atms[w], sys, flist=filterlist)
        mags[w] = pd.calc_mags(stars, starlist, lsst, filterlist)
        print w
    refmags = mags[wvs[0]]

    for w in wvs:
        dmags[w] = pd.calc_deltamags(mags[w], refmags, flist=filterlist)
        

    ylims ={}
    ylims['y3'] = [-5, 6]
    ylims['y4'] = [-5, 6]


    nwvs = numpy.array(wvs, 'float')
    (gii, wii) = numpy.meshgrid(gi, nwvs)
    dy3 = numpy.zeros([len(wvs), len(gi)], 'float')
    dy4 = numpy.zeros([len(wvs), len(gi)], 'float')
    #print numpy.shape(wii), numpy.shape(gii), numpy.shape(dy3), numpy.shape(dy4)
    #print numpy.shape(dmags[wvs[0]]['y3'])
    for i in range(len(wvs)):
        dy3[i] = dmags[wvs[i]]['y3']
        dy4[i] = dmags[wvs[i]]['y4']    
    
    y3ii = interpolate.LinearNDInterpolator((gii.ravel(), wii.ravel()), dy3.ravel())
    #y3ii = interpolate.SmoothBivariateSpline(gii.ravel(), wii.ravel(), dy3.ravel(), s=3)
    y4ii = interpolate.LinearNDInterpolator((gii.ravel(), wii.ravel()), dy4.ravel())
    #y4ii = interpolate.SmoothBivariateSpline(gii.ravel(), wii.ravel(), dy4.ravel(), s=3)



    # Plot dmag as a function of gi color, all stars and all WV's
    newfig=True
    for w in wvs:
        print w
        pd.plot_dmags(gi, dmags[w], [metallicity, logg], 'kurucz', plotfilterlist=('y3',), ylims=ylims,
                      titletext="Water Vapor", newfig=newfig)
        newfig=False


    newfig=True
    for w in wvs:
        print w
        pd.plot_dmags(gi, dmags[w], [metallicity, logg], 'kurucz', plotfilterlist=('y4',), ylims=ylims,
                      titletext="Water Vapor", newfig=newfig)
        newfig=False

    # Plot dmag as function of WV for a particular g-i color
    g2 = numpy.array([-0.5,], 'float')
    wv2 = numpy.arange(1, 8, 1.0)
    (g2ii, wv2ii) = numpy.meshgrid(g2, wv2)
    y3val = y3ii((g2ii, wv2ii))
    y4val = y4ii((g2ii, wv2ii))
    calcdy3 = numpy.transpose(y3val)[0]
    calcdy4 = numpy.transpose(y4val)[0]
    print numpy.shape(wv2), numpy.shape(y3val), numpy.shape(y3val[0]), numpy.shape(numpy.transpose(y3val))    
    pylab.figure()
    pylab.plot(wv2, calcdy3, 'b-', marker='o', label='y3')
    pylab.plot(wv2, calcdy4, 'r-', marker='o', label='y4')
    pylab.legend(loc='upper left', numpoints=1, fancybox=True, shadow=True)
    pylab.xlabel('Water Absorption, mm')
    pylab.ylabel('Change in (natural) magnitude, for g-i=%0.1f' %(g2[0]))

    
    # Generate plot of interpolated data, on grid - all dmag as heatmap with WV and g-i axes
    g2 = numpy.arange(-1, 1.6, 0.1)
    wv2 = numpy.arange(0.9, 7.1, 0.1)
    (g2ii, wv2ii) = numpy.meshgrid(g2, wv2)
    y3grid = y3ii((g2ii, wv2ii))
    #y3grid = y3ii(g2ii, wv2ii)
    pylab.figure()
    vlim = ylims['y4']
    gilim = [-1, 1.5]
    wvlim = [1, 7]
    pylab.imshow(y3grid, extent=(wvlim[0], wvlim[1], gilim[0], gilim[1]), origin='lower', vmin=vlim[0], vmax=vlim[1])
    cb = pylab.colorbar(orientation='vertical', shrink=0.6)
    cb.set_label('mmag')
    pylab.ylabel('g-i color (mag)')
    pylab.xlabel('Precipitable Water Vapor (mm)')
    pylab.title('Change in y3 as function of water vapor and SED color')

    pylab.figure()
    pylab.imshow(y3grid, extent=(wvlim[0], wvlim[1], gilim[0], gilim[1]), origin='lower')
    cb = pylab.colorbar(orientation='vertical', shrink=0.6)
    cb.set_label('mmag')
    pylab.ylabel('g-i color (mag)')
    pylab.xlabel('Precipitable Water Vapor (mm)')
    pylab.title('Change in y3 as function of water vapor and SED color')

    y4grid = y4ii((g2ii, wv2ii))
    #y4grid = y4ii(g2ii, wv2ii)
    pylab.figure()
    pylab.imshow(y4grid, extent=(wvlim[0], wvlim[1], gilim[0], gilim[1]), origin='lower', vmin=vlim[0], vmax=vlim[1])
    cb = pylab.colorbar(orientation='vertical', shrink=0.6)
    cb.set_label('mmag')
    pylab.ylabel('g-i color (mag)')
    pylab.xlabel('Precipitable Water Vapor (mm)')
    pylab.title('Change in y4 as function of water vapor and SED color')


    
    pylab.show()
