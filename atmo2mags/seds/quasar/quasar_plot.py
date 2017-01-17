import numpy
import pylab
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import os

# read filters
filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
lsst = {}
for f in filterlist:
    lsst[f] = Bandpass()
    lsst[f].readThroughput(os.path.join(filterdir, "total_" + f + ".dat"))
    
# read base quasar sed
qso = Sed()
qso.readSED_flambda("quasar.dat")

# set redshift range
redshifts = numpy.arange(0, 6, 0.2)

# set up quasar seds for this redshift range
qsoz = {}
for z in redshifts:
    qsoz[z] = Sed(qso.wavelen, flambda=qso.flambda)
    qsoz[z].redshiftSED(z, dimming=True)
    qsoz[z].flambdaTofnu()

# calculate colors
mags = {}
for f in filterlist:
    mags[f] = {}
    for z in redshifts:
        mags[f][z] = qsoz[z].calcMag(lsst[f])

# plot colorcolor plots like in sci book
colors = ['b', 'c', 'g', 'y', 'r', 'm', 'k']
pylab.figure()
pylab.subplot(221)
colorindex = 0
zprev = 0
for z in redshifts:
    pylab.plot(mags['u'][z]-mags['g'][z], mags['g'][z]-mags['r'][z], colors[colorindex]+'o')
    if (z-zprev)>=1:
        colorindex  = colorindex+1
        zprev = z
    #pylab.annotate(z, (mags['u'][z]-mags['g'][z], mags['g'][z]-mags['r'][z]+.1))
pylab.xlabel("u-g")
pylab.ylabel("g-r")
pylab.xlim(-.5, 4)
pylab.ylim(-0.5, 2)
pylab.grid()
#
pylab.subplot(222)
colorindex = 0
zprev = 0
for z in redshifts:
    pylab.plot(mags['g'][z]-mags['r'][z], mags['r'][z]-mags['i'][z], colors[colorindex]+'o')
    if (z-zprev)>=1:
        colorindex  = colorindex+1
        zprev = z
    #pylab.annotate(z, (mags['g'][z]-mags['r'][z], mags['r'][z]-mags['i'][z]+.1))
pylab.xlabel("g-r")
pylab.ylabel("r-i")
pylab.xlim(-.5, 2)
pylab.ylim(-.5, 2.5)
pylab.grid()
#
pylab.subplot(223)
colorindex = 0
zprev = 0
for z in redshifts:
    pylab.plot(mags['r'][z]-mags['i'][z], mags['i'][z]-mags['z'][z], colors[colorindex]+'o')
    if (z-zprev)>=1:
        colorindex  = colorindex+1
        zprev = z
    #pylab.annotate(z, (mags['r'][z]-mags['i'][z], mags['i'][z]-mags['z'][z]+.1))
pylab.xlabel("r-i")
pylab.ylabel("i-z")
pylab.xlim(-.5, 2.5)
pylab.ylim(-.5, 1.5)
pylab.grid()
#
pylab.subplot(224)
colorindex = 0
zprev = 0
for z in redshifts:
    pylab.plot(mags['i'][z]-mags['z'][z], mags['z'][z]-mags['y'][z], colors[colorindex]+'o')
    if (z-zprev)>=1:
        print zprev, z, colorindex, colors[colorindex]
        colorindex  = colorindex+1
        zprev = z
    #pylab.annotate(z, (mags['r'][z]-mags['i'][z], mags['i'][z]-mags['z'][z]+.1))
pylab.xlabel("i-z")
pylab.ylabel("z-y")
pylab.xlim(-.5, 1.5)
pylab.ylim(-1, 1)
pylab.grid()

pylab.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)

pylab.show()
