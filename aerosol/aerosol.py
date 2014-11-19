import os
import numpy
import pylab
from lsst.sims.catalogs.measures.photometry.Sed import Sed
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass

# Find aerosol files
tmp = os.listdir('.')
aerosolfiles = []
aerosol_mult = []
for t in tmp:
    if t.startswith('am1.2') and t.endswith('.plt'):
        aerosolfiles.append(t)
        aerosol_mult.append(float(t.split('_')[2].strip('avis').strip('.plt'))/10.0)

aerosol_mult = numpy.array(aerosol_mult)

# Read in aerosol files.
aerosols = {}
for a in aerosolfiles:
    print 'Reading %s' %(a)
    aerosols[a] = Bandpass()
    aerosols[a].readThroughput(a)


# plot files
for a in aerosolfiles:
    pylab.figure()
    pylab.plot(aerosols[a].wavelen, aerosols[a].sb)
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('Transmission')
    pylab.savefig(a+'.png', format='png')
