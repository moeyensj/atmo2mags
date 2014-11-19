import os
import numpy
import pylab
from lsst.sims.catalogs.measures.photometry.Sed import Sed
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass

# Find aerosol files
tmp = os.listdir('.')
ozonefiles = []
ozone = []
for t in tmp:
    if t.startswith('am1.2') and t.endswith('.plt'):
        ozonefiles.append(t)
        ozone.append(float(t.split('_')[2].strip('oz').strip('.plt')))

ozone = numpy.array(ozone)

# Read in aerosol files.
ozones = {}
for a in ozonefiles:
    print 'Reading %s' %(a)
    ozones[a] = Bandpass()
    ozones[a].readThroughput(a)

"""
# plot files
for a in ozonefiles:
    pylab.figure()
    pylab.plot(ozones[a].wavelen, ozones[a].sb)
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('Transmission')
    pylab.savefig(a+'.png', format='png')
"""
