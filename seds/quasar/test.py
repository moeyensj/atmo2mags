import copy
import numpy
import pylab
import lsst.sims.catalogs.measures.photometry.Sed as Sed


wavelen= numpy.arange(12, 1200, 1.0)
flat = Sed(wavelen, flambda=wavelen*0.0+1.0)
redshifts = numpy.arange(0.1, 6.0, 1.5)
#redshifts=[3.0,]
for z in redshifts:
    flat2 = copy.deepcopy(flat)
    flat2.redshiftSED(z)
    flat3 = copy.deepcopy(flat2)
    rtau = 0.50 #/ ((1+z)**(3/2))
    print z, rtau
    flat3.addIGMattenuation(z, rtau)

    pylab.plot(flat2.wavelen, flat2.flambda, 'k-')
    pylab.plot(flat3.wavelen, flat3.flambda, label='%.1f' %(z))
    pylab.axvline(91.2*(1+z), color='r', linestyle=':')
    
pylab.legend(loc='upper right')
pylab.xlabel("Wavelength (nm)")
pylab.ylabel("Flambda .. ==attenuation")
pylab.xlim(30, 1200)
pylab.show()

