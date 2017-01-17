# Python packages
import numpy as np

# LSST stack software
import lsst.sims.photUtils.Bandpass as Bandpass

WAVELENMIN = 300
WAVELENMAX = 1100
WAVELENSTEP = 0.5

class Atmo(object):
    def __init__(self, parameters, airmass, transmission, wavelength, aerosolNormCoeff, aerosolNormWavelen):
        # Airmass
    	self.X = airmass
        # List of parameters
    	self.P = parameters
        # List of wavelengths
        self.wavelen = wavelength
        # Aerosol normalization coefficient used to model aerosol transmission profile
    	self.aerosolNormCoeff = aerosolNormCoeff
        # Aerosol normalization wavelength used to model aerosol transmission profile
    	self.aerosolNormWavelen = aerosolNormWavelen
        # List of components
    	self.components = ['H2O','O2','O3','Rayleigh','Aerosol']
        # List of total transmission profiles
    	self.sb = None
        # Component-keyed dictionary of transmission profiles
    	self.sbDict = None
        # Build atmosphere and return object
    	self._buildAtmo(parameters, airmass, transmission, aerosolNormCoeff, aerosolNormWavelen)

    def _buildAtmo(self, parameters, airmass, transmission, aerosolNormCoeff, aerosolNormWavelen):
        """Builds an atmospheric transmission profile given a set of component parameters and 
        returns bandpass object. (S^{atm})"""

        parameters = np.array(parameters)

        transDict = {}

        transDict['H2O'] = transmission[airmass]['H2O']**parameters[0]
        transDict['O2'] = transmission[airmass]['O2']**parameters[1]
        transDict['O3'] = transmission[airmass]['O3']**parameters[2]
        transDict['Rayleigh'] = transmission[airmass]['Rayleigh']**parameters[3]
        transDict['Aerosol'] = self._aerosol(self.wavelen, airmass, parameters[5], aerosolNormCoeff, aerosolNormWavelen)**parameters[4]
        totalTrans = transDict['H2O']*transDict['O2']*transDict['O3']*transDict['Rayleigh']*transDict['Aerosol']

        atmo = Bandpass(wavelen=self.wavelen,sb=totalTrans)

        self.sb = atmo.sb
        self.sbDict = transDict
        
        return

    def _aerosol(self, wavelength, airmass, alpha, aerosolNormCoeff, aerosolNormWavelen):
        """Standard aerosol transmission function, returns array of transmission values over a range of
            wavelengths."""
        return np.e**(-aerosolNormCoeff * airmass * (aerosolNormWavelen / wavelength) ** alpha)