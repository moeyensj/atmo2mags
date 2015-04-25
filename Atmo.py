# Python packages
import numpy as np
import matplotlib.pyplot as plt

# LSST stack software
import lsst.sims.photUtils.Bandpass as Bandpass

class Atmo():
    def __init__(self, parameters, airmass, wavelength, transmission, aerosolNormCoeff, aerosolNormWavelen):
    	self.airmass = airmass
    	self.X = airmass
    	self.parameters = parameters
    	self.P = parameters
    	self.wavelength = wavelength
    	self.aerosolNormCoeff = aerosolNormCoeff
    	self.aerosolNormWavelen = aerosolNormWavelen
    	self.components = ['H2O','O2','O3','Rayleigh','Aerosol']
    	self.componentColors = {'H2O':'blue','O2':'green','O3':'red','Rayleigh':'purple','Aerosol':'cyan'}
    	self.componentsPlot = {'H2O':r'$t_{H_2O}$','O2':r'$t_{O_2}$','O3':r'$t_{O_3}$','Rayleigh':r'$t_{Rayleigh}$','Aerosol':r'$t_{Aerosol}$','Alpha':r'$\alpha$'}

    	self.transmission = None
    	self.sb = None
    	self.transmissionDict = None

    	self.__buildAtmo(parameters, airmass, transmission, aerosolNormCoeff, aerosolNormWavelen)

    def __buildAtmo(self, parameters, airmass, transmission, aerosolNormCoeff, aerosolNormWavelen):
        """Builds an atmospheric transmission profile given a set of component parameters and 
        returns bandpass object. (S^{atm})"""

        parameters = np.array(parameters)

        transDict = {}

        transDict['H2O'] = transmission[airmass]['H2O']**parameters[0]
        transDict['O2'] = transmission[airmass]['O2']**parameters[1]
        transDict['O3'] = transmission[airmass]['O3']**parameters[2]
        transDict['Rayleigh'] = transmission[airmass]['Rayleigh']**parameters[3]
        transDict['Aerosol'] = self.__aerosol(self.wavelength, airmass, parameters[5], aerosolNormCoeff, aerosolNormWavelen)**parameters[4]
        totalTrans = transDict['H2O']*transDict['O2']*transDict['O3']*transDict['Rayleigh']*transDict['Aerosol']

        atmo = Bandpass(wavelen=self.wavelength,sb=totalTrans)

        self.transmission = atmo.sb
        self.sb = atmo.sb
        self.transmissionDict = transDict
        
        return

    def __aerosol(self, wavelength, airmass, alpha, aerosolNormCoeff, aerosolNormWavelen):
        """Standard aerosol transmission function, returns array of transmission values over a range of
            wavelengths."""
        return np.e**(-aerosolNormCoeff * airmass * (aerosolNormWavelen / wavelength) * alpha)