import numpy as np
import os
import copy
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass

# Global wavelength variables set to MODTRAN defaults
MINWAVELEN = 300
MAXWAVELEN = 1100
WAVELENSTEP = 0.5

STDPARAMETERS = [1.0,1.0,1.0,1.0,1.0,1.7]
STDAIRMASS = 1.2
STDAEROSOLNORMCOEFF = 0.1
STDAEROSOLNORMWAVELEN = 550.0
STDAEROSOLALPHA = STDPARAMETERS[5]

FILTERLIST = ['u','g','r','i','z','y4']

class Atmo:
    def __init__(self, P, X, aerosolNormCoeff=STDAEROSOLNORMCOEFF, aerosolNormWavelen=STDAEROSOLNORMWAVELEN):
    	self.__generateAtmo(P, X, aerosolNormCoeff=STDAEROSOLNORMCOEFF, aerosolNormWavelen=STDAEROSOLNORMWAVELEN)
    	self.airmass = X
    	self.parameters = np.array(P)
    	self.transmission = None
    	self.components = None
    
    def __generateAtmo(self, P, X, aerosolNormCoeff=STDAEROSOLNORMCOEFF, aerosolNormWavelen=STDAEROSOLNORMWAVELEN):
        """Builds an atmospheric transmission profile given a set of component parameters and 
        returns bandpass object. (S^{atm})"""
        
        self.__parameterCheck(P)
        self.__airmassCheck(X)
        
        H2Ocomp = self.atmoTrans[X]['H2O']**P[0]
        O2comp = self.atmoTrans[X]['O2']**P[1]
        O3comp = self.atmoTrans[X]['O3']**P[2]
        rayleighComp = self.atmoTrans[X]['Rayleigh']**P[3]
        aerosolComp = self.aerosol(self.wavelength,X,P[5],aerosolNormCoeff,aerosolNormWavelen)**P[4]
        totalTrans = H2Ocomp*O2comp*O3comp*rayleighComp*aerosolComp
        
        return Bandpass(wavelen=self.wavelength,sb=totalTrans)

   	def __parameterCheck(self, P):
        """Checks if parameter array is of appropriate length."""
        if len(P) != 6:
            raise ValueError('Need 6 parameters to build atmosphere!')
        return

    def __airmassCheck(self, X):
        if self.airmassToString(X) not in self.airmasses:
            raise ValueError('Not a valid airmass, check MODTRAN data for valid airmasses')
        return