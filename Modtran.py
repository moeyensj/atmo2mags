import os
import copy
import numpy as np

# Global wavelength variables set to MODTRAN defaults
MINWAVELEN = 300
MAXWAVELEN = 1100
WAVELENSTEP = 0.5

STDPARAMETERS = [1.0,1.0,1.0,1.0,1.0,1.7]
STDAIRMASS = 1.2
STDAEROSOLNORMCOEFF = 0.1
STDAEROSOLNORMWAVELEN = 550.0
STDAEROSOLALPHA = STDPARAMETERS[5]

class Modtran():
    def __init__(self, modtranDir='modtran/', modtranRoot='Pachon_MODTRAN', modtranSuffix='.7sc'):
        # List of wavelengths
    	self.wavelengths = None
        # List of airmasses for which we have profiles, set in readModtranFiles
        self.airmasses = None
        # List of transmission profiles for individual airmasses
        self.transmission = None
        # List of strings containing components names to be used when plotting
        self.components = ['H2O','O2','O3','Rayleigh','Aerosol']
        # Read in MODTRAN files
    	self.__readModtranFiles()
    	
    def __readModtranFiles(self, modtranDir='modtran/', modtranRoot='Pachon_MODTRAN', modtranSuffix='.7sc'):
        """
        Reads atmospheric absorption data into an airmass-keyed directory from MODTRAN files.

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        modtranDir: (string) ['modtran/'], MODTRAN file directory
        modtranRoot: (string) ['Pachon_Modtran'], MODTRAN file root
        modtranSuffix: (string) ['.7sc'], MODTRAN file extension
        ----------------------

        * Modified from AtmoComp.py *
        """
        
        files = os.listdir(modtranDir)
        modtranFiles = []
        
        for f in files:
            if (f.startswith(modtranRoot)) & (f.endswith(modtranSuffix)):
                modtranFiles.append(f)
        
        if len(modtranFiles) > 0:
            print "Found " + str(len(modtranFiles)) + " MODTRAN files:"
        
        self.wavelengths = np.arange(MINWAVELEN, MAXWAVELEN+WAVELENSTEP, WAVELENSTEP, dtype='float')
        self.transmission = {}
        self.airmasses = []
        
        for file in modtranFiles:
            print file
            fin = open(os.path.join(modtranDir, file),'r')
            wavelenTemp = []
            transTemp = {}
            
            for comp in self.components:
                transTemp[comp] = []
            
            for lineNum,line in enumerate(fin):
                if lineNum < 4:
                    if lineNum == 1:
                        lineEle = line.split()
                        airmass = 1/np.cos((float(lineEle[2]))*np.pi/180.0)
                        airmass = round(airmass*10)/10.0
                    continue
                lineEle = line.split()
                if (float(lineEle[0]) < 0):
                    break
                if (float(lineEle[0]) > MINWAVELEN) | (float(lineEle[0]) > MAXWAVELEN):
                    wavelenTemp.append(float(lineEle[0]))
                    transTemp['H2O'].append(float(lineEle[2]))
                    transTemp['O2'].append(float(lineEle[3]))
                    transTemp['O3'].append(float(lineEle[4]))
                    transTemp['Rayleigh'].append(float(lineEle[8]))
                    transTemp['Aerosol'].append(self.__aerosol(float(lineEle[0]), airmass, STDAEROSOLALPHA, STDAEROSOLNORMCOEFF, STDAEROSOLNORMWAVELEN))
            fin.close()
            wavelenTemp = np.array(wavelenTemp, dtype='float')
            trans = {}
            templates = {}
            for comp in self.components:
                trans[comp] = np.array(transTemp[comp], dtype='float')
                trans[comp] = np.interp(self.wavelengths, wavelenTemp, transTemp[comp], left=0.0, right=0.0)
            self.airmasses.append(self.__airmassToString(airmass))
            self.transmission[airmass] = copy.deepcopy(trans)
        if self.transmission != None:
            print "MODTRAN files have been read."
        return

    def __aerosol(self, w, X, alpha, aerosolNormCoeff, aerosolNormWavelen):
        """Standard aerosol transmission function, returns array of transmission values over a range of
            wavelengths."""
        return np.e**(-aerosolNormCoeff * X * (aerosolNormWavelen / w) * alpha)

    def __airmassToString(self, airmass):
        """Converts airmass to string"""
        X = float(airmass)
        return "%.3f" % (X)
