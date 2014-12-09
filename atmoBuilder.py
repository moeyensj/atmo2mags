# Necessary imports
import numpy
import pylab
import os
import copy
from AtmoCompMod import AtmoComp as ac
import plot_dmagsMod as dm

# Global wavelength variables set to MODTRAN defaults
WMIN = 300
WMAX = 1100.5
WSTEP = 0.5

"""
    #### IMPORTANT NOTE ####
    WMIN,WMAX,WSTEP must also be set to the above in Sed.py,Bandpass.py and plot_dmagsMod.py or else the wavelengths will not
    be gridded properly and array multiplication errors will occur.
    
    The limiting factor is the MODTRAN data from which we build the standard atmosphere profile used to generate all subsequent
    atmospheres.
    
    """

class atmoBuilder:
    def __init__(self):
        # List of strings containing component names
        self.components = ['H2O','O2','O3','Rayleigh','Aerosol']
        # List of strings containing components names to be used when plotting
        self.componentsPlot = ['$H_2O$','$O_2$','$O_3$','Rayleigh','Aerosol']
        # List of parameters (H2O,O2,O3,Rayleigh,Aerosol,Alpha).
        self.parameters = [1.0,1.0,1.0,1.0,1.0,1.7]
        # List of parameters used for plotting
        self.parametersPlot = ['$t_{H_2O}$','$t_{O_2}$','$t_{O_3}$','$t_{Rayleigh}$','$t_{Aerosol}$','$alpha$']
        # List of colors for used in plotting individual absorption components
        self.componentsColor = ['blue','green','red','purple','cyan']
        # Effective wavelength range, set in readModtranFiles
        self.wavelength = None
        # Min, max values of wavelength range
        self.wavelengthRange = [WMIN,WMAX]
        # List of airmasses for which we have profiles, set in readModtranFiles
        self.airmasses = None
        # List of transmission profiles for individual airmasses
        self.atmoTrans = None
        # Filter-keyed dictionary of filter S, set in readFilters
        self.filters = None
        # Filter-keyed dictionary of filter and hardware
        self.sys = None
        self.readModtranFiles()
        self.readFilters()
        self.readHardware()
    
    def readModtranFiles(self, modtranDir='.', modtranRoot='Pachon_MODTRAN',modtranSuffix='.7sc'):
        """Reads in atmospheric absorption data from MODTRAN files into an airmass-keyed directory."""
        ### Taken from AtmoComp and modified to suit specific needs.
        
        files = os.listdir('.')
        modtranFiles = []
        
        for f in files:
            if (f.startswith(modtranRoot)) & (f.endswith(modtranSuffix)):
                modtranFiles.append(f)
        
        if len(modtranFiles) > 0:
            print "Found " + str(len(modtranFiles)) + " MODTRAN files:"
        
        self.wavelength = numpy.arange(WMIN,WMAX,WSTEP,dtype='float')
        self.atmoTemplates = {}
        self.atmoTrans = {}
        self.airmasses = []
        
        for file in modtranFiles:
            print file
            fin = open(file,'r')
            wavelenTemp = []
            transTemp = {}
            
            for comp in self.components:
                transTemp[comp] = []
            
            for lineNum,line in enumerate(fin):
                if lineNum < 4:
                    if lineNum == 1:
                        lineEle = line.split()
                        airmass = 1/numpy.cos((float(lineEle[2]))*numpy.pi/180.0)
                        airmass = round(airmass*10)/10.0
                    continue
                lineEle = line.split()
                if (float(lineEle[0]) < 0):
                    break
                if (float(lineEle[0]) > self.wavelengthRange[0]) | (float(lineEle[0]) > self.wavelengthRange[1]):
                    wavelenTemp.append(float(lineEle[0]))
                    transTemp['H2O'].append(float(lineEle[2]))
                    transTemp['O2'].append(float(lineEle[3]))
                    transTemp['O3'].append(float(lineEle[4]))
                    transTemp['Rayleigh'].append(float(lineEle[8]))
                    transTemp['Aerosol'].append(self.aerosol(float(lineEle[0]),airmass))
            fin.close()
            wavelenTemp = numpy.array(wavelenTemp,dtype='float')
            trans = {}
            templates = {}
            for comp in self.components:
                trans[comp] = numpy.array(transTemp[comp],dtype='float')
                trans[comp] = numpy.interp(self.wavelength, wavelenTemp, transTemp[comp], left=0.0, right=0.0)
            self.airmasses.append(self.airmassToString(airmass))
            self.atmoTrans[airmass] = copy.deepcopy(trans)
        if self.atmoTrans != None:
            print "MODTRAN files have been read."
        return
    
    def readHardware(self):
        """Reads in the hardware data."""
        sys = dm.read_hardware()
        self.sys = sys
        return sys
    
    def readFilters(self):
        """Reads in the filer data only."""
        filters = dm.read_filtersonly()
        self.filters = filters
        return filters
    
    def genAtmo(self,P,X=1.0,aerosolNormCoeff=0.1):
        """Builds an atmospheric transmission profile given a set of component parameters."""
        H2Ocomp = self.atmoTrans[X]['H2O']**P[0]
        O2comp = self.atmoTrans[X]['O2']**P[1]
        O3comp = self.atmoTrans[X]['O3']**(X*P[2])
        rayleighComp = self.atmoTrans[X]['Rayleigh']**(X*P[3])
        aerosolComp = self.aerosol(self.wavelength,X,alpha=P[5],aerosolNormCoeff=aerosolNormCoeff)**P[4]
        return H2Ocomp*O2comp*O3comp*rayleighComp*aerosolComp
    
    def phi(self,atmo,sys=None):
        """Calculates the normalized bandpass response function for a given sys and atmo, returns a
            filter-keyed dictionary of phi"""
        if sys == None:
            sys = self.sys
        phi = {}
        newSys = {}
        for filter in sys:
            newSys[filter] = atmo*sys[filter].sb
            
            newSys[filter]= newSys[filter]/sys[filter].wavelen
            norm = numpy.sum(newSys[filter])*WSTEP
            phi[filter] = newSys[filter]/norm
        
        return phi
    
    
    def dPhi(self,phi1,phi2):
        """Returns a filter-keyed dictionary of delta phi values"""
        phi = {}
        for p in phi1:
            phi[p]=phi1[p]-phi2[p]
        return phi
    
    ### Plotting functions
    
    def transPlot(self,P,X=1.0,aerosolNormCoeff=0.1,wavelengthRange=[WMIN,WMAX],includeStdAtmo=True,
                  stdAtmoAirmass=1.0,genAtmoColor='blue',stdAtmoColor='black',stdAtmoColorAlpha=0.5,
                  stdAtmoParameters=[1.0,1.0,1.0,1.0,1.0,1.7],stdAerosolNormCoeff=0.1):
        
        w=self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(12,6)
        
        ax.plot(w,self.genAtmo(P,X=X,aerosolNormCoeff=aerosolNormCoeff),color=genAtmoColor,label=self.labelGen(P));
        ax.set_xlabel("Wavelength, $\lambda$ (nm)")
        ax.set_ylabel("Transmission")
        ax.set_title("$S^{atm}(\lambda)$ and $S^{atm,std}(\lambda)$");
        ax.legend(loc='lower right',shadow=False)
        ax.set_xlim(wavelengthRange[0],wavelengthRange[1]);
        
        if includeStdAtmo == True:
            stdAtmoParams = [1.0,1.0,1.0,1.0,1.0,1.7]
            ax.plot(w,self.genAtmo(stdAtmoParams,X=stdAtmoAirmass,aerosolNormCoeff=stdAerosolNormCoeff),
                    label=self.labelGen(stdAtmoParams),alpha=stdAtmoColorAlpha,color=stdAtmoColor);
        
        ax.legend(loc='lower right',shadow=False)
        return

    def filterPlot(self,plotWidth=12,plotHeight=6):
        """Plots the filter response curve from LSST filter data."""
        if self.filters == None:
            self.readFilters()
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        for f in self.filters:
            ax.plot(self.filters[f].wavelen,self.filters[f].sb,label=str(f));
        
        ax.set_xlim(300,1100);
        ax.set_ylim(0,1);
        ax.set_ylabel("Transmission");
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("LSST $S^{sys}_b$ Filters Only");
        ax.legend(loc=4,shadow=False);
        return
    
    def hardwarePlot(self,plotWidth=12,plotHeight=6):
        """Plots the hardware response curve from LSST hardware data."""
        if self.sys == None:
            self.readHardware()
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        for f in self.sys:
            ax.plot(self.sys[f].wavelen,self.sys[f].sb,label=str(f));
        
        ax.set_xlim(300,1100);
        ax.set_ylim(0,1);
        ax.set_ylabel("Transmission");
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("LSST $S^{sys}_b$ Filters and Hardware");
        ax.legend(loc=4,shadow=False);
        return
    
    def phiPlot(self,phi1,phi2=None,plotWidth=12,plotHeight=6,phi2Alpha=0.5,phi2Color='black'):
        """Plots normalized bandpass response function, with the possibility to add a second function
            for comparison."""
        w = self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        if phi2 != None:
            phi2 = self.phi(self.genAtmo([1.0,1.0,1.0,1.0,1.0,1.7]))
        
        for p in phi1:
            ax.plot(w,phi1[p])
            ax.plot(w,phi2[p],alpha=phi2Alpha,color=phi2Color)
        
        ax.set_xlim(300,1100);
        ax.set_ylabel("$\phi_b$");
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("Bandpass Normalization Response");
        
        return
    
    def dPhiPlot(self,phi1,phi2,plotWidth=12,plotHeight=6):
        """Plots delta normalized bandpass response function."""
        
        w = self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        phi = self.dPhi(phi1,phi2)
        
        for p in phi:
            ax.plot(w,phi[p],label=str(p))
        
        ax.set_xlim(300,1100);
        ax.set_ylabel("$\Phi");
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("$\Delta\Phi$");
        ax.legend(loc=4,shadow=False)
        
        return

    def allPlot(self,P,X=1.0,aerosolNormCoeff=0.1,transPlot=True,phiPlot=True,dPhiPlot=True):
        """Generates an atmo and plots a bunch of functions."""
        atmo = self.genAtmo(P,X,aerosolNormCoeff)
        phi = self.phi(atmo)
        
        atmoStd = self.genAtmo([1.0,1.0,1.0,1.0,1.0,1.7])
        phiStd = self.phi(atmoStd)
        
        if transPlot:
            self.transPlot(P)
        if phiPlot:
            self.phiPlot(phi,phi2=phiStd)
        if dPhiPlot:
            self.dPhiPlot(phi,phiStd)
        
        return
    
    ### Secondary Functions
    
    def aerosol(self,w,X,alpha=1.7,aerosolNormCoeff=0.1):
        """Standard aerosol transmission function, returns array of transmission values over a range of
            wavelenghts"""
        return numpy.e**(-aerosolNormCoeff*X*(550.0/w)*alpha)
    
    def airmassToString(self,airmass):
        """Converts airmass to string"""
        X = float(airmass)
        return "%.3f" % (X)
    
    def labelGen(self,P):
        """Generates label for use in plot legends."""
        label = []
        for paramNum,param in enumerate(P):
            name = self.parametersPlot[paramNum] + ':'
            value = P
            labelEle = name + str(param)
            label.append(labelEle)
        return ' '.join(label)