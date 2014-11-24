import numpy
import pylab
import os
import copy

class atmo:   
    def __init__(self):
        # List of strings containing component names
        self.components = ['H2O','O2','O3','Rayleigh','Aerosol']
        # List of strings containing components names to be used when plotting
        self.componentsPlot = ['$H_2O$','$O_2$','$O_3$','Rayleigh','Aerosol']
        # List of parameters (H2O,O2,O3,Rayleigh,Aerosol,Alpha). 
        self.parameters = [1.0,1.0,1.0,1.0,1.0,1.7]
        # List of parameters used for plotting
        self.parametersPlot = ['$t_{H_2O}$','$t_{O_2}$','$t_{O_3}$','$t_{Rayleigh}$','$t_{Aerosol}$','$alpha$']
        self.componentsColor = ['blue','green','red','purple','cyan']
        self.wavelength = None
        self.airmasses = None
        self.atmoTrans = None  
        
    def readModtranFiles(self, modtranDir='.', modtranRoot='Pachon_MODTRAN',modtranSuffix='.7sc'):
        files = os.listdir('.')
        modtranFiles = []
        
        for f in files:
            if (f.startswith(modtranRoot)) & (f.endswith(modtranSuffix)):
                modtranFiles.append(f)
        
        if len(modtranFiles) > 0:
            print "Found " + str(len(modtranFiles)) + " MODTRAN files:"
        
        self.wavelength = numpy.arange(250,1100.5,0.5,dtype='float')
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
            self.airmasses.append(self.airmassToString(airmass))
            self.atmoTrans[airmass] = copy.deepcopy(trans)
        if self.atmoTrans != None:
            print "MODTRAN files have been read."
        return
    
    def genAtmo(self,w,P,X=1.0,aerosolNormCoeff=0.1):
        H2Ocomp = self.atmoTrans[X]['H2O']**P[0]
        O2comp = self.atmoTrans[X]['O2']**P[1]
        O3comp = self.atmoTrans[X]['O3']**(X*P[2])
        rayleighComp = self.atmoTrans[X]['Rayleigh']**(X*P[3])
        aerosolComp = self.aerosol(w,X,alpha=P[5],aerosolNormCoeff=aerosolNormCoeff)**P[4]
        return H2Ocomp*O2comp*O3comp*rayleighComp*aerosolComp
    
    def genAtmoPlot(self,P,X=1.0,aerosolNormCoeff=0.1,wavelengthRange=[300,1100],includeStdAtmo=True,
                    stdAtmoAirmass=1.0,genAtmoColor='blue',stdAtmoColor='black',stdAtmoColorAlpha=0.5,
                    stdAtmoParameters=[1.0,1.0,1.0,1.0,1.0,1.7],stdAerosolNormCoeff=0.1):
        
        w=self.wavelength
        fig, ax = pylab.subplots(1,1)
        fig.set_size_inches(12,6) 

        ax.plot(w,self.genAtmo(w,P,X,aerosolNormCoeff),color=genAtmoColor,label=self.labelGen(P));
        ax.set_xlabel("Wavelength, $\lambda$ (nm)")
        ax.set_ylabel("Transmission")
        ax.set_title("$S^{atm}(\lambda)$ and $S^{atm}_{std}(\lambda)$");
        ax.legend(loc='lower right',shadow=False) 
        ax.set_xlim(wavelengthRange[0],wavelengthRange[1]);
        
        if includeStdAtmo == True:
            stdAtmoParams = [1.0,1.0,1.0,1.0,1.0,1.7]
            ax.plot(w,self.genAtmo(w,stdAtmoParams,X=stdAtmoAirmass,aerosolNormCoeff=stdAerosolNormCoeff),
                        label=self.labelGen(stdAtmoParams),alpha=stdAtmoColorAlpha,color=stdAtmoColor);
                
        ax.legend(loc='lower right',shadow=False) 
        return
        
    def aerosol(self,w,X,alpha=1.7,aerosolNormCoeff=0.1):
        return numpy.e**(-aerosolNormCoeff*X*(550.0/w)*alpha)
    
    def airmassToString(self,airmass):
        X = float(airmass)
        return "%.3f" % (X)
        
    def labelGen(self,P):
        label = []
        for paramNum,param in enumerate(P):
            name = self.parametersPlot[paramNum] + ':'
            value = P
            labelEle = name + str(param)
            label.append(labelEle)
        return ' '.join(label)
