# Necessary imports
import numpy
import pylab
import os
import copy
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass
import lsst.sims.photUtils.photUtils as photUtils

# Global wavelength variables set to MODTRAN defaults
MINWAVELEN = 300
MAXWAVELEN = 1100
WAVELENSTEP = 0.5

"""
#### IMPORTANT NOTE ####
MINWAVELEN,MAXWAVELEN,WAVELENSTEP must also be set to the above in Sed.py,Bandpass.py and plot_dmagsMod.py or else the wavelengths will not
be gridded properly and array multiplication errors will occur.
    
The limiting factor is the MODTRAN data from which we build the standard atmosphere profile used to generate all subsequent
atmospheres.
    
"""

class AtmoBuilder:
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
        self.wavelengthRange = [MINWAVELEN,MAXWAVELEN]
        # List of airmasses for which we have profiles, set in readModtranFiles
        self.airmasses = None
        # List of transmission profiles for individual airmasses
        self.atmoTrans = None
        # Filter-keyed dictionary of filter S, set in readFilters
        self.filters = None
        # Filter-keyed dictionary of filter and hardware
        self.sys = None
        # List of filters
        self.filterlist = ['u', 'g', 'r', 'i', 'z', 'y4']
        
        # Kurucz model data
        self.stars = None
        self.starlist = None
        self.temperature = None
        self.met = None
        self.logg = None
        
        # Readers
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
        
        self.wavelength = numpy.arange(MINWAVELEN,MAXWAVELEN+WAVELENSTEP,WAVELENSTEP,dtype='float')
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
    
    def readFilters(self,shift_perc=None):
        """Reads LSST filter data only and returns a filter-keyed dictionary. (S^{filters})"""
        ### Taken from plot_dmags and modified to suit specific needs.
        # read the filter throughput curves only (called from read_hardware as well)
        # apply a shift of +shift_perc/100 * eff_wavelength to the wavelengths of the filter.
        filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
        filters = {}
        for f in self.filterlist:
            filters[f] = Bandpass()
            filters[f].readThroughput(os.path.join(filterdir, "filter_" + f + ".dat"))
            effwavelenphi, effwavelensb = filters[f].calcEffWavelen()
            if shift_perc != None:
                shift = effwavelensb * shift_perc/100.0
                print f, shift
                filters[f].wavelen = filters[f].wavelen + shift
                filters[f].resampleBandpass()
        self.filters = filters
        return
    
    def readHardware(self,shift_perc=None):
        """Reads LSST hardware data and returns a filter-keyed dictionary. (S^{sys})"""
        ### Taken from plot_dmags and modified to suit specific needs.
        # read system (hardware) transmission, return dictionary of system hardware (keyed to filter)
        filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
        hardware = ("detector.dat", "m1.dat", "m2.dat", "m3.dat", "lens1.dat", "lens2.dat", "lens3.dat")
        # Read in the standard components, but potentially shift the filter by shift_perc percent.
        self.readFilters(shift_perc=shift_perc)
        filters = self.filters
        sys = {}
        for f in self.filterlist:
            sys[f] = Bandpass()
            # put together the standard component list
            tlist = []
            for t in hardware:
                tlist.append(os.path.join(filterdir, t))
            # read in the standard components, combine into sys
            sys[f].readThroughputList(tlist)
            # multiply by the filter throughput for final hardware throughput (no atmosphere)
            sys[f].wavelen, sys[f].sb = sys[f].multiplyThroughputs(filters[f].wavelen, filters[f].sb)
        self.sys = sys
        return

    def readKurucz(self):
        """Reads Kurucz model data from LSST software stack and sets relevant class attributes."""
        ### Taken from plot_dmags and modified to suit specific needs.
        # read kurucz model MS, g40 stars SEDs
        homedir = os.getenv("SIMS_SED_LIBRARY_DIR")    # "SIMS_SED_LIBRARY_DIR"
        stardir = os.path.join(homedir, "starSED/kurucz/")  # "starSED/kurucz/"
        allfilelist = os.listdir(stardir)
        starlist = []
        # make preliminary cut for ms, g40 stars
        for filename in allfilelist:
            if filename[-11:-8] == 'g40':
                starlist.append(filename)
            if filename[-11:-8] == 'g20':
                starlist.append(filename)
        atemperature = []
        amet = []
        alogg = []
        starlist2 = []
        # make secondary cut for stars with temperature > 4000 deg
        for s in starlist:
            tmp = s.split('_')
            met = float(tmp[0][2:])
            if tmp[0][1] == 'm':
                met = -1 * met
            met = met/10.0
            temperature = float(tmp[1][:5])
            logg = float(tmp[2][1:])
            logg = logg / 10.0
            if (temperature > 4000.0):
                amet.append(met)
                atemperature.append(temperature)
                alogg.append(logg)
                starlist2.append(s)
        temperature = numpy.array(atemperature)
        met = numpy.array(amet)
        logg = numpy.array(alogg)
        starlist = starlist2
        # actually read the stars SEDS from disk
        stars = {}
        for s in starlist:
            stars[s] = Sed()
            stars[s].readSED_flambda(os.path.join(stardir,s))
        print "# Read %d MS stars from %s" %(len(starlist), stardir)
        # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
        for s in starlist:
            stars[s].synchronizeSED(wavelen_min=WMIN, wavelen_max=WMAX, wavelen_step=WSTEP)

        self.stars = stars
        self.starlist = starlist
        self.temperature = temperature
        self.met = met
        self.logg = logg

        return
    
    def genAtmo(self,P,X=1.0,aerosolNormCoeff=0.1,aerosolNormWavelength=550.0):
        """Builds an atmospheric transmission profile given a set of component parameters and 
        returns bandpass object. (S^{atm})"""
        self.parameterCheck(P)
        H2Ocomp = self.atmoTrans[X]['H2O']**P[0]
        O2comp = self.atmoTrans[X]['O2']**P[1]
        O3comp = self.atmoTrans[X]['O3']**P[2]   # linear
        rayleighComp = self.atmoTrans[X]['Rayleigh']**P[3]  # linear
        aerosolComp = self.aerosol(self.wavelength,X,alpha=P[5],aerosolNormCoeff=aerosolNormCoeff,aerosolNormWavelength=aerosolNormWavelength)**P[4]
        totalTrans = H2Ocomp*O2comp*O3comp*rayleighComp*aerosolComp
        return Bandpass(wavelen=self.wavelength,sb=totalTrans)
    
    def combineThroughputs(self,atmos,sys=None):
        """Combines atmospheric transmission profile with system responsiveness data, returns filter-keyed 
        dictionary. (S^{atm}*S^{sys})"""
        ### Taken from plot_dmags and modified to suit specific needs.
        # Set up the total throughput for this system bandpass
        if sys == None:
            sys = self.sys
        total = {}
        for f in self.filterlist:
            wavelen, sb = sys[f].multiplyThroughputs(atmos.wavelen, atmos.sb)
            total[f] = Bandpass(wavelen, sb)
            total[f].sbTophi()
        return total
    
    def phi(self,atmo,sys=None):
        """Calculates the normalized bandpass response function for a given sys and atmo, returns a
            filter-keyed dictionary of phi."""
        if sys == None:
            sys = self.sys
        phi = {}
        newSys = {}
        for filter in sys:
            newSys[filter] = atmo.sb*sys[filter].sb
            
            newSys[filter] = newSys[filter]/sys[filter].wavelen
            norm = numpy.sum(newSys[filter])*WSTEP
            phi[filter] = newSys[filter]/norm
        
        return phi
    
    def dPhi(self,phi1,phi2):
        """Returns a filter-keyed dictionary of delta phi values."""
        dphi = {}
        for p in phi1:
            dphi[p] = phi1[p] - phi2[p]
        return dphi
    
    
    def mags(self,bpDict):
        """Calculates magnitudes given a bandpass dictionary, returns filter-keyed magnitude dictionary."""
        ### Taken from plot_dmags and modified to suit specific needs.
        # calculate magnitudes for all sed objects using bpDict (a single bandpass dictionary keyed on filters)
        # pass the sedkeylist so you know what order the magnitudes are arranged in
        
        self.kuruczCheck()
        
        seds = self.stars
        sedkeylist = self.starlist
        filterlist = self.filterlist
        
        mags = {}
        for f in self.filterlist:
            mags[f] = numpy.zeros(len(sedkeylist), dtype='float')
            i = 0
            for key in sedkeylist:
                mags[f][i] = seds[key].calcMag(bpDict[f])
                if numpy.isnan(mags[f][i]):
                    print key, f, mags[f][i]
                i = i + 1
            print f, mags[f].max(), mags[f].min()
        return mags
    
    def dmags(self,mags1,mags2):
        """Returns filter-keyed dictionary of change in magnitude in millimagnitudes."""
        ### Taken from plot_dmags and modified to suit specific needs.
        dmags = {}
        for f in self.filterlist:
            # difference, in millimags
            dmags[f] = (mags1[f] - mags2[f]) * 1000.0
        return dmags
    
    def gi(self,mags_std):
        """Returns standard color temperature given standard magnitude dictionary keyed on filters."""
        ### Taken from plot_dmags and modified to suit specific needs.
        # calculate some colors in the standard atmosphere, should be also standard bandpass, not shifted)
        gi = mags_std['g'] - mags_std['i']
        return gi
    
    ### Plotting functions
    
    def transPlot(self,P,X=1.0,aerosolNormCoeff=0.1,aerosolNormWavelength=550.0,wavelengthRange=[WMIN,WMAX],includeStdAtmo=True,
                  stdAtmoAirmass=1.0,genAtmoColor='blue',stdAtmoColor='black',stdAtmoColorAlpha=0.5,
                  stdAtmoParameters=[1.0,1.0,1.0,1.0,1.0,1.7],stdAerosolNormCoeff=0.1,stdAerosolNormWavelength=550.0,figName=None):
        """Plots atmospheric transmission profile given a parameter array."""
        
        w=self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(12,6)
        
        atmo = self.genAtmo(P,X=X,aerosolNormCoeff=aerosolNormCoeff,aerosolNormWavelength=aerosolNormWavelength)
        
        ax.plot(w,atmo.sb,color=genAtmoColor,label=self.labelGen(P,X));
        ax.set_xlabel("Wavelength, $\lambda$ (nm)")
        ax.set_ylabel("Transmission")
        ax.set_title("$S^{atm}(\lambda)$ and $S^{atm,std}(\lambda)$");
        ax.legend(loc='lower right',shadow=False)
        ax.set_xlim(wavelengthRange[0],wavelengthRange[1]);
        
        if includeStdAtmo == True:
            stdAtmoParams = [1.0,1.0,1.0,1.0,1.0,1.7]
            atmoStd = self.genAtmo(stdAtmoParams,X=stdAtmoAirmass,aerosolNormCoeff=stdAerosolNormCoeff,aerosolNormWavelength=stdAerosolNormWavelength)
            ax.plot(w,atmoStd.sb,
                    label=self.labelGen(stdAtmoParams,X=stdAtmoAirmass),alpha=stdAtmoColorAlpha,color=stdAtmoColor);
        
        ax.legend(loc='lower right',shadow=False)
        
        if figName != None:
            title = figName + "_transPlot.png"
            pylab.savefig(title,format='png')
        return
    
    def filterPlot(self,plotWidth=12,plotHeight=6):
        """Plots the filter response curve from LSST filter data."""
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        for f in self.filterlist:
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
    
    def phiPlot(self,phi1,phi2=None,plotWidth=12,plotHeight=6,phi2Alpha=0.5,phi2Color='black',figName=None):
        """Plots normalized bandpass response function, with the possibility to add a second function
            for comparison."""
        w = self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        if phi2 != None:
            phi2 = self.phi(self.genAtmo([1.0,1.0,1.0,1.0,1.0,1.7]))
        
        for f in self.filterlist:
            ax.plot(w,phi1[f],label=str(f))
            ax.plot(w,phi2[f],alpha=phi2Alpha,color=phi2Color)
        
        ax.set_xlim(300,1100);
        ax.set_ylabel("$\phi_b^{obs}(\lambda)$",fontsize=15);
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("Normalized Bandpass Response");
        ax.legend(loc=4,shadow=False);
        
        if figName != None:
            title = figName + "_phiPlot.png"
            pylab.savefig(title,format='png')
        return
    
    def dPhiPlot(self,phi1,phi2,plotWidth=12,plotHeight=6,figName=None):
        """Plots change in normalized bandpass response function given two phi functions."""
        
        w = self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        phi = self.dPhi(phi1,phi2)
        
        for f in self.filterlist:
            ax.plot(w,phi[f],label=str(f))
        
        ax.set_xlim(300,1100);
        ax.set_ylabel("$\Delta\phi_b^{obs}(\lambda)$",fontsize=15);
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("Change in Normalized Bandpass Response");
        ax.legend(loc=4,shadow=False)
        
        if figName != None:
            title = figName + "_dPhiPlot.png"
            pylab.savefig(title,format='png')
        
        return
    
    def colorcolorPlot(self, bpDict, newfig=True, titletext=None):
        
        seds = self.stars
        sedkeylist = self.starlist
        sedcolorkey = self.met
        
        mags = self.mags(bpDict)
        gi = self.gi(mags)
        magcolors = {}
        colorlabels = numpy.chararray(len(self.filterlist))
        for i in range(1-7):
            colorlabels[i] = self.filterlist[i-1] + '-' + self.filterlist[i]
            magscolors[colorlabels[i]] = mags[self.filterlist[i-1]] - mags[self.filterlist[i]]
        # colors = ug, gr, ri, iz, zy
        metallicity = sedcolorkey
        metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
        metbinsize = abs(sedcolorkey.min() - sedcolorkey.max())/6.0
        metbins = numpy.arange(sedcolorkey.min(), sedcolorkey.max() + metbinsize, metbinsize)
        print metbinsize, metbins
        i = 1
        # use a different subplot for each color/color combo
        for i in range(len(colorlabels-1)):
            ax = pylab.subplot(3,2,i)
            for metidx in range(len(metbins)):
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
                mcolor = metcolors[metidx]
                pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], mcolor+'.')
            i = i + 1
        ax = pylab.subplot(3,2,7)
        for metidx in range(len(metbins)):
            condition = ((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
            mcolor = metcolors[metidx]
            pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], mcolor+'.')
        # set up generic items
        for i in range(1, 7):
            f = filterlist[i-1]
            ax = pylab.subplot(3,2,i)
            #pylab.xlabel("g-i")
            #pylab.ylabel(r"$\Delta$ %s (mmag)" %(f))
            def axis_formatter(x, pos):
                return "%.1f" %(x)
            formatter = pylab.FuncFormatter(axis_formatter)
            ax.yaxis.set_major_formatter(formatter)
            # set axes limits
            if ylims == None:
                pass
            else:
                try:
                    pylab.ylim(ylims[f][0], ylims[f][1])
                except KeyError:
                    pass
        # put a grid in the background
        if newfig:
            for i in range(1, 7):
                ax = pylab.subplot(3, 2, i)
                pylab.grid(True)
                #pylab.suptitle(titletext)
                #pylab.savefig("delta_mags2.eps", format='eps')
        return


    def dmagsPlot(self, gi, dmags, titletext=None, ylims=None, xlims=None, newfig=True, figName=None,verbose=False):
        """Plots dmags with each filter in its own subplot."""
        ### Taken from plot_dmags and modified to suit specific needs.
        sedcolorkey = [self.met,self.logg]
        plotfilterlist = self.filterlist
        
        # make figure of change in magnitude
        if newfig:
            pylab.figure(figsize=(10,15))
        yplots = 3
        xplots = 2
        if len(plotfilterlist) == 1:
            yplots = 1
            xplots = 1
        pylab.subplots_adjust(top=0.93, wspace=0.32, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
        # For Kurucz models
        # set colors of data points based on their metallicity
        metallicity = numpy.array(sedcolorkey[0])
        logg = numpy.array(sedcolorkey[1])
        metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
        metbinsize = abs(metallicity.min() - metallicity.max())/6.0
        metbins = numpy.arange(metallicity.min(), metallicity.max() + metbinsize, metbinsize)
        if verbose:
            print metbinsize, metbins
        
        i = 1
        # for each filter, use a different subplot
        plot_logg2 = False
        for f in plotfilterlist:
            ax = pylab.subplot(yplots, xplots,i)
            for metidx in range(len(metbins)):
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize) \
                            & (logg>3.5))
                mcolor = metcolors[metidx]
                pylab.plot(gi[condition], dmags[f][condition], mcolor+'.')
                if plot_logg2:
                        condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize) \
                            & (logg<2.5))
                        mcolor = metcolors[metidx]
                        pylab.plot(gi[condition], dmags[f][condition], mcolor+'x')
            i = i + 1
        # set up generic items
        for i in range(0, len(plotfilterlist)):
            f = plotfilterlist[i]
            ax = pylab.subplot(yplots,xplots,i+1)
            pylab.xlabel("g-i")
            pylab.ylabel(r"$\Delta$ %s (mmag)" %(f))
            def axis_formatter(x, pos):
                return "%.1f" %(x)
            formatter = pylab.FuncFormatter(axis_formatter)
            ax.yaxis.set_major_formatter(formatter)
            # set axes limits
            if ylims == None:
                pass
            else:
                try:
                    pylab.ylim(ylims[f][0], ylims[f][1])
                except KeyError:
                    pass
            if xlims == None:
                pass
            else:
                try:
                    pylab.xlim(xlims[f][0], xlims[f][1])
                except KeyError:
                    pass
        # put a grid in the background
        if newfig:
            for i in range(0, len(plotfilterlist)):
                ax = pylab.subplot(yplots,xplots, i+1)
                pylab.grid(True)
            if titletext!=None:
                pylab.suptitle("$\Delta$mmags for each LSST filter")
                
        if figName != None:
            title = figName+"_dMagsPlot.png"
            pylab.savefig(title, format='png')
        
        return
    
    def allPlot(self,P,X=1.0,aerosolNormCoeff=0.1,aerosolNormWavelength=550.0,transPlot=True,phiPlot=True,dPhiPlot=True,dmagsPlot=True,saveFig=False,figName=None):
        """Generates an atmosphere with given parameters and plots appropriate functions."""
        
        atmo = self.genAtmo(P,X,aerosolNormCoeff,aerosolNormWavelength)
        
        phi = self.phi(atmo)
        atmoStd = self.genAtmo([1.0,1.0,1.0,1.0,1.0,1.7])
        phiStd = self.phi(atmoStd)
        
        if saveFig == True:
            if figName != None:
                figName='X'+str(int(X*10))+'_'+self.pToString(P)+'_'+figName
            else:
                figName='X'+str(int(X*10))+'_'+self.pToString(P)
        else:
            figName = None
        
        if transPlot:
            self.transPlot(P,X=X,aerosolNormCoeff=aerosolNormCoeff,aerosolNormWavelength=aerosolNormWavelength,figName=figName)
        if phiPlot:
            self.phiPlot(phi,phi2=phiStd,figName=figName)
        if dPhiPlot:
            self.dPhiPlot(phi,phiStd,figName=figName)
        if dmagsPlot:
            bp = self.combineThroughputs(atmo)
            bpStd = self.combineThroughputs(atmoStd)
            
            mag = self.mags(bp)
            magStd = self.mags(bpStd)
            
            dmags = self.dmags(mag,magStd)
            gi = self.gi(magStd)
            
            self.dmagsPlot(gi,dmags,figName=figName)
        return

    ### Secondary Functions
    
    def aerosol(self,w,X,alpha=1.7,aerosolNormCoeff=0.1,aerosolNormWavelength=550.0):
        """Standard aerosol transmission function, returns array of transmission values over a range of
            wavelengths."""
        return numpy.e**(-aerosolNormCoeff*X*(aerosolNormWavelength/w)*alpha)
    
    def airmassToString(self,airmass):
        """Converts airmass to string"""
        X = float(airmass)
        return "%.3f" % (X)
    
    def labelGen(self,P,X=1.0):
        """Generates label for use in plot legends."""
        label = []
        for paramNum,param in enumerate(P):
            name = self.parametersPlot[paramNum] + ':'
            labelEle = name + str(param)
            label.append(labelEle)
        return ' '.join(label) + ' $X$:' + str(X)
    
    def parameterCheck(self,P):
        """Checks if parameter array is of appropriate length."""
        if len(P) != 6:
            print "Need 6 parameters to build atmosphere!"
        return
    
    def pToString(self,P):
        """Returns string version of parameter array."""
        stringP = "P"
        for i in P:
            if i < 1.0:
                stringP+="0"+str(int(i*10))
            else:
                stringP+=str(int(i*10))
        return stringP
    
    def kuruczCheck(self):
        """Checks if Kurucz model data has been read in."""
        if self.stars == None:
            print "No Kurucz model data found, please run self.readKurucz()"
        return
