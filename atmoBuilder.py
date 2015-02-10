# Necessary imports
import numpy
import pylab
import os
import copy
import matplotlib.patches as mp
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass
import lsst.sims.photUtils.photUtils as photUtils

from astroML.plotting.mcmc import convert_to_stdev
from astroML.decorators import pickle_results

# Global wavelength variables set to MODTRAN defaults
MINWAVELEN = 300
MAXWAVELEN = 1100
WAVELENSTEP = 0.5

STDPARAMETERS = [1.0,1.0,1.0,1.0,1.0,1.7]
STDAIRMASS = 1.4
STDAEROSOLNORMCOEFF = 0.1
STDAEROSOLNORMWAVELEN = 550.0
STDAEROSOLALPHA = 1.7

"""
#### IMPORTANT NOTE ####
MINWAVELEN,MAXWAVELEN,WAVELENSTEP must also be set to the above in Sed.py and Bandpass.py or else the wavelengths will not
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
        # List of parameters (H2O,O2,O3,Rayleigh,Aerosol,Alpha)
        self.parameters = [1.0,1.0,1.0,1.0,1.0,1.7]
        # List of parameters used for plotting
        self.parametersPlot = [r'$t_{H_2O}$',r'$t_{O_2}$',r'$t_{O_3}$',r'$t_{Rayleigh}$',r'$t_{Aerosol}$',r'$\alpha$']
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
        # List of filter colors
        self.filtercolors = ['blue', 'green', 'red', 'cyan', 'purple', 'yellow']
        
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
    
    def readModtranFiles(self, modtranDir='.', modtranRoot='Pachon_MODTRAN', modtranSuffix='.7sc'):
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
                    transTemp['Aerosol'].append(self.aerosol(float(lineEle[0]),airmass,STDAEROSOLALPHA,STDAEROSOLNORMCOEFF,STDAEROSOLNORMWAVELEN))
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
    
    def readFilters(self, shift_perc=None):
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
    
    def readHardware(self, shift_perc=None):
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
            stars[s].synchronizeSED(wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP)

        self.stars = stars
        self.starlist = starlist
        self.temperature = temperature
        self.met = met
        self.logg = logg

        return
    
    def genAtmo(self, P, X, aerosolNormCoeff=STDAEROSOLNORMCOEFF, aerosolNormWavelen=STDAEROSOLNORMWAVELEN):
        """Builds an atmospheric transmission profile given a set of component parameters and 
        returns bandpass object. (S^{atm})"""
        
        self.parameterCheck(P)
        self.airmassCheck(X)

        P = numpy.array(P)
        
        H2Ocomp = self.atmoTrans[X]['H2O']**P[0]
        O2comp = self.atmoTrans[X]['O2']**P[1]
        O3comp = self.atmoTrans[X]['O3']**P[2]   # linear
        rayleighComp = self.atmoTrans[X]['Rayleigh']**P[3]  # linear
        aerosolComp = self.aerosol(self.wavelength,X,P[5],aerosolNormCoeff,aerosolNormWavelen)**P[4]
        totalTrans = H2Ocomp*O2comp*O3comp*rayleighComp*aerosolComp
        
        return Bandpass(wavelen=self.wavelength,sb=totalTrans)
    
    def combineThroughputs(self, atmo, sys=None):
        """Combines atmospheric transmission profile with system responsiveness data, returns filter-keyed 
        dictionary. (S^{atm}*S^{sys})"""
        ### Taken from plot_dmags and modified to suit specific needs.
        # Set up the total throughput for this system bandpass
        if sys == None:
            sys = self.sys
        total = {}
        for f in self.filterlist:
            wavelen, sb = sys[f].multiplyThroughputs(atmo.wavelen, atmo.sb)
            total[f] = Bandpass(wavelen, sb)
            total[f].sbTophi()
        return total
    
    def mags(self, bpDict, verbose=False):
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
            if verbose == True:
                print f, mags[f].max(), mags[f].min()
        return mags
    
    def dmags(self, mags1, mags2):
        """Returns filter-keyed dictionary of change in magnitude in millimagnitudes."""
        ### Taken from plot_dmags and modified to suit specific needs.
        dmags = {}
        for f in self.filterlist:
            # difference, in millimags
            dmags[f] = (mags1[f] - mags2[f]) * 1000.0
        return dmags
    
    def gi(self, mags_std):
        """Returns standard color temperature given standard magnitude dictionary keyed on filters."""
        ### Taken from plot_dmags and modified to suit specific needs.
        # calculate some colors in the standard atmosphere, should be also standard bandpass, not shifted)
        gi = mags_std['g'] - mags_std['i']
        return gi

    ### Regression Functions

    def compute_logL(self, P, X, err, f, magsStd):
        """Return logL for a given array of parameters P, airmass X, error, a filter and the magnitudes of a standard atmosphere."""
        atmo = self.genAtmo(P,X)
        throughputAtmo = self.combineThroughputs(atmo)
        mags = self.mags(throughputAtmo)
    
        return -numpy.sum(0.5 * ((mags[f] - magsStd[f]) / err) ** 2)

    def compute_mag_color_nonlinear(self, comp1, comp2, P=STDPARAMETERS, X=STDAIRMASS, err=0.02, Nbins=50, generateFig=True, pickleString=None, filters=None, verbose=True):

        range1, pNum1 = self.componentCheck(comp1,Nbins)
        range2, pNum2 = self.componentCheck(comp2,Nbins)
        
        if filters == None:
            filters = self.filterlist

        if verbose:
            print 'Computing nonlinear regression for ' + comp1 + ' and ' + comp2 + '.'
            print 'Arbitrary Atmosphere Parameters = ' + str(P)

        if pickleString != None:
            pickleString = self.pickleNameGen(comp1, comp2, P, X, Nbins) + '_' + pickleString + '.pkl'
        else:
            pickleString = self.pickleNameGen(comp1, comp2, P, X, Nbins) + '.pkl'

        P_std = copy.deepcopy(P)

        std = self.genAtmo(P,X)
        throughputStd = self.combineThroughputs(std)
        magsStd = self.mags(throughputStd)
        giStd = self.gi(magsStd)
                                    
        @pickle_results(pickleString)
        def run_regression(comp1, comp2):
            
            logL = {}
            whr = {}
            comp1best = {}
            comp2best = {}

            for f in filters:
                logL[f] = numpy.empty((Nbins, Nbins))
                for i in range(len(range1)):
                    for j in range(len(range2)):
                        P[pNum1] = range1[i]
                        P[pNum2] = range2[j]
                        logL[f][i, j] = self.compute_logL(P,X,err,f,magsStd)

                print 'Completed ' + f + ' filter.'

            for f in filters:
                logL[f] -= numpy.max(logL[f].T)
                whr[f] = numpy.where(logL[f] == numpy.max(logL[f]))
                comp1best[f] = range1[whr[f][0][0]]
                comp2best[f] = range2[whr[f][1][0]]

            return comp1best, comp2best, logL

        comp1best, comp2best, logL  = run_regression(comp1, comp2)

        if verbose:
            print 'filter, ' + comp1 + ', ' + comp2
            for f in filters:
                print f, comp1best[f], comp2best[f]

        if generateFig == True:
            self.regressionPlot(comp1, comp1best, comp2, comp2best, logL, P_std, X, pNum1=pNum1, pNum2=pNum2,
                                comp1range=range1, comp2range=range2, Nbins=Nbins, figName=pickleString, filters=filters, verbose=verbose)

        return range1, range2, comp1best, comp2best, logL

    def regressionPlot(self, comp1, comp1best, comp2, comp2best, logL, P, X, pNum1=None, pNum2=None,
                       comp1range=None, comp2range=None, Nbins=50, figName=None, filters=None , verbose=True):
        """Plots dmags with each filter in its own subplot."""
        ### Taken from plot_dmags and modified to suit specific needs.
        sedcolorkey = [self.met,self.logg]

        if any([pNum1,pNum2,comp1range,comp2range]) == None:
            comp1range, pNum1 = self.componentCheck(comp1, Nbins)
            comp2range, pNum2 = self.componentCheck(comp2, Nbins)
            
        if filters == None:
            filters = self.filterlist
        
        fig, ax = pylab.subplots(len(filters),3)
        fig.suptitle(r'$\Delta$mmags, Regression Contours and dPhi for each LSST filter', fontsize=14)
        fig.set_size_inches(12,len(filters)*4)
        fig.subplots_adjust(top=0.93, wspace=0.20, hspace=0.20, bottom=0.09, left=0.10, right=0.96)
        
        metallicity = numpy.array(sedcolorkey[0])
        logg = numpy.array(sedcolorkey[1])
        metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
        metbinsize = abs(metallicity.min() - metallicity.max())/6.0
        metbins = numpy.arange(metallicity.min(), metallicity.max() + metbinsize, metbinsize)

        comp1Std = P[pNum1]
        comp2Std = P[pNum2]

        # Create arbitrary standard atmosphere
        std = self.genAtmo(P,X)
        throughputStd = self.combineThroughputs(std)
        magsStd = self.mags(throughputStd)
        gi = self.gi(magsStd)

        # Create atmosphere at mean best fit parameters
        w = self.wavelength
        P1_mean = self.meanParameter(comp1best, filters)
        P2_mean = self.meanParameter(comp2best, filters)
        P_mean = copy.deepcopy(P)
        P_mean[pNum1] = P1_mean
        P_mean[pNum2] = P2_mean
        atmoMean = self.genAtmo(P_mean,X)
        throughputMean = self.combineThroughputs(atmoMean)
        magsMean = self.mags(throughputMean)
        dmagsMean = self.dmags(magsMean,magsStd)

        if verbose == True:
            print 'Mean atmosphere parameter array: ' + str(P_mean)

        # For each filter plot dmags, regression contours and dphis
        for i,f in enumerate(filters):
            # Set component parameters to best fit parameters
            P[pNum1] = comp1best[f]
            P[pNum2] = comp2best[f]

            # Create atmosphere at best fit parameters
            atmo = self.genAtmo(P,X)
            throughputAtmo = self.combineThroughputs(atmo)
            mags = self.mags(throughputAtmo)
            dmags = self.dmags(mags,magsStd)

            # Plot dmmag plots
            # Initialize values to keep track of dmag range
            dmag_max = 0;
            dmag_min = 0;
            dmag_range_max = 0;

            dmag_maxMean = 0;
            dmag_minMean = 0;
            dmag_range_maxMean = 0;

            dmagFit = 0;
            dmagMean = 0;

            for metidx in range(len(metbins)):
                # Make cut of stars
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize) \
                        & (logg>3.5))
                mcolor = metcolors[metidx]

                # Find minimum, max dmag values that fit condition
                minv = numpy.min(dmags[f][condition])
                maxv = numpy.max(dmags[f][condition])

                minvMean = numpy.min(dmagsMean[f][condition])
                maxvMean = numpy.max(dmagsMean[f][condition])

                # Save lowest and highest dmag
                if minv < dmag_min:
                    dmag_min = copy.deepcopy(minv)
                if maxv > dmag_max:
                    dmag_max = copy.deepcopy(maxv)
                if minvMean < dmag_minMean:
                    dmag_minMean = copy.deepcopy(minvMean)
                if maxvMean > dmag_maxMean:
                    dmag_maxMean = copy.deepcopy(maxvMean)
                    
                dmag_range = (dmag_max - dmag_min)/2.0
                dmag_rangeMean = (dmag_maxMean - dmag_minMean)/2.0
                
                if dmag_range_max < dmag_range:
                    dmag_range_max = copy.deepcopy(dmag_range)
                if dmag_range_maxMean < dmag_rangeMean:
                    dmag_range_maxMean = copy.deepcopy(dmag_rangeMean)

                if i == 0 and metidx == 6:
                    ax[i][0].plot(gi[condition], dmags[f][condition], mcolor+'.', label='Fit')
                    ax[i][0].scatter(gi[condition], dmagsMean[f][condition], color='gray', alpha=0.5, label='Mean')
                else:
                    ax[i][0].plot(gi[condition], dmags[f][condition], mcolor+'.')
                    ax[i][0].scatter(gi[condition], dmagsMean[f][condition], color='gray', alpha=0.5)

            # If dmag range exceeds 2.0, plot dashed lines at +-2dmmags
            if dmag_range_max > 2.0 or dmag_range_maxMean > 2.0:
                ax[i][0].axhline(2,color='black',linestyle='--')
                ax[i][0].axhline(-2,color='black',linestyle='--')
            else:
                ax[i][0].set_ylim(-2,2)

            # Label axes and add grid
            ax[i][0].set_xlabel("g-i")
            ax[i][0].set_ylabel(r"$\Delta$ %s (mmag)" %(f))
            ax[i][0].grid()

            # Plot parameter space regression plots
            # Plot contours and true values
            ax[i][1].contour(comp1range, comp2range, convert_to_stdev(logL[f]), levels=(0.683, 0.955, 0.997),colors='k')
            ax[i][1].scatter(P1_mean, P2_mean, color='gray', alpha=0.8, label='Mean')
            ax[i][1].scatter(comp1Std, comp2Std, label='Std')

            # Plot dashed lines at best fit parameters
            ax[i][1].axvline(comp1best[f], color='black', linestyle='--', label='Fit')
            ax[i][1].axhline(comp2best[f], color='black', linestyle='--')

            # Set y-axis, x-axis limits
            ax[i][1].set_xlim(min(comp1range), max(comp1range))
            ax[i][1].set_ylim(min(comp2range), max(comp2range))

            # Label axes and show legend
            ax[i][1].set_xlabel(comp1)
            ax[i][1].set_ylabel(comp2)

            # Plot dphi plots
            # Plot dphis for each filter
            ax[i][2].plot(w, throughputAtmo[f].phi - throughputStd[f].phi, color=self.filtercolors[i], label='Fit - Std')
            ax[i][2].plot(w, throughputMean[f].phi - throughputStd[f].phi, alpha=0.8, color='gray', label='Mean - Std')

            # Format axes and add labels, legend
            ax[i][2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax[i][2].yaxis.set_label_position('right')
            ax[i][2].set_ylabel(r'$\Delta\phi_' + f + r'^{obs}(\lambda)$');
            ax[i][2].set_xlabel('Wavelength, $\lambda$ (nm)');
            ax[i][2].grid()

            if i == 0:
                ax[i][0].legend(loc='upper center', bbox_to_anchor=(0.5,1.25), ncol=2)
                ax[i][1].legend(loc='upper center', bbox_to_anchor=(0.5,1.25), ncol=3)
                ax[i][2].legend(loc='upper center', bbox_to_anchor=(0.5,1.25), ncol=2)

        if figName != None:
            title = figName+"_regressionPlot.png"
            pylab.savefig(title, format='png')

        return
    
    ### Plotting Functions
    
    def transPlot(self, P1, X1, P2=None, X2=None, includeStdAtmo=True, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN],
                  aerosolNormCoeff1=STDAEROSOLNORMCOEFF, aerosolNormWavelen1=STDAEROSOLNORMWAVELEN,
                  aerosolNormCoeff2=STDAEROSOLNORMCOEFF, aerosolNormWavelen2=STDAEROSOLNORMWAVELEN,
                  atmoColor1='blue', atmo2Color='black', atmo2Alpha=0.5, figName=None):
        """Plots atmospheric transmission profile given a parameter array."""
        
        w=self.wavelength

        atmo1 = self.genAtmo(P1, X1, aerosolNormCoeff=aerosolNormCoeff1, aerosolNormWavelen=aerosolNormWavelen1)
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth,plotHeight)
        
        ax.plot(w, atmo1.sb, color=atmoColor1, label=self.labelGen(P1,X1));
        ax.set_xlabel("Wavelength, $\lambda$ (nm)")
        ax.set_ylabel("Transmission")
        ax.set_title("$S^{atm}(\lambda)$ and $S^{atm,std}(\lambda)$");
        ax.legend(loc='lower right', shadow=False)
        ax.set_xlim(wavelengthRange[0], wavelengthRange[1]);
        
        if (P2 != None) & (X2 != None):
            self.parameterCheck(P2)
            self.airmassCheck(X2)
            atmo2 = self.genAtmo(P2, X2, aerosolNormCoeff=aerosolNormCoeff2, aerosolNormWavelen=aerosolNormWavelen2)
            ax.plot(w,atmo2.sb,label=self.labelGen(P2, X2), alpha=atmo2Alpha, color=atmo2Color)
        elif includeStdAtmo == True:
            atmo2 = self.genAtmo(STDPARAMETERS, STDAIRMASS, aerosolNormCoeff=STDAEROSOLNORMCOEFF, aerosolNormWavelen=STDAEROSOLNORMWAVELEN)
            ax.plot(w, atmo2.sb, label=self.labelGen(STDPARAMETERS, STDAIRMASS), alpha=atmo2Alpha, color=atmo2Color);
        
        ax.legend(loc='lower right', shadow=False)
        
        if figName != None:
            title = figName + "_transPlot.png"
            pylab.savefig(title, format='png')
        return
    
    def filterPlot(self, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN]):
        """Plots the filter response curve from LSST filter data."""
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth, plotHeight)
        
        for f in self.filterlist:
            ax.plot(self.filters[f].wavelen, self.filters[f].sb,label=str(f));
        
        ax.set_xlim(wavelengthRange[0], wavelengthRange[1]);
        ax.set_ylim(0,1);
        ax.set_ylabel("Transmission");
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("LSST $S^{sys}_b$ Filters Only");
        ax.legend(loc=4, shadow=False);
        return
    
    def hardwarePlot(self, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN]):
        """Plots the hardware response curve from LSST hardware data."""
        if self.sys == None:
            self.readHardware()
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth, plotHeight)
        
        for f in self.sys:
            ax.plot(self.sys[f].wavelen, self.sys[f].sb,label=str(f));
        
        ax.set_xlim(wavelengthRange[0], wavelengthRange[1]);
        ax.set_ylim(0,1);
        ax.set_ylabel("Transmission");
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("LSST $S^{sys}_b$ Filters and Hardware");
        ax.legend(loc=4, shadow=False);
        return
    
    def phiPlot(self, bpDict1, bpDict2=None, includeStdAtmo=True, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN],
                phi2Alpha=0.5, phi2Color='black', figName=None):
        """Plots normalized bandpass response function, with the possibility to add a second function
            for comparison."""
        
        w = self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth, plotHeight)
        
        for f in self.filterlist:
            ax.plot(w, bpDict1[f].phi, label=str(f))
            ax.plot(w, bpDict2[f].phi, alpha=phi2Alpha, color=phi2Color)
        
        ax.set_xlim(wavelengthRange[0], wavelengthRange[1]);
        ax.set_ylabel("$\phi_b^{obs}(\lambda)$", fontsize=15);
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("Normalized Bandpass Response");
        ax.legend(loc=4, shadow=False);
        
        if figName != None:
            title = figName + "_phiPlot.png"
            pylab.savefig(title, format='png')
        return
    
    def dphiPlot(self, bpDict1, bpDict2, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN], figName=None):
        """Plots change in normalized bandpass response function given two phi functions."""
        
        w = self.wavelength
        
        fig,ax = pylab.subplots(1,1)
        fig.set_size_inches(plotWidth, plotHeight)
        
        for f in self.filterlist:
            ax.plot(w, bpDict1[f].phi - bpDict2[f].phi, label=str(f))
        
        ax.set_xlim(wavelengthRange[0], wavelengthRange[1]);
        ax.set_ylabel("$\Delta\phi_b^{obs}(\lambda)$", fontsize=15);
        ax.set_xlabel("Wavelength, $\lambda$ (nm)");
        ax.set_title("Change in Normalized Bandpass Response");
        ax.legend(loc=4, shadow=False)
        
        if figName != None:
            title = figName + "_dphiPlot.png"
            pylab.savefig(title, format='png')
        
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


    def dmagsPlot(self, gi, dmags, titletext=None, ylims=None, xlims=None, newfig=True, figName=None, verbose=False):
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
        if verbose == True:
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
            title = figName+"_dmagsPlot.png"
            pylab.savefig(title, format='png')
        
        return
    
    def allPlot(self, P1, X1, P2=None, X2=None, includeStdAtmo=True, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN],
                aerosolNormCoeff1=STDAEROSOLNORMCOEFF, aerosolNormWavelen1=STDAEROSOLNORMWAVELEN,
                aerosolNormCoeff2=STDAEROSOLNORMCOEFF, aerosolNormWavelen2=STDAEROSOLNORMWAVELEN,
                transPlot=True, phiPlot=True, dphiPlot=True, dmagsPlot=True, saveFig=False, figName=None):
        """Plots transmission profile, normalized bandpass response function, change in normalized bandpass response function and change in magnitude 
        against standard atmosphere given an array of parameters and a specific airmass. User can provide array of parameters and airmass in order to 
        override the comparison to the standard atmosphere."""
        
        atmo1 = self.genAtmo(P1, X1, aerosolNormCoeff1, aerosolNormWavelen1)
        bpDict1 = self.combineThroughputs(atmo1)

        if (P2 == None) & (X2 == None):
            if includeStdAtmo == True:
                atmo2 = self.genAtmo(STDPARAMETERS, STDAIRMASS, STDAEROSOLNORMCOEFF, STDAEROSOLNORMWAVELEN)
                bpDict2 = self.combineThroughputs(atmo2)
                figName = self.figNameGen(saveFig, figName, P1, X1, STDPARAMETERS, STDAIRMASS)
        else:
            atmo2 = self.genAtmo(P2, X2, aerosolNormCoeff=aerosolNormCoeff2, aerosolNormWavelen=aerosolNormWavelen2)
            bpDict2 = self.combineThroughputs(atmo2)
            figName = self.figNameGen(saveFig, figName, P1, X1, P2, X2)
        
        if transPlot:
            self.transPlot(P1, X1, P2=P2, X2=X2, includeStdAtmo=includeStdAtmo, plotWidth=plotWidth, plotHeight=plotHeight, wavelengthRange=wavelengthRange,
                           aerosolNormCoeff1=aerosolNormCoeff1, aerosolNormWavelen1=aerosolNormWavelen1,
                           aerosolNormCoeff2=aerosolNormCoeff2, aerosolNormWavelen2=aerosolNormWavelen2, figName=figName)
        if phiPlot:
            self.phiPlot(bpDict1, bpDict2, plotWidth=plotWidth, plotHeight=plotHeight, wavelengthRange=wavelengthRange, figName=figName)
        if dphiPlot:
            self.dphiPlot(bpDict1, bpDict2, plotWidth=plotWidth, plotHeight=plotHeight, wavelengthRange=wavelengthRange, figName=figName)
        if dmagsPlot:
            mag1 = self.mags(bpDict1)
            mag2 = self.mags(bpDict2)
            
            dmags = self.dmags(mag1, mag2)
            gi = self.gi(mag2)
            
            self.dmagsPlot(gi, dmags, figName=figName)
        return

    ### Secondary Functions
    
    def aerosol(self, w, X, alpha, aerosolNormCoeff, aerosolNormWavelen):
        """Standard aerosol transmission function, returns array of transmission values over a range of
            wavelengths."""
        return numpy.e**(-aerosolNormCoeff * X * (aerosolNormWavelen / w) * alpha)
    
    def airmassToString(self, airmass):
        """Converts airmass to string"""
        X = float(airmass)
        return "%.3f" % (X)

    def meanParameter(self, compDict, filters=None):
        """Given a filter-keyed dictionary of best fit values, returns mean value"""
        if filters == None:
            filters = self.filterlist

        meanValue = 0;
        for f in filters:
            meanValue += compDict[f]
            
        return meanValue/float(len(filters))

    def pToString(self, P):
        """Returns string version of parameter array."""
        stringP = "P"
        for i in P:
            if i < 1.0:
                stringP+="0"+str(int(i * 10))
            else:
                stringP+=str(int(i * 10))
        return stringP    
    
    def labelGen(self, P, X):
        """Generates label for use in plot legends."""
        label = []
        for paramNum,param in enumerate(P):
            name = self.parametersPlot[paramNum] + ':'
            labelEle = name + str(param)
            label.append(labelEle)
        return ' '.join(label) + ' $X$:' + str(X)

    def figNameGen(self, saveFig, figName, P1, X1, P2, X2):
        """Generates a string for figure names: 'X1_P1_X2_P2_figName' """
        if saveFig == True:
            if figName != None:
                figName='X'+str(int(X1*10))+'_'+self.pToString(P1)+'_'+'X'+str(int(X2*10))+'_'+self.pToString(P2)+'_'+figName
            else:
                figName='X'+str(int(X1*10))+'_'+self.pToString(P1)+'_'+'X'+str(int(X2*10))+'_'+self.pToString(P2)
        else:
            figName = None
        return figName

    def pickleNameGen(self, comp1, comp2, P, X, Nbins):
        """Generates a string for pickle files. """
        return 'X' + str(int(X*10)) + '_' + self.pToString(P) + '_' + comp1 + '_' + comp2 + '_' + str(Nbins) + 'bins'

    def componentCheck(self,comp,Nbins):
        """Returns a range of values of length Nbins for a given component."""
        if comp == 'H2O':
            return numpy.linspace(0.0,5.0,Nbins), 0
        elif comp == 'O2':
            return numpy.linspace(0.0,5.0,Nbins), 1
        elif comp == 'O3':
            return numpy.linspace(0.0,5.0,Nbins), 2
        elif comp == 'Rayleigh':
            return numpy.linspace(0.0,5.0,Nbins), 3
        elif comp == 'Aerosol':
            return numpy.linspace(0.0,5.0,Nbins), 4
        elif comp == 'Alpha':
            return numpy.linspace(0.0,5.0,Nbins), 5
        else:
            raise ValueError(comp + ' is not a valid component')
        
    def parameterCheck(self, P):
        """Checks if parameter array is of appropriate length."""
        if len(P) != 6:
            raise ValueError('Need 6 parameters to build atmosphere!')
        return

    def airmassCheck(self, X):
        if self.airmassToString(X) not in self.airmasses:
            raise ValueError('Not a valid airmass, check self.airmasses for valid airmasses')
        return
    
    def kuruczCheck(self):
        """Checks if Kurucz model data has been read in."""
        if self.stars == None:
            raise ValueError('No Kurucz model data found, please run self.readKurucz()')
        return
