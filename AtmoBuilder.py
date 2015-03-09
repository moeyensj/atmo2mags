### Necessary imports
import numpy
import os
import copy
import matplotlib.pyplot as plt
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
STDAIRMASS = 1.2
STDAEROSOLNORMCOEFF = 0.1
STDAEROSOLNORMWAVELEN = 550.0
STDAEROSOLALPHA = STDPARAMETERS[5]

SIMSSEDLIBRARY = "SIMS_SED_LIBRARY_DIR"
SEDTYPES = ['kurucz','quasar','galaxy','sn','wd','mlt']

FILTERLIST = ['u','g','r','i','z','y4']

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
        self.filterlist = FILTERLIST
        # List of filter colors
        self.filtercolors = ['blue', 'green', 'red', 'cyan', 'purple', 'yellow']
        
        # Kurucz model data
        self.stars = None
        self.starlist = None
        self.temperature = None
        self.met = None
        self.logg = None

        # White Dwarf Model data
        self.wds = None
        self.wdslist = None
        self.wdslist_H = None
        self.wdslist_He = None
        self.wdtemperature = None
        self.wdlogg = None

        # Galaxy model data
        self.gals = None
        self.gallist = None        
        self.galredshifts = None

        # MLT model data
        self.mlts = None
        self.mltlist = None
        self.mlist = None
        self.llist = None
        self.tlist = None

        # Quasar model data
        self.quasars = None
        self.quasarRedshifts = None

        # Supernova model data
        self.sns = None
        self.snList = None
        self.snDays = None
        self.snRedshifts = None
        
        # Readers
        self.readModtranFiles()
        self.readFilters()
        self.readHardware()

### Reading Functions
    
    def readModtranFiles(self, modtranDir='modtran/', modtranRoot='Pachon_MODTRAN', modtranSuffix='.7sc'):
        """Reads in atmospheric absorption data from MODTRAN files into an airmass-keyed directory."""
        ### Taken from AtmoComp and modified to suit specific needs.
        
        files = os.listdir(modtranDir)
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
            fin = open(os.path.join(modtranDir, file),'r')
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
        homedir = os.getenv(SIMSSEDLIBRARY)  
        stardir = os.path.join(homedir, "starSED/kurucz/")
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

    def readWhiteDwarf(self):
        # read white dwarf bergeron models
        homedir = os.getenv(SIMSSEDLIBRARY)
        whitedwarfdir = os.path.join(homedir, "starSED/wDs/")
        
        allfilelist = os.listdir(whitedwarfdir)
        hlist = []
        helist = []
        temperatures = []
        loggs = []

        for filename in allfilelist:
            tmp = filename.split('_')
            if len(tmp) == 4: # H dwarfs
                temperature = float(tmp[1])
                logg = float(tmp[2].split('.')[0])
                logg = logg/10.0
                if (logg > 7.0) & (temperature>5000):
                    hlist.append(filename)
                    temperatures.append(temperature)
                    loggs.append(logg)
	                
            if len(tmp) == 5: # He dwarfs
                tmp = filename.split('_')
                temperature = float(tmp[2])
                logg = float(tmp[3].split('.')[0])
                logg = logg/10.0
                if (logg > 7.0) & (temperature>5000):
                    helist.append(filename)                        
                    temperatures.append(temperature)
                    loggs.append(logg)

        temperatures = numpy.array(temperatures)
        loggs = numpy.array(loggs)
        wdlist = hlist + helist
        wds = {}
        for w in wdlist:
            wds[w] = Sed()
            if w in hlist:
                wds[w].readSED_flambda(os.path.join(whitedwarfdir, w))
            if w in helist:
                wds[w].readSED_flambda(os.path.join(whitedwarfdir, w))

        print "# Read %d white dwarfs from %s" %(len(wdlist), whitedwarfdir)
        # synchronize seds for faster mag calcs later
        for w in wdlist:
            wds[w].synchronizeSED(wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP)

        self.wds = wds
        self.wdslist = wdlist
        self.wdslist_H = hlist
        self.wdslist_He = helist
        self.wdtemperature = temperatures
        self.wdlogg = loggs

        return 

    def readGalaxies(self, redshiftRange=[0,3.0], redshiftStep=0.5):
        # read sn spectra and redshift
        homedir = os.getenv(SIMSSEDLIBRARY)
        galdir = os.path.join(homedir, "galaxySED/")
        allfilelist = os.listdir(galdir)
        gallist_base = []
        metal = ['002Z', '04Z', '25Z']
        gtype = ['Const', 'Inst', 'Burst'] # removed Exp
        redshifts= numpy.arange(redshiftRange[0], redshiftRange[1]+redshiftStep, redshiftStep)
        # pull out the filenames we want
        for filename in allfilelist:
            if filename.endswith('.spec.gz'):
                tmp = filename.split('.')
                metallicity = tmp[2]
                galaxytype = tmp[0]
                if (metallicity in metal) & (galaxytype in gtype):
                    gallist_base.append(filename)
        # read base SEDs for these  galaxies
        gals_base = {}
        for g in gallist_base:
            gals_base[g] = Sed()
            gals_base[g].readSED_flambda(os.path.join(galdir, g))
        # and redshift
        gals = {}
        gallist = []
        for g in gallist_base:
            for z in redshifts:
                gal_name = "%s_%.1f" %(g, z)
                wavelen, flambda = gals_base[g].redshiftSED(z, wavelen=gals_base[g].wavelen,
                                                            flambda=gals_base[g].flambda)
                gals[gal_name] = Sed(wavelen=wavelen, flambda=flambda)
                gallist.append(gal_name)
        print "# Generated %d galaxies at redshifts between %f and %f" %(len(gallist),
                                                                         redshifts.min(), redshifts.max())
        # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
        for g in gallist:
            gals[g].synchronizeSED(wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP)
        # add dust
        ax, bx = gals[gallist[0]].setupCCMab()
        for g in gallist:
            gals[g].addCCMDust(ax, bx, A_v=0.02)

        self.gals = gals
        self.gallist = gallist
        self.galredshifts = redshifts 

        return 

    def readMLT(self):
        # read mlt stars - only keep 'm's
        # find the filenames and mark 'm', 'l', 't' stars separately
        homedir = os.getenv(SIMSSEDLIBRARY)
        mltdir = os.path.join(homedir, "starSED/mlt/")
        allfilelist = os.listdir(mltdir)
        mltlist = []
        mlist = []
        llist = []
        tlist = []
        for filename in allfilelist:
            if filename.startswith('m'):
                mlist.append(filename)
            elif filename.startswith('L'):
                llist.append(filename)
            elif filename.startswith('burrows'):
                tlist.append(filename)
        mltlist = mlist # + llist + tlist
        # read the mlt seds from disk
        mlts = {}
        for s in mltlist:
            mlts[s] = Sed()
            mlts[s].readSED_flambda(os.path.join(mltdir, s))
        print "# Read %d mlt stars from %s" %(len(mltlist), mltdir)
        # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
        for s in mltlist:
            mlts[s].synchronizeSED(wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP)

        self.mlts = mlts
        self.mltlist = mltlist
        self.mlist = mlist
        self.llist = llist
        self.tlist = tlist
        
        return 

    def readQuasar(self, redshiftRange=[0,7.5], redshiftStep=0.1):
        homedir = os.getenv("HOME")
        quasardir = os.path.join(homedir, "atmo2mags/seds/quasar")
        # read zero redshift quasar
        base = Sed()
        base.readSED_flambda(os.path.join(quasardir, "quasar.dat"))
        # redshift 
        redshifts= numpy.arange(redshiftRange[0], redshiftRange[1]+redshiftStep, redshiftStep)
        quasars = {}
        for z in redshifts:
            wavelen, flambda = base.redshiftSED(z, wavelen=base.wavelen, flambda=base.flambda)
            quasars[z] = Sed(wavelen=wavelen, flambda=flambda)
        print "# Generated %d quasars at redshifts between %f and %f" %(len(redshifts), redshifts.min(), redshifts.max())
        # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
        for z in redshifts:
            quasars[z].synchronizeSED(wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP)

        self.quasars = quasars
        self.quasarRedshifts = redshifts

        return

    def readSNes(self, redshiftRange=[0,1.2], redshiftStep=0.1, days=['0', '20', '40']):
        # read sn spectra and redshift
        homedir = os.getenv("HOME")
        sndir = os.path.join(homedir, "atmo2mags/seds/sn")
        allfilelist = os.listdir(sndir)
        snlist = []
        redshifts= numpy.arange(redshiftRange[0], redshiftRange[1]+redshiftStep, redshiftStep)
        # pull out the filenames we want
        for filename in allfilelist:
            if filename.endswith('.dat') & filename.startswith('sn1a_'):
                tmp = filename.split('_')
                day = tmp[1].split('.')[0]
                if day in days:
                    snlist.append(filename)
        # read base SEDs for these days
        sns_base = {}
        for r in zip(snlist, days):
            day = r[1]
            sns_base[day] = Sed()
            sns_base[day].readSED_flambda(os.path.join(sndir, r[0]))

        # and redshift
        sns = {}
        snlist = []
        for d in days:        
            for z in redshifts:
                sn_name = "%d_%.1f" %(int(d), z)
                wavelen, flambda = sns_base[d].redshiftSED(z, wavelen=sns_base[d].wavelen, flambda=sns_base[d].flambda)
                sns[sn_name] = Sed(wavelen=wavelen, flambda=flambda)
                snlist.append(sn_name)

        print "# Generated %d sn's at redshifts between %f and %f on days %s" %(len(snlist),redshifts.min(), redshifts.max(), days)
        # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
        for s in snlist:
            sns[s].synchronizeSED(wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP)

        self.sns = sns
        self.snList = snlist
        self.snDays = days
        self.snRedshifts = redshifts
        
        return

    def readAll(self, kurucz=True, galaxies=True, whiteDwarfs=True, mltDwarfs=True, quasars=True, SNes=True):
        if kurucz:
            self.readKurucz()
        if whiteDwarfs:
            self.readWhiteDwarf()
        if mltDwarfs:
            self.readMLT()
        if galaxies:
            self.readGalaxies()
        if quasars:
            self.readQuasar()
        if SNes:
            self.readSNes()
        return

### Calculator / Generator Functions
        
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
    
    def mags(self, bpDict, seds=None, sedkeylist=None, filters=None, verbose=False):
        """Calculates magnitudes given a bandpass dictionary, returns filter-keyed magnitude dictionary. If seds and sedkeylist are not none
        returns mags for Kurucz model MS stars."""
        ### Taken from plot_dmags and modified to suit specific needs.
        # calculate magnitudes for all sed objects using bpDict (a single bandpass dictionary keyed on filters)
        # pass the sedkeylist so you know what order the magnitudes are arranged in
        if filters == None:
            filters = self.filterlist

        if filters == 'y4':
            filters = ['y4']

        mags = {}
        for f in filters:
            mags[f] = numpy.zeros(len(sedkeylist), dtype='float')
            for i,key in enumerate(sedkeylist):
                mags[f][i] = seds[key].calcMag(bpDict[f])
                if numpy.isnan(mags[f][i]):
                    print key, f, mags[f][i]
            if verbose == True:
                print f, mags[f].max(), mags[f].min()
        return mags
    
    def dmags(self, mags1, mags2, filters=None):
        """Returns filter-keyed dictionary of change in magnitude in millimagnitudes."""
        ### Taken from plot_dmags and modified to suit specific needs.
        if filters == None:
            filters = self.filterlist

        if filters == 'y4':
            filters = ['y4']
        
        dmags = {}
        for f in filters:
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

    def compute_logL(self, P, X, err, f, mags_obs, mags_std, seds, sedkeylist, deltaGrey):
        """Return logL for a given array of parameters P, airmass X, error, a filter and the magnitudes of a standard atmosphere."""
        atmo = self.genAtmo(P,X)
        throughputAtmo = self.combineThroughputs(atmo)
        mags_fit = self.mags(throughputAtmo, seds=seds, sedkeylist=sedkeylist, filters=f)
        
        dmags_fit = self.dmags(mags_fit, mags_std, filters=f)
        dmags_obs = self.dmags(mags_obs, mags_std, filters=f)

        if deltaGrey >= 0.0:
            dmags_fit[f] = dmags_fit[f] - deltaGrey
            dmags_obs[f] = dmags_obs[f] - deltaGrey
        else: # (deltaGrey < 0.0)
            dmags_fit[f] = dmags_fit[f] - numpy.mean(dmags_fit[f])
            dmags_obs[f] = dmags_obs[f] - numpy.mean(dmags_obs[f])
    
        return -numpy.sum(0.5 * ((dmags_fit[f] - dmags_obs[f]) / err) ** 2)

    def compute_mag_color_nonlinear(self, comp1, comp2, P_obs, X_obs, err=0.005, Nbins=50, deltaGrey=0.0, regressionSed='kurucz', 
        comparisonSeds=SEDTYPES, generateFig=True, generateDphi=True, pickleString=None, filters=None, verbose=True):
        # Insure valid parameters, airmass and sedtypes are given
        self.parameterCheck(P_obs)
        self.airmassCheck(X_obs)
        self.sedTypeCheck(regressionSed)

        # Find range over which to vary parameter and the parameter number for comp1, comp2
        range1, pNum1 = self.componentCheck(comp1,Nbins)
        range2, pNum2 = self.componentCheck(comp2,Nbins)

        # Find seds and sedkeylist for sedtype
        seds, sedkeylist = self.sedFinder(regressionSed)
        
        if filters == None:
            filters = self.filterlist

        if verbose:
            print 'Computing nonlinear regression for ' + comp1 + ' and ' + comp2 + '.'
            print 'Observed atmosphere parameters: ' + str(P_obs)
            print 'Observed atmosphere airmass:    ' + str(X_obs)
            print 'Standard atmosphere parameters: ' + str(STDPARAMETERS)
            print 'Standard atmosphere airmass:    ' + str(STDAIRMASS)
            print 'Observed atmosphere parameter for ' + comp1 + ': ' + str(P_obs[pNum1])
            print 'Observed atmosphere parameter for ' + comp2 + ': ' + str(P_obs[pNum2])
            print ''
        
        P_fit = copy.deepcopy(P_obs)
        X_fit = copy.deepcopy(X_obs)

        # Create observed atmosphere
        obs = self.genAtmo(P_obs,X_obs)
        throughput_obs = self.combineThroughputs(obs)
        mags_obs = self.mags(throughput_obs, seds=seds, sedkeylist=sedkeylist, filters=filters)

        # Create standard atmosphere
        std = self.genAtmo(STDPARAMETERS,STDAIRMASS)
        throughput_std = self.combineThroughputs(std)
        mags_std = self.mags(throughput_std, seds=seds, sedkeylist=sedkeylist, filters=filters)

        logL = {}
        whr = {}
        comp1best = {}
        comp2best = {}

        for f in filters:
            if pickleString != None:
                pickString_temp = self.pickleNameGen(comp1, comp2, P_obs, X_obs, Nbins, regressionSed, deltaGrey, f=f) + '_' + pickleString + '.pkl'
            else:
                pickleString_temp = self.pickleNameGen(comp1, comp2, P_obs, X_obs, Nbins, regressionSed, deltaGrey, f=f) + '.pkl'
                    
            @pickle_results(pickleString_temp)
            def run_regression(comp1, comp2, f):
                
                logL = []
                whr = []
                comp1best = []
                comp2best = []

                print 'Calculating best parameters for ' + f + ' filter...'
                logL = numpy.empty((Nbins, Nbins))
                for i in range(len(range1)):
                    for j in range(len(range2)):
                        P_fit[pNum1] = range1[i]
                        P_fit[pNum2] = range2[j]
                        logL[i, j] = self.compute_logL(P_fit, X_fit, err, f, mags_obs, mags_std, seds, sedkeylist, deltaGrey)

                print 'Completed ' + f + ' filter.'

                for f in filters:
                    logL -= numpy.max(logL)
                    whr = numpy.where(logL == numpy.max(logL))
                    comp1best = range1[whr[0][0]]
                    comp2best = range2[whr[1][0]]

                return comp1best, comp2best, logL

            comp1best[f], comp2best[f], logL[f]  = run_regression(comp1, comp2, f)
            pickleString_temp = ''

        if pickleString != None:
            pickleString = self.pickleNameGen(comp1, comp2, P_obs, X_obs, Nbins, regressionSed, deltaGrey) + '_' + pickleString
        else:
            pickleString = self.pickleNameGen(comp1, comp2, P_obs, X_obs, Nbins, regressionSed, deltaGrey)

        if generateDphi:
            self.dphiPlot(throughput_obs, throughput_std, figName=pickleString)

        if verbose:
        	print ''
        	print 'Best fit parameters (Filter, ' + comp1 + ', ' + comp2 + '):'
        	for f in filters:
 				print '%s %.2f %.2f' % (f, comp1best[f], comp2best[f])

        if generateFig == True:
            self.regressionPlot(comp1, comp1best, comp2, comp2best, logL, P_obs, X_obs, pNum1=pNum1, pNum2=pNum2,
                                comp1_range=range1, comp2_range=range2, Nbins=Nbins, figName=pickleString, 
                                regressionSed=regressionSed, comparisonSeds=comparisonSeds, filters=filters, verbose=verbose)

        return range1, range2, comp1best, comp2best, logL

### Plotting Functions

    def regressionPlot(self, comp1, comp1_best, comp2, comp2_best, logL, P_obs, X_obs, pNum1=None, pNum2=None,
        comp1_range=None, comp2_range=None, Nbins=50, regressionSed='kurucz', comparisonSeds=SEDTYPES, plotDifference=True, 
        figName=None, filters=None , verbose=True):
        """Plots dmags with each filter in its own subplot."""
        ### Taken from plot_dmags and modified to suit specific needs.

        if any([pNum1,pNum2,comp1_range,comp2_range]) == None:
            comp1_range, pNum1 = self.componentCheck(comp1, Nbins)
            comp2_range, pNum2 = self.componentCheck(comp2, Nbins)
            
        if filters == None:
            filters = self.filterlist

        seds, sedkeylist = self.sedFinder(regressionSed)
        
        fig, ax = plt.subplots(len(filters),3)
        fig.suptitle(r'$\Delta$mmags, Regression Contours $\Delta\Delta$mmags for each LSST filter', fontsize=14)
        fig.set_size_inches(15,len(filters)*5)
        fig.subplots_adjust(top=0.93, wspace=0.20, hspace=0.20, bottom=0.09, left=0.10, right=0.96)

        # Save observed parameters
        comp1_obs = P_obs[pNum1]
        comp2_obs = P_obs[pNum2]
    
        # Create observed atmosphere
        obs = self.genAtmo(P_obs,X_obs)
        throughput_obs = self.combineThroughputs(obs)
        mags_obs = self.mags(throughput_obs, seds=seds, sedkeylist=sedkeylist, filters=filters)

        # Create standard atmosphere
        std = self.genAtmo(STDPARAMETERS,STDAIRMASS)
        throughput_std = self.combineThroughputs(std)
        mags_std = self.mags(throughput_std, seds=seds, sedkeylist=sedkeylist, filters=filters)

        P_fit = copy.deepcopy(P_obs)
        X_fit = copy.deepcopy(X_obs)

        # For each filter plot dmags and regression contours
        for i,f in enumerate(filters):
            # Set component parameters to best fit parameters
            P_fit[pNum1] = comp1_best[f]
            P_fit[pNum2] = comp2_best[f]

            # Create atmosphere at best fit parameters
            fit = self.genAtmo(P_fit,X_fit)
            throughput_fit = self.combineThroughputs(fit)

            self.dmagSED(ax[i][0], f, throughput_fit, throughput_std, regressionSed)
            self.dmagSED(ax[i][0], f, throughput_obs, throughput_std, regressionSed, truth=True)

            # Plot parameter space regression plots
            # Plot contours and true values
            ax[i][1].contour(comp1_range, comp2_range, convert_to_stdev(logL[f].T), levels=(0.683, 0.955, 0.997), colors='k')
            ax[i][1].contour(comp1_range, comp2_range, logL[f].T, colors='k')
            ax[i][1].scatter(comp1_obs, comp2_obs, marker='o', s=25, facecolors='none', edgecolors='b', label='Truth')

            # Plot dashed lines at best fit parameters
            ax[i][1].axvline(comp1_best[f], color='black', linestyle='--', label='Fit')
            ax[i][1].axhline(comp2_best[f], color='black', linestyle='--')

            # Set y-axis, x-axis limits
            ax[i][1].set_xlim(min(comp1_range), max(comp1_range))
            ax[i][1].set_ylim(min(comp2_range), max(comp2_range))

            # Label axes
            str1 = r'%s (fit: %.2f, truth: %.2f)' % (comp1, comp1_best[f], comp1_obs)
            str2 = r'%s (fit: %.2f, truth: %.2f)' % (comp2, comp2_best[f], comp2_obs)
            ax[i][1].set_xlabel(str1)
            ax[i][1].set_ylabel(str2)

            # Plot dmags for other SEDS:
            
            if plotDifference == False:
                for s in comparisonSeds:
                    if s != regressionSed:
                        self.dmagSED(ax[i][2], f, throughput_fit, throughput_std, s, dmaglimit=False)

                for s in comparisonSeds:
                    if s != regressionSed:
                        self.dmagSED(ax[i][2], f, throughput_obs, throughput_std, s, dmaglimit=False, truth=True)
            else:
                for s in comparisonSeds:
                    if s != regressionSed:
                        self.dmagSED(ax[i][2], f, throughput_fit, throughput_std, s, bpDict2=throughput_obs)

            if i == 0:
                ax[i][0].legend(loc='upper center', bbox_to_anchor=(0.5,1.25), ncol=2)
                ax[i][1].legend(loc='upper center', bbox_to_anchor=(0.5,1.25), ncol=3)
                ax[i][2].legend(loc='upper center', bbox_to_anchor=(0.5,1.25), ncol=3)
            
        if figName != None:
            title = figName+"_regressionPlot.png"
            plt.savefig(title, format='png')

        return

    def dmagSED(self, ax, f, bpDict1, bpDict_std, sedtype, bpDict2=None, truth=False, dmaglimit=True):
        # Check if valid sedtype, check if sed data read:
        self.sedTypeCheck(sedtype)
        self.sedReadCheck(sedtype)

        # Label axes, add grid
        ax.set_xlabel("g-i")
        ax.set_ylabel(r"$\Delta$ %s (mmag)" %(f))
        ax.grid(b=True)

        if sedtype == 'kurucz':
            mags = self.mags(bpDict1, seds=self.stars, sedkeylist=self.starlist)
            mags_std = self.mags(bpDict_std, seds=self.stars, sedkeylist=self.starlist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=self.stars, sedkeylist=self.starlist)
                dmags2 = self.dmags(mags2, mags_std)

            metallicity = numpy.array(self.met)
            logg = numpy.array(self.logg)
            metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
            metbinsize = abs(metallicity.min() - metallicity.max())/6.0
            metbins = numpy.arange(metallicity.min(), metallicity.max() + metbinsize, metbinsize)

            for metidx in range(len(metbins)):
                # Make cut of stars
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize) \
                        & (logg>3.5))
                mcolor = metcolors[metidx]

                if bpDict2 != None:
                    ax.plot(gi[condition], dmags[f][condition]-dmags2[f][condition], mcolor+'.')
                else:
                    if truth == True:
                        if metidx == len(metbins)-1:
                            ax.plot(gi[condition], dmags[f][condition], mcolor+'.', label='Truth')
                        else:
                            ax.plot(gi[condition], dmags[f][condition], mcolor+'.')
                    else:
                        if metidx == len(metbins)-1:
                            ax.plot(gi[condition], dmags[f][condition], mcolor+'.', color='gray', label='Fit')
                        else:
                            ax.plot(gi[condition], dmags[f][condition], mcolor+'.', color='gray')

        elif sedtype == 'quasar':
            mags = self.mags(bpDict1, seds=self.quasars, sedkeylist=self.quasarRedshifts)
            mags_std = self.mags(bpDict_std, seds=self.quasars, sedkeylist=self.quasarRedshifts)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=self.quasars, sedkeylist=self.quasarRedshifts)
                dmags2 = self.dmags(mags2, mags_std)

            redshift = self.quasarRedshifts
            redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
            redbinsize = 0.5
            redbins = numpy.arange(0.0, 3.0+redbinsize, redbinsize)
            for redidx in range(len(redbins)):
                condition =((redshift>=redbins[redidx]) & (redshift<=redbins[redidx]+redbinsize))
                rcolor = redcolors[redidx]

                if bpDict2 != None:
                    ax.plot(gi[condition], dmags[f][condition]-dmags2[f][condition], rcolor+'.')
                else:
                    if truth == True:
                        ax.plot(gi[condition], dmags[f][condition], rcolor+'.')
                    else:
                        ax.plot(gi[condition], dmags[f][condition], rcolor+'.', color='gray')
        
        elif sedtype == 'galaxy':
            mags = self.mags(bpDict1, seds=self.gals, sedkeylist=self.gallist)
            mags_std = self.mags(bpDict_std, seds=self.gals, sedkeylist=self.gallist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=self.gals, sedkeylist=self.gallist)
                dmags2 = self.dmags(mags2, mags_std)

            gallist = self.gallist
            redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
            redbinsize = 0.5
            redbins = numpy.arange(0.0, 3.0+redbinsize, redbinsize)

            for i,g in enumerate(gallist):
                galbase, redshift = g.split('_')
                redshift = float(redshift)
                redidx = int(redshift / redbinsize)

                if bpDict2 != None:
                    ax.plot(gi[i], dmags[f][i]-dmags2[f][i], redcolors[redidx]+'.')
                else:
                    if truth == True:
                        ax.plot(gi[i], dmags[f][i], redcolors[redidx]+'.')
                    else:
                        ax.plot(gi[i], dmags[f][i], redcolors[redidx]+'.', color='gray')

        elif sedtype == 'mlt':
            mags = self.mags(bpDict1, seds=self.mlts, sedkeylist=self.mltlist)
            mags_std = self.mags(bpDict_std, seds=self.mlts, sedkeylist=self.mltlist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=self.mlts, sedkeylist=self.mltlist)
                dmags2 = self.dmags(mags2, mags_std)

            mltlist = self.mltlist
            mlist = self.mlist
            llist = self.llist
            tlist = self.tlist
        
            for j in range(len(mltlist)):
                if bpDict2 != None:
                    if (mltlist[j] in mlist):
                        ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'bx')
                    elif (mltlist[j] in llist):
                        ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'gx')
                    elif (mltlist[j] in tlist):
                        ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'mx')
                else:
                    if truth == True:
                        if (mltlist[j] in mlist):
                            ax.plot(gi[j], dmags[f][j], 'bx')
                        elif (mltlist[j] in llist):
                            ax.plot(gi[j], dmags[f][j], 'gx')
                        elif (mltlist[j] in tlist):
                            ax.plot(gi[j], dmags[f][j], 'mx')
                    else:
                        if (mltlist[j] in mlist):
                            ax.plot(gi[j], dmags[f][j], marker='x', color='gray')
                        elif (mltlist[j] in llist):
                            ax.plot(gi[j], dmags[f][j], marker='x', color='gray')
                        elif (mltlist[j] in tlist):
                            ax.plot(gi[j], dmags[f][j], marker='x', color='gray')

        elif sedtype == 'wd':
            mags = self.mags(bpDict1, seds=self.wds, sedkeylist=self.wdslist)
            mags_std = self.mags(bpDict_std, seds=self.wds, sedkeylist=self.wdslist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=self.wds, sedkeylist=self.wdslist)
                dmags2 = self.dmags(mags2, mags_std)

            wdslist = self.wdslist
            hlist = self.wdslist_H
            helist = self.wdslist_He

            for j in range(len(wdslist)):
                if bpDict2 != None:
                    if (wdslist[j] in hlist):
                        ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'y+')
                    elif (wdslist[j] in helist):
                        ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'y+')
                else:
                    if truth == True:
                        if (wdslist[j] in hlist):
                            ax.plot(gi[j], dmags[f][j], 'y+')
                        elif (wdslist[j] in helist):
                            ax.plot(gi[j], dmags[f][j], 'y+')
                    else:
                        if (wdslist[j] in hlist):
                            ax.plot(gi[j], dmags[f][j], marker='+', color='gray')
                        elif (wdslist[j] in helist):
                            ax.plot(gi[j], dmags[f][j], marker='+', color='gray')

        elif sedtype == 'sn':
            mags = self.mags(bpDict1, seds=self.sns, sedkeylist=self.snList)
            mags_std = self.mags(bpDict_std, seds=self.sns, sedkeylist=self.snList)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=self.sns, sedkeylist=self.snList)
                dmags2 = self.dmags(mags2, mags_std)
 
            snlist = self.snList
            redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
            redbinsize = 0.2
            redbins = numpy.arange(0.0, 1.2+redbinsize, redbinsize)
            day_symbol = {'0':'s', '20':'s', '40':'s'}

            for j,s in enumerate(snlist):
                day, redshift = s.split('_')
                redshift = float(redshift)
                redidx = int(redshift / redbinsize)

                if bpDict2 != None:
                    ax.plot(gi[j], dmags[f][j]-dmags2[f][j], redcolors[redidx]+day_symbol[day])
                else:
                    if truth == True:
                        ax.plot(gi[j], dmags[f][j], redcolors[redidx]+day_symbol[day])
                    else:
                        ax.plot(gi[j], dmags[f][j], redcolors[redidx]+day_symbol[day], color='gray')

        # Add appropriate y-axis limits
        if dmaglimit == True:
            self.dmagLimit(ax, f, dmags)
            if bpDict2 != None:
                self.dmagLimit(ax, f, dmags2)

        return

    def transPlot(self, P1, X1, P2=None, X2=None, includeStdAtmo=True, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN],
        aerosolNormCoeff1=STDAEROSOLNORMCOEFF, aerosolNormWavelen1=STDAEROSOLNORMWAVELEN,
        aerosolNormCoeff2=STDAEROSOLNORMCOEFF, aerosolNormWavelen2=STDAEROSOLNORMWAVELEN,
        atmoColor1='blue', atmo2Color='black', atmo2Alpha=0.5, figName=None):
        """Plots atmospheric transmission profile given a parameter array."""
        
        w=self.wavelength

        atmo1 = self.genAtmo(P1, X1, aerosolNormCoeff=aerosolNormCoeff1, aerosolNormWavelen=aerosolNormWavelen1)
        
        fig,ax = plt.subplots(1,1)
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
            plt.savefig(title, format='png')
        return
    
    def filterPlot(self, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN]):
        """Plots the filter response curve from LSST filter data."""
        
        fig,ax = plt.subplots(1,1)
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
        
        fig,ax = plt.subplots(1,1)
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
        
        fig,ax = plt.subplots(1,1)
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
            plt.savefig(title, format='png')
        return
    
    def dphiPlot(self, bpDict1, bpDict2, plotWidth=12, plotHeight=6, wavelengthRange=[MINWAVELEN,MAXWAVELEN], figName=None):
        """Plots change in normalized bandpass response function given two phi functions."""
        
        w = self.wavelength
        
        fig,ax = plt.subplots(1,1)
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
            plt.savefig(title, format='png')
        
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
            ax = plt.subplot(3,2,i)
            for metidx in range(len(metbins)):
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
                mcolor = metcolors[metidx]
                plt.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], mcolor+'.')
            i = i + 1
        ax = plt.subplot(3,2,7)
        for metidx in range(len(metbins)):
            condition = ((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
            mcolor = metcolors[metidx]
            plt.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], mcolor+'.')
        # set up generic items
        for i in range(1, 7):
            f = filterlist[i-1]
            ax = plt.subplot(3,2,i)
            #plt.xlabel("g-i")
            #plt.ylabel(r"$\Delta$ %s (mmag)" %(f))
            def axis_formatter(x, pos):
                return "%.1f" %(x)
            formatter = plt.FuncFormatter(axis_formatter)
            ax.yaxis.set_major_formatter(formatter)
            # set axes limits
            if ylims == None:
                pass
            else:
                try:
                    plt.ylim(ylims[f][0], ylims[f][1])
                except KeyError:
                    pass
        # put a grid in the background
        if newfig:
            for i in range(1, 7):
                ax = plt.subplot(3, 2, i)
                plt.grid(True)
                #plt.suptitle(titletext)
                #plt.savefig("delta_mags2.eps", format='eps')
        return

    def dmagsPlot(self, gi, dmags, titletext=None, ylims=None, xlims=None, newfig=True, figName=None, verbose=False):
        """Plots dmags with each filter in its own subplot."""
        ### Taken from plot_dmags and modified to suit specific needs.
        sedcolorkey = [self.met,self.logg]
        plotfilterlist = self.filterlist
        
        # make figure of change in magnitude
        if newfig:
            plt.figure(figsize=(10,15))
        yplots = 3
        xplots = 2
        if len(plotfilterlist) == 1:
            yplots = 1
            xplots = 1
        plt.subplots_adjust(top=0.93, wspace=0.32, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
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
            ax = plt.subplot(yplots, xplots,i)
            for metidx in range(len(metbins)):
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize) \
                            & (logg>3.5))
                mcolor = metcolors[metidx]
                plt.plot(gi[condition], dmags[f][condition], mcolor+'.')
                if plot_logg2:
                        condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize) \
                            & (logg<2.5))
                        mcolor = metcolors[metidx]
                        plt.plot(gi[condition], dmags[f][condition], mcolor+'x')
            i = i + 1
        # set up generic items
        for i in range(0, len(plotfilterlist)):
            f = plotfilterlist[i]
            ax = plt.subplot(yplots,xplots,i+1)
            plt.xlabel("g-i")
            plt.ylabel(r"$\Delta$ %s (mmag)" %(f))
            def axis_formatter(x, pos):
                return "%.1f" %(x)
            formatter = plt.FuncFormatter(axis_formatter)
            ax.yaxis.set_major_formatter(formatter)
            # set axes limits
            if ylims == None:
                pass
            else:
                try:
                    plt.ylim(ylims[f][0], ylims[f][1])
                except KeyError:
                    pass
            if xlims == None:
                pass
            else:
                try:
                    plt.xlim(xlims[f][0], xlims[f][1])
                except KeyError:
                    pass
        # put a grid in the background
        if newfig:
            for i in range(0, len(plotfilterlist)):
                ax = plt.subplot(yplots,xplots, i+1)
                plt.grid(True)
            if titletext!=None:
                plt.suptitle("$\Delta$mmags for each LSST filter")
                
        if figName != None:
            title = figName+"_dmagsPlot.png"
            plt.savefig(title, format='png')
        
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

    def pToString(self, P):
        """Returns string version of parameter array."""
        stringP = "P"
        for i in P:
            if i < 1.0:
                stringP+="0"+str(int(i * 10))
            else:
                stringP+=str(int(i * 10))
        return stringP

    def meanParameter(self, compDict, filters=None):
        """Given a filter-keyed dictionary of best fit values, returns mean value"""
        if filters == None:
            filters = self.filterlist

        meanValue = 0;
        for f in filters:
            meanValue += compDict[f]
            
        return meanValue/float(len(filters))

    def dmagLimit(self, ax, f, dmags):
        # Initialize values to keep track of dmag range
        dmag_max = numpy.max(dmags[f]);
        dmag_min = numpy.min(dmags[f]);
        dmag_range = abs((dmag_max - dmag_min)/2.0);

        # If dmag range exceeds 2.0, plot dashed lines at +-2dmmags
        if dmag_range > 2.0:
            ax.axhline(2,color='black',linestyle='--')
            ax.axhline(-2,color='black',linestyle='--')
        else:
            ax.set_ylim(-2,2)

        return
    
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

    def pickleNameGen(self, comp1, comp2, P, X, Nbins, regressionSed, deltaGrey, f=None):
        """Generates a string for pickle files. """
        s1 = 'X' + str(int(X*10))
        s2 = self.pToString(P)
        s3 = comp1 + '_' + comp2
        s4 = 'XSTD' + str(int(STDAIRMASS*10))
        s5 = 'DG' + str(int(deltaGrey*10.0))
        s6 = ''  
        s7 = str(Nbins) + 'bins'

        if f != None:
            s6 = regressionSed + '_' + f
        else:
            s6 = regressionSed 

        return '%s_%s_%s_%s_%s_%s_%s' % (s1, s2, s3, s4, s5, s6, s7)

    def componentCheck(self, comp, Nbins):
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

    def sedTypeCheck(self, sedtype):
        sedtypes = SEDTYPES
        if sedtype not in sedtypes:
            raise ValueError(str(sedtype) + ' is not a valid SED type, valid SED types: ' + str(sedtypes))
    
    def sedReadCheck(self, sedtype):
        """Checks if Kurucz model data has been read in."""
        if sedtype == 'kurucz':
            if self.stars == None:
                raise ValueError('No Kurucz model data found, please run self.readKurucz() or self.readAll()')
        elif sedtype == 'quasar':
            if self.quasars == None:
                raise ValueError('No quasar data found, please run self.readQuasars() or self.readAll()')
        elif sedtype == 'galaxy':
            if self.gals == None:
                raise ValueError('No galaxy data found, please run self.readGalaxies() or self.readAll()')
        elif sedtype == 'wds':
            if self.wds == None:
                raise ValueError('No white dwarf data found, please run self.readWhiteDwarf or self.readAll()')
        elif sedtype == 'mlt':
            if self.mlts == None:
                raise ValueError('No mlt dwarf data found, please run self.readMLT or self.readAll()')
        elif sedtype == 'sn':
            if self.sns == None:
                raise ValueError('No supernova data found, please run self.readSNes or self.readAll()')
        return

    def sedFinder(self, sedtype):
        """Returns seds and sedkeylist given an sedtype."""
        seds = []
        sedkeylist = []

        if sedtype == 'kurucz':
            seds = self.stars
            sedkeylist = self.starlist
        elif sedtype == 'quasar':
            seds = self.quasars
            sedkeylist = self.quasarRedshifts
        elif sedtype == 'galaxy':
            seds = self.gals
            sedkeylist = self.gallist
        elif sedtype == 'wds':
            seds = self.wds
            sedkeylist = self.wdslist
        elif sedtype == 'mlt':
            seds = self.mlts
            sedkeylist = self.mltlist
        elif sedtype == 'sn':
            seds = self.sns
            sedkeylist = self.snList

        return seds, sedkeylist




