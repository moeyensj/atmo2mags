### Necessary imports
import numpy as np
import os
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass

from Atmo import Atmo

from astroML.plotting.mcmc import convert_to_stdev
from astroML.decorators import pickle_results

# Global wavelength variables set to MODTRAN defaults
WAVELENMIN = 300
WAVELENMAX = 1100
WAVELENSTEP = 0.5

STDPARAMETERS = [1.0,1.0,1.0,1.0,1.0,1.7]
STDAIRMASS = 1.2
STDAEROSOLNORMCOEFF = 0.1
STDAEROSOLNORMWAVELEN = 550.0
STDAEROSOLALPHA = STDPARAMETERS[5]

SIMSSEDLIBRARY = 'SIMS_SED_LIBRARY_DIR'
PICKLEDIRECTORY = 'pickles/'
PLOTDIRECTORY = 'plots/'
LOGLDIRECTORY = 'logls/'

SEDTYPES = ['mss','qsos','gals','sns','wds','mlts']
FILTERLIST = ['u','g','r','i','z','y4']
COLORS = ['g-i','u-g','g-r','r-i','i-z','z-y','z-y4']

FIGUREWIDTH = 10
FIGUREHEIGHT = 7
TITLESIZE = 14
LABELSIZE = 13

class AtmoBuilder(object):
    """
    Functions:
    ----------------------
    function: 
        description

    ## Reading Functions ###

    readModtranFiles: 
        Reads atmospheric absorption data into an airmass-keyed directory from MODTRAN files.

    readFilters: 
        Reads LSST filter data and returns a filter-keyed dictionary. (S^{filters})
    
    readHardware: 
        Reads LSST hardware data and returns a filter-keyed dictionary. (S^{sys})

    readMSs: 
        Reads Kurucz main sequence star model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.mss            # list of main sequence stars (SED objects)
        self.msList         # list of main sequence stars
        self.msMet          # list of main sequence star metallicities
        self.msTemp         # list of main sequence star temperatures
        self.msLogg         # list of main sequence star log surface gravities

    readWDs: 
        Reads white dwarf model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.wds            # list of white dwarfs (SED objects)
        self.wdList         # list of white dwarfs
        self.wdListH        # list of hydrogen white dwarfs 
        self.wdListHe       # list of helium white dwarfs
        self.wdTemp         # list of white dwarf temperatures
        self.wdLogg         # list of white dwarf surface gravities

    readGals: 
        Reads galaxy model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.gals           # list of galaxies (SED objects)
        self.galList        # list of galaxies      
        self.galRedshifts   # list of galaxy redshifts

    readMLTs: 
        Reads MLT dwarf model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.mlts           # list of mlt dwarfs (SED objects)
        self.mltList        # list of mlt dwarfs
        self.mList          # list of m dwarfs
        self.lList          # list of l dwarfs
        self.tList          # list of t dwarfs

    readQsos: 
        Reads quasar model data and sets relevant class attributes.

        The following attributes will be set:
        self.qsos           # list of quasars (SED objects)    
        self.qsoRedshifts   # list of quasar redshifts

    readSns: 
        Reads supernova model data and sets relevant class attributes.

        The following attributes will be set:
        self.sns            # list of supernova (SED objects)
        self.snList         # list of supernova
        self.snDays         # list of supernova days
        self.snRedshifts    # list of supernova redshifts

    readAll: 
        Read all (or subset of) SED model data.


    ### Calculator / Generator Functions ###

    buildAtmo: 
        Builds an atmospheric transmission profile (as an atmo object) given a set of component parameters 
        and an airmass. (S^{atm})

    combineThroughputs: 
        Combines atmospheric transmission profile with system responsiveness data, returns filter-keyed 
        dictionary. (S^{atm}*S^{sys})

    mags: 
        Calculates magnitudes given a bandpass dictionary, returns filter-keyed magnitude dictionary. If seds and 
        sedkeylist are not none returns mags for Kurucz model MS stars.

    dmags: 
        Returns filter-keyed dictionary of change in magnitude in millimagnitudes.

    gi: 
        Returns standard color temperature given standard magnitude dictionary keyed on filters.


    ### Regression Functions ###

    computeAtmoFit: 
        Computes the best fit atmospheric parameters for two given components and an observed atmosphere. 
        Requires the SED data for the specified regression and comparison SEDs to be read in. 


    ### Plotting Functions ###

    regressionPlot: 
        Plots regression data with each filter in its own row of subplots. Requires the 
        SED data for the specified regression and comparison SEDs to be read in.

    transPlot: 
        Plots atmospheric transmission profile given an atmosphere object.

    throughputPlot: 
        Plots combined throughput given appropriate filter-keyed bandpass dictionary.

    filterPlot: 
        Plots the filter response curve from LSST filter data.

    hardwarePlot: 
        Plots the hardware response curve from LSST hardware data.

    phiPlot: 
        Plots normalized bandpass response function.

    dphiPlot: 
        Plots change in normalized bandpass response function given two filter-keyed 
        bandpass dictionaries.

    ddphiPlot: 
        Plots change in normalized bandpass response function given two filter-keyed bandpass 
        dictionaries and a standard filter-keyed bandpass dictionary.

    dmagPlot: 
        Given two filter-keyed bandpass dictionaries and a valid SED type, will plot dmags. 
    """
    def __init__(self):
        # List of strings containing component names
        self.components = ['H2O','O2','O3','Rayleigh','Aerosol']
        # List of strings containing components names to be used when plotting
        self.componentsPlot = ['$H_2O$','$O_2$','$O_3$','Rayleigh','Aerosol']
        # List of parameters (H2O,O2,O3,Rayleigh,Aerosol,Alpha)
        self.parameters = [1.0,1.0,1.0,1.0,1.0,1.7]
        # List of parameters used for plotting
        self.parametersPlot = [r'$t_{H_2O}$',r'$t_{O_2}$',r'$t_{O_3}$',r'$t_{Rayleigh}$',r'$t_{Aerosol}$',r'$t_{\alpha}$']
        # List of colors for used in plotting individual absorption components
        self.componentColors = {'H2O':'blue','O2':'green','O3':'red','Rayleigh':'purple','Aerosol':'cyan'}
        # Effective wavelength range, set in readModtranFiles
        self.wavelen = None
        # Min, max values of wavelength range
        self.wavelenRange = [WAVELENMIN,WAVELENMAX]
        # List of airmasses for which we have profiles, set in readModtranFiles
        self.airmasses = None
        # List of transmission profiles for individual airmasses
        self.transDict = None
        # Filter-keyed dictionary of filter S, set in readFilters
        self.filters = None
        # Filter-keyed dictionary of filter and hardware
        self.sys = None
        # List of filters
        self.filterlist = FILTERLIST
        # List of filter colors
        self.filtercolors = {'u': 'b', 'g': 'm', 'r': 'r','i': 'g','z': 'y','y4': 'k'}
        
        # SED data set with read functions:
        # Kurucz main sequence stars model data
        self.mss = None
        self.msList = None
        self.msMet = None
        self.msTemp = None
        self.msLogg = None

        # White dwarfs model data
        self.wds = None
        self.wdList = None
        self.wdListH = None
        self.wdListHe = None
        self.wdTemp = None
        self.wdLogg = None

        # Galaxy model data
        self.gals = None
        self.galList = None        
        self.galRedshifts = None

        # MLT dwarf model data
        self.mlts = None
        self.mltList = None
        self.mList = None
        self.lList = None
        self.tList = None

        # Quasar model data
        self.qsos = None
        self.qsoRedshifts = None

        # Supernova model data
        self.sns = None
        self.snList = None
        self.snDays = None
        self.snRedshifts = None
        
        # Readers
        self.readModtranFiles()
        print ''
        self.readHardware()

### Reading Functions
    
    def readModtranFiles(self, modtranDir='modtran/', modtranRoot='Pachon_MODTRAN', modtranSuffix='.7sc'):
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
        
        self.wavelen = np.arange(WAVELENMIN, WAVELENMAX+WAVELENSTEP, WAVELENSTEP, dtype='float')
        self.transDict = {}
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
                if (float(lineEle[0]) > self.wavelenRange[0]) | (float(lineEle[0]) > self.wavelenRange[1]):
                    wavelenTemp.append(float(lineEle[0]))
                    transTemp['H2O'].append(float(lineEle[2]))
                    transTemp['O2'].append(float(lineEle[3]))
                    transTemp['O3'].append(float(lineEle[4]))
                    transTemp['Rayleigh'].append(float(lineEle[8]))
                    #transTemp['Aerosol'].append(self.aerosol(float(lineEle[0]), airmass, STDAEROSOLALPHA, STDAEROSOLNORMCOEFF, STDAEROSOLNORMWAVELEN))
            fin.close()
            wavelenTemp = np.array(wavelenTemp, dtype='float')
            trans = {}
            templates = {}
            modtranComponents = ['H2O','O2','O3','Rayleigh']
            for comp in modtranComponents:
                trans[comp] = np.array(transTemp[comp], dtype='float')
                trans[comp] = np.interp(self.wavelen, wavelenTemp, transTemp[comp], left=0.0, right=0.0)
            self.airmasses.append(self._airmassToString(airmass))
            self.transDict[airmass] = copy.deepcopy(trans)
        if self.transDict != None:
            print "MODTRAN files have been read."
        return

    def readFilters(self, shiftPercent=None):
        """
        Reads LSST filter data and returns a filter-keyed dictionary. (S^{filters})

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        shiftPercent: (float) [None], percentage value to shift effective wavelength of filter throughput
        ----------------------

        * Modified from plot_dmags.py *
        """
        
        filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
        filters = {}

        for f in self.filterlist:
            filters[f] = Bandpass()
            filters[f].readThroughput(os.path.join(filterdir, "filter_" + f + ".dat"))
            effwavelenphi, effwavelensb = filters[f].calcEffWavelen()
            if shiftPercent != None:
                shift = effwavelensb * shiftPercent/100.0
                print f, shift
                filters[f].wavelen = filters[f].wavelen + shift
                filters[f].resampleBandpass()

        self.filters = filters

        print 'Read filter data from LSST software stack.'
        print 'Filters: ' + str(self.filterlist)

        return
    
    def readHardware(self, shiftPercent=None):
        """
        Reads LSST hardware data and returns a filter-keyed dictionary. (S^{sys})

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        shiftPercent: (float) [None], percentage value to shift effective wavelength of filter throughput
        ----------------------

        * Modified from plot_dmags.py *
        """
        
        filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
        hardware = ("detector.dat", "m1.dat", "m2.dat", "m3.dat", "lens1.dat", "lens2.dat", "lens3.dat")
        
        self.readFilters(shiftPercent=shiftPercent)

        filters = self.filters
        sys = {}

        for f in self.filterlist:
            sys[f] = Bandpass()
            tlist = []
            for t in hardware:
                tlist.append(os.path.join(filterdir, t))
            sys[f].readThroughputList(tlist)
            sys[f].wavelen, sys[f].sb = sys[f].multiplyThroughputs(filters[f].wavelen, filters[f].sb)

        self.sys = sys

        print 'Read hardware data from LSST software stack.'

        return

    def readMSs(self, sedDirectory=SIMSSEDLIBRARY, subDirectory='starSED/kurucz/'):
        """
        Reads Kurucz main sequence star model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.mss            # list of main sequence stars (SED objects)
        self.msList         # list of main sequence stars
        self.msMet          # list of main sequence star metallicities
        self.msTemp         # list of main sequence star temperatures
        self.msLogg         # list of main sequence star log surface gravities

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        sedDirectory: (string) [SIMSSEDLIBRARY], LSST SED library directory
        subDirectory: (string) ['starSED/kurucz/'], SED sub directory within sedDirectory
        ----------------------

        * Modified from plot_dmags.py *
        """
        ### Taken from plot_dmags and modified to suit specific needs.
        # read kurucz model MS, g40 stars SEDs
        homedir = os.getenv(sedDirectory)  
        stardir = os.path.join(homedir, subDirectory)
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
        temperature = np.array(atemperature)
        met = np.array(amet)
        logg = np.array(alogg)
        starlist = starlist2
        # actually read the stars SEDS from disk
        stars = {}
        for s in starlist:
            stars[s] = Sed()
            stars[s].readSED_flambda(os.path.join(stardir,s))
        print "# Read %d MS stars from %s" %(len(starlist), stardir)
        # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
        for s in starlist:
            stars[s].synchronizeSED(wavelen_min=WAVELENMIN, wavelen_max=WAVELENMAX, wavelen_step=WAVELENSTEP)

        self.mss = stars
        self.msList = starlist
        self.msTemp = temperature
        self.msMet = met
        self.msLogg = logg

        return

    def readWDs(self, sedDirectory=SIMSSEDLIBRARY, subDirectory='starSED/wDs/'):
        """
        Reads white dwarf model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.wds            # list of white dwarfs (SED objects)
        self.wdList         # list of white dwarfs
        self.wdListH        # list of hydrogen white dwarfs 
        self.wdListHe       # list of helium white dwarfs
        self.wdTemp         # list of white dwarf temperatures
        self.wdLogg         # list of white dwarf surface gravities

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        sedDirectory: (string) [SIMSSEDLIBRARY], LSST SED library directory
        subDirectory: (string) ['starSED/wDs/'], SED sub directory within sedDirectory
        ----------------------

        * Modified from plot_dmags.py *
        """
        # read white dwarf bergeron models
        homedir = os.getenv(sedDirectory)
        whitedwarfdir = os.path.join(homedir, subDirectory)
        
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

        temperatures = np.array(temperatures)
        loggs = np.array(loggs)
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
        for wd in wdlist:
            wds[wd].synchronizeSED(wavelen_min=WAVELENMIN, wavelen_max=WAVELENMAX, wavelen_step=WAVELENSTEP)

        self.wds = wds
        self.wdList = wdlist
        self.wdListH = hlist
        self.wdListHe = helist
        self.wdTemp = temperatures
        self.wdLogg = loggs

        return 

    def readGals(self, redshiftRange=[0,3.0], redshiftStep=0.5, sedDirectory=SIMSSEDLIBRARY, subDirectory='galaxySED/'):
        """
        Reads galaxy model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.gals           # list of galaxies (SED objects)
        self.galList        # list of galaxies      
        self.galRedshifts   # list of galaxy redshifts

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        redshiftRange: (list of floats) [0.0,3.0], range over which to redshift SEDs
        redshiftStep: (float) [0.5], redshift step
        sedDirectory: (string) [SIMSSEDLIBRARY], LSST SED library directory
        subDirectory: (string) ['galaxySED/'], SED sub directory within sedDirectory
        ----------------------

        * Modified from plot_dmags.py *
        """
        # read sn spectra and redshift
        homedir = os.getenv(sedDirectory)
        galdir = os.path.join(homedir, subDirectory)
        allfilelist = os.listdir(galdir)
        gallist_base = []
        metal = ['002Z', '04Z', '25Z']
        gtype = ['Const', 'Inst', 'Burst'] # removed Exp
        redshifts= np.arange(redshiftRange[0], redshiftRange[1]+redshiftStep, redshiftStep)
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
            gals[g].synchronizeSED(wavelen_min=WAVELENMIN, wavelen_max=WAVELENMAX, wavelen_step=WAVELENSTEP)
        # add dust
        ax, bx = gals[gallist[0]].setupCCMab()
        for g in gallist:
            gals[g].addCCMDust(ax, bx, A_v=0.02)

        self.gals = gals
        self.galList = gallist
        self.galRedshifts = redshifts 

        return 

    def readMLTs(self, sedDirectory=SIMSSEDLIBRARY, subDirectory='starSED/mlt/'):
        """
        Reads MLT dwarf model data from LSST software stack and sets relevant class attributes.

        The following attributes will be set:
        self.mlts           # list of mlt dwarfs (SED objects)
        self.mltList        # list of mlt dwarfs
        self.mList          # list of m dwarfs
        self.lList          # list of l dwarfs
        self.tList          # list of t dwarfs

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        sedDirectory: (string) [SIMSSEDLIBRARY], LSST SED library directory
        subDirectory: (string) ['starSED/mlt/'], SED sub directory within sedDirectory
        ----------------------

        * Modified from plot_dmags.py *
        """
        # read mlt stars - only keep 'm's
        # find the filenames and mark 'm', 'l', 't' stars separately
        homedir = os.getenv(sedDirectory)
        mltdir = os.path.join(homedir, subDirectory)
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
            mlts[s].synchronizeSED(wavelen_min=WAVELENMIN, wavelen_max=WAVELENMAX, wavelen_step=WAVELENSTEP)

        self.mlts = mlts
        self.mltList = mltlist
        self.mList = mlist
        self.lList = llist
        self.tList = tlist
        
        return 

    def readQsos(self, redshiftRange=[0,7.5], redshiftStep=0.1, sedDirectory='HOME', subDirectory='atmo2mags/seds/quasar/'):
        """
        Reads quasar model data and sets relevant class attributes.

        The following attributes will be set:
        self.qsos           # list of quasars (SED objects)    
        self.qsoRedshifts   # list of quasar redshifts

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        redshiftRange: (list of floats) [0.0,7.5], range over which to redshift SEDs
        redshiftStep: (float) [0.1], redshift step
        sedDirectory: (string) ['HOME'], LSST SED library directory
        subDirectory: (string) ['galaxySED/'], SED sub directory within sedDirectory
        ----------------------

        * Modified from plot_dmags.py *
        """
        homedir = os.getenv(sedDirectory)
        quasardir = os.path.join(homedir, subDirectory)
        # read zero redshift quasar
        base = Sed()
        base.readSED_flambda(os.path.join(quasardir, "quasar.dat"))
        # redshift 
        redshifts= np.arange(redshiftRange[0], redshiftRange[1]+redshiftStep, redshiftStep)
        quasars = {}
        for z in redshifts:
            wavelen, flambda = base.redshiftSED(z, wavelen=base.wavelen, flambda=base.flambda)
            quasars[z] = Sed(wavelen=wavelen, flambda=flambda)
        print "# Generated %d quasars at redshifts between %f and %f" %(len(redshifts), redshifts.min(), redshifts.max())
        # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
        for z in redshifts:
            quasars[z].synchronizeSED(wavelen_min=WAVELENMIN, wavelen_max=WAVELENMAX, wavelen_step=WAVELENSTEP)

        self.qsos = quasars
        self.qsoRedshifts = redshifts

        return

    def readSNs(self, redshiftRange=[0,1.2], redshiftStep=0.1, days=['0', '20', '40'], sedDirectory='HOME', 
        subDirectory='atmo2mags/seds/sn/'):
        """
        Reads supernova model data and sets relevant class attributes.

        The following attributes will be set:
        self.sns            # list of supernova (SED objects)
        self.snList         # list of supernova
        self.snDays         # list of supernova days
        self.snRedshifts    # list of supernova redshifts

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        redshiftRange: (list of floats) [0.0,1.2], range over which to redshift SEDs
        redshiftStep: (float) [0.1], redshift step
        days: (list of strings) ['0','20','40'], supernova days
        sedDirectory: (string) ['HOME'], LSST SED library directory
        subDirectory: (string) ['atmo2mags/seds/sn/'], SED sub directory within sedDirectory
        ----------------------

        * Modified from plot_dmags.py *
        """
        # read sn spectra and redshift
        homedir = os.getenv(sedDirectory)
        sndir = os.path.join(homedir, subDirectory)
        allfilelist = os.listdir(sndir)
        snlist = []
        redshifts= np.arange(redshiftRange[0], redshiftRange[1]+redshiftStep, redshiftStep)
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
            sns[s].synchronizeSED(wavelen_min=WAVELENMIN, wavelen_max=WAVELENMAX, wavelen_step=WAVELENSTEP)

        self.sns = sns
        self.snList = snlist
        self.snDays = days
        self.snRedshifts = redshifts
        
        return

    def readAll(self, mss=True, gals=True, wds=True, mlts=True, qsos=True, sns=True):
        """
        Read all (or subset of) SED model data.

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        mss: (boolean) [True], read kurucz main sequence star model data
        gals: (boolean) [True], read galaxy model data
        wds: (boolean) [True], read white dwarf model data
        mlts: (boolean) [True], read mlt dwarf model data
        qsos: (boolean) [True], read quasar model data
        sns: (boolean) [True], read supernova data
        ----------------------
        """
        if mss:
            self.readMSs()
        if wds:
            self.readWDs()
        if mlts:
            self.readMLTs()
        if gals:
            self.readGals()
        if qsos:
            self.readQsos()
        if sns:
            self.readSNs()

        return

### Calculator / Generator Functions
        
    def buildAtmo(self, P, X, aerosolNormCoeff=STDAEROSOLNORMCOEFF, aerosolNormWavelen=STDAEROSOLNORMWAVELEN):
        """
        Builds an atmospheric transmission profile (as an atmo object) given a set of component parameters 
        and an airmass. (S^{atm})

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        P: (list of floats), parameter array of t_k values for each atmospheric component
        X: (float), airmass
        aerosolNormCoeff: (float) [STDAEROSOLNORMCOEFF], aerosol power law normalization coefficient
        aerosolNormWavelen: (float) [STDAEROSOLNORMWAVELEN], aerosol normalization wavelength
        ----------------------
        """

        self._parameterCheck(P)
        self._airmassCheck(X)
        
        return Atmo(P, X, self.transDict, self.wavelen, aerosolNormCoeff, aerosolNormWavelen)
    
    def combineThroughputs(self, atmo, sys=None, filters=FILTERLIST):
        """
        Combines atmospheric transmission profile with system responsiveness data, returns filter-keyed 
        dictionary. (S^{atm}*S^{sys})

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        atmo: (atmo object), atmosphere
        sys: (dictionary) [None], filter-keyed bandpass dictionary of system responsiveness (if none
            will set to self.sys)
        filters: (list of strings) [FILTERLIST], list of filters
        ----------------------

        * Taken from plot_dmags and modified to suit specific needs. *
        """

        if filters == 'y4':
            filters = ['y4']

        if sys == None:
            sys = self.sys

        total = {}
        for f in filters:
            wavelen, sb = sys[f].multiplyThroughputs(atmo.wavelen, atmo.sb)
            total[f] = Bandpass(wavelen, sb, wavelen_min=WAVELENMIN, wavelen_max=WAVELENMAX, wavelen_step=WAVELENSTEP)
            total[f].sbTophi()
        return total
    
    def mags(self, bpDict, seds=None, sedkeylist=None, filters=FILTERLIST, verbose=False):
        """
        Calculates magnitudes given a bandpass dictionary, returns filter-keyed magnitude dictionary. If seds and sedkeylist are not none
        returns mags for Kurucz model MS stars.

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        bpDict: (dictionary), filter-keyed throughput dictionary
        seds: (list) [None], list of SED objects (if None, will set to main sequence stars' SED list)
        sedkeylist: (list) [None], list of SEDs (if None, will set to main sequence stars' SED key list)
        filters: (list of strings) [FILTERLIST], list of strings
        verbose: (boolean) [False], print out verbose statements
        ----------------------

        * Taken from plot_dmags and modified to suit specific needs. *
        """
        # calculate magnitudes for all sed objects using bpDict (a single bandpass dictionary keyed on filters)
        # pass the sedkeylist so you know what order the magnitudes are arranged in

        if filters == 'y4':
            filters = ['y4']

        if seds == None:
            seds = self.mss
        if sedkeylist == None:
            sedkeylist = self.msList

        mags = {}
        for f in filters:
            mags[f] = np.zeros(len(sedkeylist), dtype='float')
            for i,key in enumerate(sedkeylist):
                mags[f][i] = seds[key].calcMag(bpDict[f])
                if np.isnan(mags[f][i]):
                    print key, f, mags[f][i]
            if verbose == True:
                print f, mags[f].max(), mags[f].min()
        return mags
    
    def dmags(self, mags, mags_std, filters=FILTERLIST):
        """
        Returns filter-keyed dictionary of change in magnitude in millimagnitudes.

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        mags: (dictionary), filter-keyed dictionary of magnitudes
        mags_std: (dictionary), filter-keyed dictionary of magnitudes created at standard atmosphere
        filters: (list of strings) [FILTERLIST], list of filters
        ----------------------

        * Taken from plot_dmags and modified to suit specific needs. *
        """

        if filters == 'y4':
            filters = ['y4']

        dmags = {}
        for f in filters:
            # difference, in millimags
            dmags[f] = (mags[f] - mags_std[f]) * 1000.0
        return dmags
    
    def gi(self, mags_std):
        """
        Returns standard color temperature given standard magnitude dictionary keyed on filters.

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        mags_std: (dictionary), filter-keyed dictionary of magnitudes created at standard atmosphere
        ----------------------       

        * Taken from plot_dmags and modified to suit specific needs. *
        """
        # calculate some colors in the standard atmosphere, should be also standard bandpass, not shifted)
        gi = mags_std['g'] - mags_std['i']
        return gi

### Regression Functions

    def _computeLogL(self, P, X, err, f, dmags_obs, mags_std, seds, sedkeylist, deltaGrey):
        """Regression Function: returns log-likelihood given P, X, error, mags_obs, mags_std, seds, sedkeylist and deltaGrey."""

        atmo_fit = self.buildAtmo(P,X)
        throughput_fit = self.combineThroughputs(atmo_fit, filters=f)
        mags_fit = self.mags(throughput_fit, seds=seds, sedkeylist=sedkeylist, filters=f)
        dmags_fit = self.dmags(mags_fit, mags_std, filters=f)
        dmags_fit[f] -= deltaGrey
    
        return -np.sum(0.5 * ((dmags_fit[f] - dmags_obs[f]) / err) ** 2)

    def _computeChiSquared(self, f, dmags_fit, dmags_obs):
        """Returns array of chi squared values"""
        return np.sum((dmags_fit[f] - dmags_obs[f])**2 / dmags_obs[f])

    def computeAtmoFit(self, comp1, comp2, atmo_obs, err=0.005, componentBins=50, deltaGrey=0.0, deltaGreyBins=50, deltaGreyRange=[-50.0,50.0],
        regressionSed='mss', comparisonSeds=SEDTYPES, generateFig=True, generateDphi=True, saveLogL=True, useLogL=False, plotLogL=False, plotBoth=True,
        normalize=True, includeColorBar=False, pickleString='', filters=FILTERLIST, dmagLimit=True, returnData=False, verbose=True):
        """
        Computes the best fit atmospheric parameters for two given components and an observed atmosphere. Requires the 
        SED data for the specified regression and comparison SEDs to be read in. 

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        comp1: (string), name of component to regress
        comp2: (string), name of component to regress
        atmo_obs: (atmo object), observed atmosphere
        err: (float) [0.0005], percent error
        componentBins: (int) [50], number of bins for regression
        deltaGrey: (float) [0.0], adds extinction factor due to clouds (if less than 0 will subract mean dmags, 
            if greater than zero will subtract as mmag value from delta magnitudes during regression)
        deltaGreyBins: (int) [50], number of bins for regression over deltaGrey space
        deltaGreyRange: (list of ints), min and max deltaGrey value between which to regress
        regressionSed: (string) ['mss'], SED type to run regress over
        comparisonSeds: (list of strings) [SEDTYPES], 
        generateFig: (boolean) [True], generate a regression plot
        generateDphi: (boolean) [True], generate a dphi and ddphi plot
        saveLogL: (boolean) [True], save logL as txt file
        useLogL: (boolean) [False], use LogL to replace contour plots
        plotLogL: (boolean) [False], plot individual logLs and contours seperately
        plotBoth: (boolean) [False], plot both logLs and contours
        normalize: (boolean) [True], normalize logL by median when plotting
        includeColorBar: (boolean) [False], include logL color bar (requires useLogL to be True)
        pickleString: (string) [''], add custom string to plot titles
        filters: (list of strings) [FILTERLIST], list of filters
        dmagLimit: (boolean) [True], create +-2 mmags axis lines if certain axis requirements
            are met. 
        returnData: (boolean) [False], return data elements
        verbose: (boolean) [True], print out verbose statements
        ----------------------
        """
        # Insure valid parameters, airmass and sedtypes are given
        self._sedTypeCheck(regressionSed)

        # Find range over which to vary parameter and the parameter number for comp1, comp2
        range1, pNum1 = self._componentCheck(comp1,componentBins)
        range2, pNum2 = self._componentCheck(comp2,componentBins)

        # Find seds and sedkeylist for sedtype
        seds, sedkeylist = self._sedFinder(regressionSed)

        if verbose:
            print 'Computing nonlinear regression for ' + comp1 + ' and ' + comp2 + '.'
            print 'Observed atmosphere parameters: ' + str(atmo_obs.P)
            print 'Observed atmosphere airmass:    ' + str(atmo_obs.X)
            print 'Standard atmosphere parameters: ' + str(STDPARAMETERS)
            print 'Standard atmosphere airmass:    ' + str(STDAIRMASS)
            print 'Observed atmosphere parameter for ' + comp1 + ': ' + str(atmo_obs.P[pNum1])
            print 'Observed atmosphere parameter for ' + comp2 + ': ' + str(atmo_obs.P[pNum2])
            print ''
            print 'Fitting for %s between %.2f and %.2f in %s bins.' % (comp1, min(range2), max(range1), componentBins)
            print 'Fitting for %s between %.2f and %.2f in %s bins.' % (comp2, min(range2), max(range2), componentBins)
            print ''

            total = componentBins*componentBins

            if deltaGrey != 0.0:

                dgrange = np.linspace(deltaGreyRange[0], deltaGreyRange[1], deltaGreyBins)

                print 'Non-zero deltaGrey detected.'
                print 'Fitting for deltaGrey between %.2f and %.2f mmags in %s bins.' % (min(dgrange), max(dgrange), deltaGreyBins)
                print ''
                total = total*deltaGreyBins

            print 'Regressing %s parameter combinations per filter...' % (total)
            print ''
        
        P_fit = copy.deepcopy(atmo_obs.P)
        X_fit = copy.deepcopy(atmo_obs.X)

        # Create standard atmosphere and magnitudes
        std = self.buildAtmo(STDPARAMETERS,STDAIRMASS)
        throughput_std = self.combineThroughputs(std, filters=filters)
        mags_std = self.mags(throughput_std, seds=seds, sedkeylist=sedkeylist, filters=filters)

        # Create observed atmosphere and magnitudes
        throughput_obs = self.combineThroughputs(atmo_obs)
        mags_obs = self.mags(throughput_obs, seds=seds, sedkeylist=sedkeylist, filters=filters)
        dmags_obs = self.dmags(mags_obs, mags_std, filters=filters)

        logL = {}
        whr = {}
        comp1best = {}
        comp2best = {}
        dgbest = {}

        figName = self._regressionNameGen(comp1, comp2, atmo_obs, componentBins, err, regressionSed, deltaGrey, add=pickleString)

        for f in filters:

            dmags_obs[f] -= deltaGrey

            pickleString_temp = self._regressionNameGen(comp1, comp2, atmo_obs, componentBins, err, regressionSed, deltaGrey, add=pickleString, pickle=True, f=f)
                    
            print 'Calculating best fit parameters for ' + f + ' filter...'

            @pickle_results(os.path.join(PICKLEDIRECTORY, pickleString_temp))
            def run_regression(comp1, comp2, f):
                
                logL = []
                whr = []
                comp1best = []
                comp2best = []

                if deltaGrey != 0.0:
                    logL = np.ndarray([componentBins,componentBins,deltaGreyBins])

                    for d,dg in enumerate(dgrange):
                        for i in range(len(range1)):
                            for j in range(len(range2)):
                                P_fit[pNum1] = range1[i]
                                P_fit[pNum2] = range2[j]
                                logL[i][j][d] = self._computeLogL(P_fit, X_fit, err, f, dmags_obs, mags_std, seds, sedkeylist, dg)

                    logL -= np.amax(logL)
                    whr = np.where(logL == np.amax(logL))
                    comp1best = range1[whr[0][0]]
                    comp2best = range2[whr[1][0]]
                    dgbest = dgrange[whr[2][0]]
                else:
                    logL = np.ndarray([componentBins,componentBins,1])
                    for i in range(len(range1)):
                        for j in range(len(range2)):
                            P_fit[pNum1] = range1[i]
                            P_fit[pNum2] = range2[j]
                            logL[i][j][0] = self._computeLogL(P_fit, X_fit, err, f, dmags_obs, mags_std, seds, sedkeylist, deltaGrey)

                    logL -= np.amax(logL)
                    whr = np.where(logL == np.amax(logL))
                    comp1best = range1[whr[0][0]]
                    comp2best = range2[whr[1][0]]
                    dgbest = deltaGrey

                return comp1best, comp2best, dgbest, logL

            comp1best[f], comp2best[f], dgbest[f], logL[f]   = run_regression(comp1, comp2, f)

            print 'Completed ' + f + ' filter.'

            if saveLogL and deltaGrey != 0.0:
                name = self._regressionNameGen(comp1, comp2, atmo_obs, componentBins, err, regressionSed, deltaGrey, add=pickleString, f=f)
                np.savetxt(os.path.join(LOGLDIRECTORY, name + '_logL.txt'), logL[f][:][:][np.where(dgrange == dgbest[f])[0][0]])
                print 'Saved LogL at best fit deltaGrey for ' + f + ' filter.'
            elif saveLogL and deltaGrey == 0.0:
                name = self._regressionNameGen(comp1, comp2, atmo_obs, componentBins, err, regressionSed, deltaGrey, add=pickleString, f=f)
                np.savetxt(os.path.join(LOGLDIRECTORY, name + '_logL.txt'), logL[f][:][:][0])
                print 'Saved LogL for ' + f + ' filter.'

        if verbose and deltaGrey == 0.0:
            print ''
            print r'Best fit parameters (Filter, %s, %s):' % (comp1, comp2)
            for f in filters:
                print '%s %.2f %.2f' % (f, comp1best[f], comp2best[f])
        elif verbose and deltaGrey != 0.0:
            print ''
            print r'Best fit parameters (Filter, %s, %s, deltaGrey):' % (comp1, comp2)
            for f in filters:
                print '%s %.2f %.2f %.2f' % (f, comp1best[f], comp2best[f], dgbest[f])

        if generateDphi == True:

            throughput_fit = {}

            for f in filters:
                P_fit[pNum1] = comp1best[f]
                P_fit[pNum2] = comp2best[f]
                atmo_fit = self.buildAtmo(P_fit,X_fit)
                throughput_fit[f] = self.combineThroughputs(atmo_fit,filters=f)[f]

            self.dphiPlot(throughput_obs, throughput_std, bpDict2=throughput_fit, filters=filters, regression=True, figName=figName)
            self.ddphiPlot(throughput_obs, throughput_fit, throughput_std, filters=filters, regression=True, figName=figName)

        if generateFig == True:
            self.regressionPlot(comp1, comp1best, comp2, comp2best, dgbest, logL, atmo_obs, bins=componentBins, figName=figName, deltaGrey=deltaGrey,
                                regressionSed=regressionSed, comparisonSeds=comparisonSeds, useLogL=useLogL, dmagLimit=dmagLimit,
                                includeColorBar=includeColorBar, normalize=normalize, plotBoth=plotBoth, filters=filters, verbose=verbose)
        ###if returnData == True:
        return comp1best, comp2best, dgbest, logL
        ##else:
           # return

### Plotting Functions

    def regressionPlot(self, comp1, comp1_best, comp2, comp2_best, dgbest, logL, atmo_obs, bins=50, regressionSed='kurucz', comparisonSeds=SEDTYPES, plotDifference=True, 
        useLogL=False, includeColorBar=False, plotBoth=False, normalize=True, deltaGrey=0.0, dmagLimit=True, filters=FILTERLIST, verbose=True, figName=None,):
        """
        Plots regression data with each filter in its own row of subplots. Requires the 
        SED data for the specified regression and comparison SEDs to be read in.

        Recommendation: Use computeAtmoFit to call this function, will save a lot of work! 

        Parameters:
        ----------------------
        comp2_best: (dictionary), filter-keyed dictionary of best fit values
        parameter: (dtype) [default (if optional)], information

        comp1: (string), name of component to regress
        comp1_best: (dictionary), filter-keyed dictionary of best fit values
        comp2: (string), name of component to regress
        dgbest: (dictionary), filter-keyed dictionary of best fit deltaGrey values
        logL: (dictionary), filter-keyed dictionary of logL arrays
        atmo_obs: (atmo object), observed atmosphere
        bins: (int) [50], number of bins for regression
        regressionSed: (string) ['mss'], SED type to run regress over
        comparisonSeds: (list of strings) [SEDTYPES], 
        plotDifference: (boolean) [True], will plot difference between fit and truth for 
            comparisonSeds
        useLogL: (boolean) [False], use LogL to replace contour plots
        includeColorBar: (boolean) [False], include logL color bar (requires useLogL to be True)
        plotBoth: (boolean) [False], plot both logLs and contours
        normalize: (boolean) [True], normalize logL by median when plotting
        deltaGrey: (float) [0.0], adds extinction factor due to clouds (if less than 0 will subract mean dmags, 
            if greater than zero will subtract as mmag value from delta magnitudes during regression)
        dmagLimit: (boolean) [True], create +-2 mmags axis lines if certain axis requirements
            are met. 
        filters: (list of strings) [FILTERLIST], list of filters
        verbose: (boolean) [True], print out verbose statements
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """

        comp1_range, pNum1 = self._componentCheck(comp1, bins)
        comp2_range, pNum2 = self._componentCheck(comp2, bins)
        dgrange, dgnum = self._componentCheck('deltaGrey',bins)
            
        seds, sedkeylist = self._sedFinder(regressionSed)
        
        fig, ax = plt.subplots(len(filters),3)

        if deltaGrey == 0.0:
            if useLogL:
                title = r'$\Delta$mmags, Log-Likelihood, $\Delta\Delta$mmags for each LSST filter'
                fig.suptitle(title, fontsize=TITLESIZE)
            else:
                title = r'$\Delta$mmags, Regression Contours, $\Delta\Delta$mmags for each LSST filter'
                fig.suptitle(title, fontsize=TITLESIZE)
        else:
            if useLogL:
                title = r'$\Delta$mmags, Log-Likelihood, $\Delta\Delta$mmags for each LSST filter ($\delta$Grey: %s)' % (deltaGrey)
                fig.suptitle(title, fontsize=TITLESIZE)
            else:
                title = r'$\Delta$mmags, Regression Contours, $\Delta\Delta$mmags for each LSST filter ($\delta$Grey: %s)' % (deltaGrey)
                fig.suptitle(title, fontsize=TITLESIZE)

        fig.set_size_inches(15,len(filters)*5)
        fig.subplots_adjust(top=0.93, wspace=0.20, hspace=0.20, bottom=0.09, left=0.10, right=0.96)

        # Save observed parameters
        comp1_obs = atmo_obs.P[pNum1]
        comp2_obs = atmo_obs.P[pNum2]
    
        # Create observed throughput
        throughput_obs = self.combineThroughputs(atmo_obs)
        mags_obs = self.mags(throughput_obs, seds=seds, sedkeylist=sedkeylist, filters=filters)

        # Create standard atmosphere
        std = self.buildAtmo(STDPARAMETERS,STDAIRMASS)
        throughput_std = self.combineThroughputs(std)
        mags_std = self.mags(throughput_std, seds=seds, sedkeylist=sedkeylist, filters=filters)

        P_fit = copy.deepcopy(atmo_obs.P)
        X_fit = copy.deepcopy(atmo_obs.X)

        # For each filter plot dmags and regression contours
        for i,f in enumerate(filters):
            # Set component parameters to best fit parameters
            P_fit[pNum1] = comp1_best[f]
            P_fit[pNum2] = comp2_best[f]

            # Create atmosphere at best fit parameters
            fit = self.buildAtmo(P_fit,X_fit)
            throughput_fit = self.combineThroughputs(fit)

            self._dmagSED(ax[i][0], f, throughput_fit, throughput_std, regressionSed)
            self._dmagSED(ax[i][0], f, throughput_obs, throughput_std, regressionSed, truth=True)

            # Plot parameter space regression plots
            # Plot contours and true values
            if useLogL:
                self._logL(fig, ax[i][1], logL[f], 'imshow', comp1, comp1_obs, comp1_best[f], comp2, comp2_obs, 
                    comp2_best[f], deltaGrey, dgbest[f], bins, normalize=normalize, includeColorBar=includeColorBar)
            elif plotBoth:
                self._logL(fig, ax[i][1], logL[f], 'both', comp1, comp1_obs, comp1_best[f], comp2, comp2_obs, 
                    comp2_best[f], deltaGrey, dgbest[f], bins, normalize=normalize, includeColorBar=includeColorBar)
            else:
                self._logL(fig, ax[i][1], logL[f], 'contour', comp1, comp1_obs, comp1_best[f], comp2, comp2_obs, 
                    comp2_best[f], deltaGrey, dgrange[f], bins, normalize=normalize, includeColorBar=includeColorBar)

            # Plot dmags for other SEDS:
            if plotDifference == False:
                for s in comparisonSeds:
                    if s != regressionSed:
                        self._dmagSED(ax[i][2], f, throughput_fit, throughput_std, s, comparisonSed=True, dmagLimit=False)

                for s in comparisonSeds:
                    if s != regressionSed:
                        self._dmagSED(ax[i][2], f, throughput_obs, throughput_std, s, comparisonSed=True, dmagLimit=False, truth=True)
            else:
                for s in comparisonSeds:
                    if s != regressionSed:
                        self._dmagSED(ax[i][2], f, throughput_fit, throughput_std, s, comparisonSed=True, bpDict2=throughput_obs)

            if dmagLimit:
                self._axisLimiter(ax[i][0], [-2.0,2.0])
                self._axisLimiter(ax[i][2], [-2.0,2.0])

        label = self._sedLabelGen(regressionSed)
        col1Title = r'%s $\Delta$mmags' % (label)
        col2Title = 'Log-Likelihood'
        col3Title = r'$\Delta\Delta$mmags (Fit - Truth)'
        ax[0][0].set_title(col1Title, y=1.20, fontsize=LABELSIZE)
        ax[0][2].set_title(col3Title, y=1.20, fontsize=LABELSIZE)

        ax[0][0].legend(loc='upper center', bbox_to_anchor=(0.5,1.15), ncol=2, fontsize=LABELSIZE)
        if includeColorBar:
            ax[0][1].legend(loc='upper center', bbox_to_anchor=(0.5,1.27), ncol=2, fontsize=LABELSIZE)
            ax[0][1].set_title(col2Title, y=1.33, fontsize=LABELSIZE)
        else:
            ax[0][1].legend(loc='upper center', bbox_to_anchor=(0.5,1.15), ncol=2, fontsize=LABELSIZE)
            ax[0][1].set_title(col2Title, y=1.20, fontsize=LABELSIZE)
            
        ax[i][2].legend(loc='upper center', bbox_to_anchor=(-0.70,-0.2), ncol=len(comparisonSeds), 
            fontsize=LABELSIZE, title='Comparison SEDs')
            
        if figName != None:
            title = figName+"_regressionPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')

        return

    def _dmagSED(self, ax, f, bpDict1, bpDict_std, sedtype, bpDict2=None, truth=False, comparisonSed=False, dmagLimit=True):
        """Plots dmags for a specific filter to a given axis given appropriate filter-keyed bandpass dictionaries."""
        # Check if valid sedtype, check if sed data read:
        self._sedTypeCheck(sedtype)
        self._sedReadCheck(sedtype)

        # Label axes, only label y if not comparison sed
        ax.set_xlabel("g-i", fontsize=LABELSIZE)
        if comparisonSed == False:
            ax.set_ylabel(r"$\Delta$ %s (mmag)" %(f), fontsize=LABELSIZE)
        
        # Add grid
        ax.grid(b=True)

        if comparisonSed == True:
            label = self._sedLabelGen(sedtype)

        seds, sedkeylist = self._sedFinder(sedtype)

        if sedtype == 'mss':
            mags = self.mags(bpDict1, seds=seds, sedkeylist=sedkeylist)
            mags_std = self.mags(bpDict_std, seds=seds, sedkeylist=sedkeylist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=seds, sedkeylist=sedkeylist)
                dmags2 = self.dmags(mags2, mags_std)

            metallicity = np.array(self.msMet)
            logg = np.array(self.msLogg)
            metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
            metbinsize = abs(metallicity.min() - metallicity.max())/6.0
            metbins = np.arange(metallicity.min(), metallicity.max() + metbinsize, metbinsize)

            for metidx in range(len(metbins)):
                # Make cut of stars
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize) \
                        & (logg>3.5))
                mcolor = metcolors[metidx]

                if bpDict2 != None:
                    if metidx == len(metbins)-1:
                        ax.plot(gi[condition], dmags[f][condition]-dmags2[f][condition], mcolor+'.', label=label)
                    else:
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

        elif sedtype == 'qsos':
            mags = self.mags(bpDict1, seds=seds, sedkeylist=sedkeylist)
            mags_std = self.mags(bpDict_std, seds=seds, sedkeylist=sedkeylist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=seds, sedkeylist=sedkeylist)
                dmags2 = self.dmags(mags2, mags_std)

            redshift = self.qsoRedshifts
            redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
            redbinsize = 0.5
            redbins = np.arange(0.0, 3.0+redbinsize, redbinsize)
            for redidx in range(len(redbins)):
                condition =((redshift>=redbins[redidx]) & (redshift<=redbins[redidx]+redbinsize))
                rcolor = redcolors[redidx]

                if bpDict2 != None:
                    if redidx == len(redbins)-1:
                        ax.plot(gi[condition], dmags[f][condition]-dmags2[f][condition], rcolor+'o', label=label)
                    else:
                        ax.plot(gi[condition], dmags[f][condition]-dmags2[f][condition], rcolor+'o')
                else:
                    if truth == True:
                        if redidx == len(redbins)-1:
                            ax.plot(gi[condition], dmags[f][condition], rcolor+'o', label='Truth')
                        else:
                            ax.plot(gi[condition], dmags[f][condition], rcolor+'o')
                    else:
                        if redidx == len(redbins)-1:
                            ax.plot(gi[condition], dmags[f][condition], rcolor+'o', color='gray', label='Fit')
                        else: 
                            ax.plot(gi[condition], dmags[f][condition], rcolor+'o', color='gray')
        
        elif sedtype == 'gals':
            mags = self.mags(bpDict1, seds=seds, sedkeylist=sedkeylist)
            mags_std = self.mags(bpDict_std, seds=seds, sedkeylist=sedkeylist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=seds, sedkeylist=sedkeylist)
                dmags2 = self.dmags(mags2, mags_std)

            gallist = self.galList
            redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
            redbinsize = 0.5
            redbins = np.arange(0.0, 3.0+redbinsize, redbinsize)

            for i,g in enumerate(gallist):
                galbase, redshift = g.split('_')
                redshift = float(redshift)
                redidx = int(redshift / redbinsize)

                if bpDict2 != None:
                    if i == len(gallist)-1:
                        ax.plot(gi[i], dmags[f][i]-dmags2[f][i], redcolors[redidx]+'.', label=label)
                    else:
                        ax.plot(gi[i], dmags[f][i]-dmags2[f][i], redcolors[redidx]+'.')
                else:
                    if truth == True:
                        if i == len(gallist)-1:
                            ax.plot(gi[i], dmags[f][i], redcolors[redidx]+'.', label='Truth')
                        else:
                            ax.plot(gi[i], dmags[f][i], redcolors[redidx]+'.')
                    else:
                        if i == len(gallist)-1:
                            ax.plot(gi[i], dmags[f][i], redcolors[redidx]+'.', color='gray', label='Fit')
                        else:
                            ax.plot(gi[i], dmags[f][i], redcolors[redidx]+'.', color='gray')

        elif sedtype == 'mlts':
            mags = self.mags(bpDict1, seds=seds, sedkeylist=sedkeylist)
            mags_std = self.mags(bpDict_std, seds=seds, sedkeylist=sedkeylist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=seds, sedkeylist=sedkeylist)
                dmags2 = self.dmags(mags2, mags_std)

            mltlist = self.mltList
            mlist = self.mList
            llist = self.lList
            tlist = self.tList
        
            for j in range(len(mltlist)):
                if bpDict2 != None:
                    if j == len(mltlist)-1:
                        if (mltlist[j] in mlist):
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'bx', label=label)
                        elif (mltlist[j] in llist):
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'gx', label=label)
                        elif (mltlist[j] in tlist):
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'mx', label=label)
                    else:
                        if (mltlist[j] in mlist):
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'bx')
                        elif (mltlist[j] in llist):
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'gx')
                        elif (mltlist[j] in tlist):
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'mx')
                else:
                    if truth == True:
                        if j == len(mltlist)-1:
                            if (mltlist[j] in mlist):
                                ax.plot(gi[j], dmags[f][j], 'bx', label='Truth')
                            elif (mltlist[j] in llist):
                                ax.plot(gi[j], dmags[f][j], 'gx', label='Truth')
                            elif (mltlist[j] in tlist):
                                ax.plot(gi[j], dmags[f][j], 'mx', label='Truth')
                        else:
                            if (mltlist[j] in mlist):
                                ax.plot(gi[j], dmags[f][j], 'bx')
                            elif (mltlist[j] in llist):
                                ax.plot(gi[j], dmags[f][j], 'gx')
                            elif (mltlist[j] in tlist):
                                ax.plot(gi[j], dmags[f][j], 'mx')
                    else:
                        if j == len(mltlist)-1:
                            if (mltlist[j] in mlist):
                                ax.plot(gi[j], dmags[f][j], marker='x', color='gray', label='Fit')
                            elif (mltlist[j] in llist):
                                ax.plot(gi[j], dmags[f][j], marker='x', color='gray', label='Fit')
                            elif (mltlist[j] in tlist):
                                ax.plot(gi[j], dmags[f][j], marker='x', color='gray', label='Fit')
                        else:
                            if (mltlist[j] in mlist):
                                ax.plot(gi[j], dmags[f][j], marker='x', color='gray')
                            elif (mltlist[j] in llist):
                                ax.plot(gi[j], dmags[f][j], marker='x', color='gray')
                            elif (mltlist[j] in tlist):
                                ax.plot(gi[j], dmags[f][j], marker='x', color='gray')

        elif sedtype == 'wds':
            mags = self.mags(bpDict1, seds=seds, sedkeylist=sedkeylist)
            mags_std = self.mags(bpDict_std, seds=seds, sedkeylist=sedkeylist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=seds, sedkeylist=sedkeylist)
                dmags2 = self.dmags(mags2, mags_std)

            wdslist = self.wdList
            hlist = self.wdListH
            helist = self.wdListHe

            for j in range(len(wdslist)):
                if bpDict2 != None:
                    if (wdslist[j] in hlist):
                        if j == len(wdslist)-1:
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'y+', label=label)
                        else:
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'y+')
                    elif (wdslist[j] in helist):
                        if j == len(wdslist)-1:
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'y+', label=label)
                        else:
                            ax.plot(gi[j], dmags[f][j]-dmags2[f][j], 'y+')
                else:
                    if truth == True:
                        if (wdslist[j] in hlist):
                            if j == len(wdslist)-1:
                                ax.plot(gi[j], dmags[f][j], 'y+', label='Truth')
                            else:
                                ax.plot(gi[j], dmags[f][j], 'y+')
                        elif (wdslist[j] in helist):
                            if j == len(wdslist)-1:
                                ax.plot(gi[j], dmags[f][j], 'y+', label='Truth')
                            else:
                                ax.plot(gi[j], dmags[f][j], 'y+')
                    else:
                        if (wdslist[j] in hlist):
                            if j == len(wdslist)-1:
                                ax.plot(gi[j], dmags[f][j], marker='+', color='gray', label='Fit')
                            else:
                                ax.plot(gi[j], dmags[f][j], marker='+', color='gray')
                        elif (wdslist[j] in helist):
                            if j == len(wdslist)-1:
                                ax.plot(gi[j], dmags[f][j], marker='+', color='gray', label='Fit')
                            else:
                                ax.plot(gi[j], dmags[f][j], marker='+', color='gray')

        elif sedtype == 'sns':
            mags = self.mags(bpDict1, seds=seds, sedkeylist=sedkeylist)
            mags_std = self.mags(bpDict_std, seds=seds, sedkeylist=sedkeylist)
            gi = self.gi(mags_std)
            dmags = self.dmags(mags, mags_std)

            if bpDict2 != None:
                mags2 = self.mags(bpDict2, seds=seds, sedkeylist=sedkeylist)
                dmags2 = self.dmags(mags2, mags_std)
 
            snlist = self.snList
            redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
            redbinsize = 0.2
            redbins = np.arange(0.0, 1.2+redbinsize, redbinsize)
            day_symbol = {'0':'s', '20':'s', '40':'s'}

            for j,s in enumerate(snlist):
                day, redshift = s.split('_')
                redshift = float(redshift)
                redidx = int(redshift / redbinsize)

                if bpDict2 != None:
                    if j == len(snlist)-1:
                        ax.plot(gi[j], dmags[f][j]-dmags2[f][j], redcolors[redidx]+day_symbol[day], label=label)
                    else:
                        ax.plot(gi[j], dmags[f][j]-dmags2[f][j], redcolors[redidx]+day_symbol[day])
                else:
                    if truth == True:
                        if j == len(snlist)-1:
                            ax.plot(gi[j], dmags[f][j], redcolors[redidx]+day_symbol[day], label='Truth')
                        else:
                            ax.plot(gi[j], dmags[f][j], redcolors[redidx]+day_symbol[day])
                    else:
                        if j == len(snlist)-1:
                            ax.plot(gi[j], dmags[f][j], redcolors[redidx]+day_symbol[day], color='gray', label='Fit')
                        else:
                            ax.plot(gi[j], dmags[f][j], redcolors[redidx]+day_symbol[day], color='gray')

        return

    def transPlot(self, atmo1, atmo2=None, includeStdAtmo=False, includeComponents=False, figName=None):
        """
        Plots atmospheric transmission profile given an atmosphere object.
        
        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        atmo1: (atmo object), atmosphere
        atmo2: (atmo object) [None], optional comparison atmosphere
        includeStdAtmo: (boolean) [False], add comparison standard atmosphere
        includeComponents: (boolean) [False], include atmospheric components
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """ 
        fig,ax = plt.subplots(1,1)
        fig.set_size_inches(FIGUREWIDTH, FIGUREHEIGHT)
        
        ax.set_xlabel(r'Wavelength, $\lambda$ (nm)', fontsize=LABELSIZE)
        ax.set_ylabel(r'Transmission', fontsize=LABELSIZE)
        ax.set_title(r'$S^{atm}(\lambda)$', fontsize=TITLESIZE)
        
        if atmo2 != None:
            ax.plot(atmo2.wavelen,atmo2.sb, label=self._labelGen(atmo2.P, atmo2.X), alpha=0.5, color='black')
        elif includeStdAtmo:
            atmo2 = self.buildAtmo(STDPARAMETERS, STDAIRMASS)
            ax.plot(atmo2.wavelen, atmo2.sb, label=self._labelGen(atmo2.P, atmo2.X), alpha=0.5, color='black');
            ax.set_title(r'$S^{atm}(\lambda)$ and $S^{atm,std}(\lambda)$', fontsize=TITLESIZE)

        if includeComponents:
            ax.plot(atmo1.wavelen, atmo1.sb, color='black', label='Total', linestyle='--', alpha=0.7);
            for i,comp in enumerate(self.components):
                ax.plot(atmo1.wavelen,atmo1.sbDict[comp], color=self.componentColors[comp], label=self.componentsPlot[i])
                if atmo2 != None:
                    ax.plot(atmo2.wavelen,atmo2.sbDict[comp], alpha=0.5, color='black')
        else: 
            ax.plot(atmo1.wavelen, atmo1.sb, color='blue', label=self._labelGen(atmo1.P,atmo1.X));

        ax.legend(loc='lower right', shadow=False)
        
        if figName != None:
            title = figName + "_transPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')
        return

    def throughputPlot(self, bpDict1, bpDict2=None, includeStdAtmo=False, wavelenRange=[WAVELENMIN,WAVELENMAX], 
        filters=FILTERLIST, figName=None):
        """
        Plots combined throughput given appropriate filter-keyed bandpass dictionary.
        
        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        bpDict1: (dictionary), filter-keyed throughput dictionary
        bpDict2: (dictionary) [None], optional filter-keyed throughput dictionary
        includeStdAtmo: (boolean) [False], add standard atmosphere throughput dictionary
        wavelenRange: (list of ints) [WAVELENMIN, WAVELENMAX], wavelength plot range
        filters: (list of strings) [FILTERLIST], list of filters to plot
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """ 
        fig,ax = plt.subplots(1,1)
        fig.set_size_inches(FIGUREWIDTH, FIGUREHEIGHT)

        if includeStdAtmo:
            atmo_std = self.buildAtmo(STDPARAMETERS,STDAIRMASS)
            bpDict_std = self.combineThroughputs(atmo_std)

        for f in filters: 
            ax.plot(bpDict1[f].wavelen, bpDict1[f].sb, label=str(f), color=self.filtercolors[f])
            if includeStdAtmo:
                ax.plot(bpDict_std[f].wavelen, bpDict_std[f].sb, alpha=0.5, color='black')
            if bpDict2 != None:
                ax.plot(bpDict2[f].wavelen, bpDict2[f].sb, alpha=0.5, color='black')

        ax.set_title(r'$S^{atm}(\lambda)S_{b}^{sys}(\lambda)$', fontsize=TITLESIZE)
        if includeStdAtmo:
            ax.set_title(r'$S^{atm}(\lambda)S_{b}^{sys}(\lambda)$ and $S^{atm, std}(\lambda)S_{b}^{sys}(\lambda)$', fontsize=TITLESIZE)

        ax.set_xlabel(r'Wavelength, $\lambda$ (nm)', fontsize=LABELSIZE)
        ax.set_ylabel(r'Throughput', fontsize=LABELSIZE)
        ax.set_xlim(wavelenRange[0], wavelenRange[1]);
        ax.legend(loc='best')

        if figName != None:
            title = figName + "_throughputPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')

        return

    def filterPlot(self, filters=FILTERLIST, wavelenRange=[WAVELENMIN,WAVELENMAX], figName=None):
        """
        Plots the filter response curve from LSST filter data.
        
        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        filters: (list of strings) [FILTERLIST], list of filters to plot
        wavelenRange: (list of ints) [WAVELENMIN, WAVELENMAX], wavelength plot range     
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """ 
        
        fig,ax = plt.subplots(1,1)
        fig.set_size_inches(FIGUREWIDTH, FIGUREHEIGHT)
        
        for f in filters:
            ax.plot(self.filters[f].wavelen, self.filters[f].sb, label=str(f), color=self.filtercolors[f]);
        
        ax.set_xlim(wavelenRange[0], wavelenRange[1]);
        ax.set_ylim(0,1);
        ax.set_ylabel(r'Transmission', fontsize=LABELSIZE);
        ax.set_xlabel(r'Wavelength, $\lambda$ (nm)', fontsize=LABELSIZE);
        ax.set_title(r'LSST $S^{sys}_b$ Filters Only', fontsize=TITLESIZE);
        ax.legend(loc='best', shadow=False);

        if figName != None:
            title = figName + "_filterPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')

        return
    
    def hardwarePlot(self, filters=FILTERLIST, wavelenRange=[WAVELENMIN,WAVELENMAX], figName=None):
        """
        Plots the hardware response curve from LSST hardware data.
        
        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        filters: (list of strings) [FILTERLIST], list of filters to plot
        wavelenRange: (list of ints) [WAVELENMIN, WAVELENMAX], wavelength plot range     
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """ 

        fig,ax = plt.subplots(1,1)
        fig.set_size_inches(FIGUREWIDTH, FIGUREHEIGHT)
        
        for f in filters:
            ax.plot(self.sys[f].wavelen, self.sys[f].sb, label=str(f), color=self.filtercolors[f]);
        
        ax.set_xlim(wavelenRange[0], wavelenRange[1]);
        ax.set_ylim(0,1);
        ax.set_ylabel(r'Transmission', fontsize=LABELSIZE);
        ax.set_xlabel(r'Wavelength, $\lambda$ (nm)', fontsize=LABELSIZE);
        ax.set_title(r'LSST $S^{sys}_b$ Filters and Hardware', fontsize=TITLESIZE);
        ax.legend(loc='best', shadow=False);

        if figName != None:
            title = figName + "_hardwarePlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')
        return
    
    def phiPlot(self, bpDict1, bpDict2=None, filters=FILTERLIST, wavelenRange=[WAVELENMIN,WAVELENMAX], figName=None):
        """
        Plots normalized bandpass response function.

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        bpDict1: (dictionary), filter-keyed bandpass dictionary
        bpDict2: (dictionary) [None], optional comparison filter-keyed bandpass dictionary
        filters: (list of strings) [FILTERLIST], list of filters
        wavelenRange: (list of ints) [WAVELENMIN, WAVELENMAX], wavelength plot range     
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """
        
        fig,ax = plt.subplots(1,1)
        fig.set_size_inches(FIGUREWIDTH, FIGUREHEIGHT)
        
        for f in filters:
            ax.plot(bpDict1[f].wavelen, bpDict1[f].phi, label=str(f))
            if bpDict2 != None:
                ax.plot(bpDict2[f].wavelen, bpDict2[f].phi, alpha=0.5, color='black')
        
        ax.set_xlim(wavelenRange[0], wavelenRange[1]);
        ax.set_ylabel(r'$\phi_b^{obs}(\lambda)$', fontsize=LABELSIZE);
        ax.set_xlabel(r'Wavelength, $\lambda$ (nm)', fontsize=LABELSIZE);
        ax.set_title(r'Normalized Bandpass Response Function', fontsize=TITLESIZE);
        ax.legend(loc='best', shadow=False);
        
        if figName != None:
            title = figName + "_phiPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')
        return
    
    def dphiPlot(self, bpDict1, bpDict_std, bpDict2=None, filters=FILTERLIST, wavelenRange=[WAVELENMIN,WAVELENMAX], regression=False, figName=None):
        """
        Plots change in normalized bandpass response function given two filter-keyed bandpass dictionaries.
        
        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        bpDict1: (dictionary), filter-keyed bandpass dictionary
        bpDict_std: (dictionary), filter-keyed bandpass dictionary created at standard atmosphere
        bpDict2: (dictionary) [None], optional comparison filter-keyed bandpass dictionary
        filters: (list of strings) [FILTERLIST], list of filters
        wavelenRange: (list of ints) [WAVELENMIN, WAVELENMAX], wavelength plot range
        regression: (boolean) [False], if True will adjust y-labels for regression plotting 
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """

        fig,ax = plt.subplots(1,1)
        fig.set_size_inches(FIGUREWIDTH, FIGUREHEIGHT)

        for f in filters:
            ax.plot(self.wavelen, bpDict1[f].phi - bpDict_std[f].phi, color=self.filtercolors[f], label=str(f))
            if bpDict2 != None:
                ax.plot(self.wavelen, bpDict2[f].phi - bpDict_std[f].phi, color='gray', alpha=0.5)
        
        ax.set_xlim(wavelenRange[0], wavelenRange[1]);
        if regression:
            ax.set_ylabel(r"$\Delta\phi_b^{truth}(\lambda)$ and $\Delta\phi_b^{fit}(\lambda)$", fontsize=LABELSIZE);
        else:
            ax.set_ylabel(r"$\Delta\phi_b^{obs-std}(\lambda)$", fontsize=LABELSIZE);
        ax.set_xlabel(r"Wavelength, $\lambda$ (nm)", fontsize=LABELSIZE);
        ax.set_title(r"Change in Normalized Bandpass Response Function", fontsize=TITLESIZE);
        ax.legend(loc='best', shadow=False)
        
        if figName != None:
            title = figName + "_dphiPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')
        
        return

    def ddphiPlot(self, bpDict1, bpDict2, bpDict_std, filters=FILTERLIST, wavelenRange=[WAVELENMIN,WAVELENMAX], regression=False, figName=None):
        """
        Plots change in normalized bandpass response function given two filter-keyed bandpass dictionaries and a standard
        filter-keyed bandpass dictionary.

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        bpDict1: (dictionary), filter-keyed bandpass dictionary
        bpDict2: (dictionary), filter-keyed bandpass dictionary
        bpDict_std: (dictionary), filter-keyed bandpass dictionary created at standard atmosphere
        filters: (list of strings) [FILTERLIST], list of filters
        wavelenRange: (list of ints) [WAVELENMIN, WAVELENMAX], wavelength plot range    
        regression: (boolean) [False], if True will adjust y-labels for regression plotting 
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """

        fig,ax = plt.subplots(1,1)
        fig.set_size_inches(FIGUREWIDTH, FIGUREHEIGHT)

        for f in filters:
            ddphi = (bpDict2[f].phi - bpDict_std[f].phi) - (bpDict1[f].phi - bpDict_std[f].phi)
            ax.plot(self.wavelen, ddphi, color=self.filtercolors[f], label=str(f))

        ax.set_xlim(wavelenRange[0], wavelenRange[1]);
        if regression:
            ax.set_ylabel(r"$\Delta\phi_b^{fit}(\lambda) - \Delta\phi_b^{truth}(\lambda)$", fontsize=LABELSIZE);
        else:
            ax.set_ylabel(r"$\Delta\phi_b^{obs1-std}(\lambda) - \Delta\phi_b^{obs2-std}(\lambda)$", fontsize=LABELSIZE);
        ax.set_xlabel(r"Wavelength, $\lambda$ (nm)", fontsize=LABELSIZE);
        ax.set_title(r"Change in Normalized Bandpass Response Comparison", fontsize=TITLESIZE);
        ax.legend(loc='best', shadow=False)

        if figName != None:
            title = figName + "_ddphiPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')

        return

    def dmagPlot(self, bpDict1, bpDict_std, sedtype, filters=FILTERLIST, dmagLimit=True, figName=None):
        """
        Given two filter-keyed bandpass dictionaries and a valid SED type, will plot dmags. 

        Parameters:
        ----------------------
        parameter: (dtype) [default (if optional)], information

        bpDict1: (dictionary), filter-keyed bandpass dictionary
        bpDict_std: (dictionary), filter-keyed bandpass dictionary created at standard atmosphere
        sedtype: (string), name of SED type to plot
        filters: (list of strings) [FILTERLIST], list of filters
        dmagLimit: (boolean) [True], create +-2 mmags axis lines if certain axis requirements
            are met.   
        figName: (string) [None], if passed a string will save figure with string as title
        ----------------------
        """
        rows, columns = self._subplotFinder(filters)
        self._sedTypeCheck(sedtype)

        fig,ax = plt.subplots(rows, columns)
        fig.set_size_inches(10,len(filters)*2.5)
        fig.subplots_adjust(top=0.93, wspace=0.20, hspace=0.20, bottom=0.09, left=0.10, right=0.96)
        fig.suptitle(r'$\Delta$mmags for ' + self._sedLabelGen(sedtype), fontsize=TITLESIZE)

        filters = np.reshape(filters, (rows,columns))

        for i in range(rows):
            for j in range(columns):
                self._dmagSED(ax[i][j], filters[i][j], bpDict1, bpDict_std, sedtype, truth=True, dmagLimit=dmagLimit)
                if dmagLimit:
                    self._axisLimiter(ax[i][j],[-2.0,2.0])

        if figName != None:
            title = figName + "_dmagPlot.png"
            plt.savefig(os.path.join(PLOTDIRECTORY, title), format='png')

        return

    def _logL(self, fig, ax, logL, plotType, comp1, comp1_obs, comp1_best, comp2, comp2_obs, comp2_best, deltaGrey, dgbest, bins,
        normalize=True, includeColorBar=False):
        """Plots desired logL plot type given figure and axis object along with appropriate data."""
        comp1_range, pNum1 = self._componentCheck(comp1,bins)
        comp2_range, pNum2 = self._componentCheck(comp2,bins)
        dgrange, dgnum = self._componentCheck('deltaGrey', bins)

        if deltaGrey != 0.0:
            logL = logL[:][:][np.where(dgrange == dgbest)[0][0]]
        else:
            logL = logL[:][:][0]

        if normalize: 
            logL = logL / np.median(-logL)
        else:
            logL = logL

        if plotType == 'contour':
            contour = ax.contour(comp1_range, comp2_range, convert_to_stdev(logL.T), levels=(0.683, 0.955, 0.997), colors='k')
            ax.scatter(comp1_obs, comp2_obs, marker='o', s=25, facecolors='none', edgecolors='b', label='Truth')
            ax.clabel(contour, fontsize=9, inline=1)

        elif plotType == 'imshow':
            im = ax.imshow(logL.T, interpolation='nearest', cmap=plt.cm.bone, extent=(0.0,5.0,0.0,5.0), origin='lower')
            ax.scatter(comp1_obs, comp2_obs, marker='o', s=25, facecolors='none', edgecolors='b', label='Truth')
            if includeColorBar:
                fig.colorbar(im, ax=ax, format='%.0e')

        elif plotType == 'both':
            contour = ax.contour(comp1_range, comp2_range, convert_to_stdev(logL.T), levels=(0.683, 0.955, 0.997), colors='k')
            ax.scatter(comp1_obs, comp2_obs, marker='o', s=25, facecolors='none', edgecolors='b', label='Truth')
            ax.clabel(contour, fontsize=9, inline=1)
            im = ax.imshow(logL.T, interpolation='nearest', cmap=plt.cm.bone, extent=(0.0,5.0,0.0,5.0), origin='lower')
            if includeColorBar:
                fig.colorbar(im, ax=ax, format='%.0e')

        # Plot dashed lines at best fit parameters
        ax.axvline(comp1_best, color='black', linestyle='--', label='Fit')
        ax.axhline(comp2_best, color='black', linestyle='--')

        # Set y-axis, x-axis limits
        ax.set_xlim(min(comp1_range), max(comp1_range))
        ax.set_ylim(min(comp2_range), max(comp2_range))

        # Label axes
        if deltaGrey != 0.0:
            str1 = r'%s (fit: %.2f, truth: %.2f, $\delta$G: %.2f)' % (comp1, comp1_best, comp1_obs, dgbest)
            str2 = r'%s (fit: %.2f, truth: %.2f, $\delta$G: %.2f)' % (comp2, comp2_best, comp2_obs, dgbest)
        else:
            str1 = r'%s (fit: %.2f, truth: %.2f)' % (comp1, comp1_best, comp1_obs)
            str2 = r'%s (fit: %.2f, truth: %.2f)' % (comp2, comp2_best, comp2_obs)
        ax.set_xlabel(str1, fontsize=LABELSIZE)
        ax.set_ylabel(str2, fontsize=LABELSIZE)

        return


### Secondary Functions
    
    def _airmassToString(self, airmass):
        """Converts airmass to string"""
        X = float(airmass)
        return "%.3f" % (X)

    def _pToString(self, P):
        """Returns string version of parameter array."""
        stringP = "P"
        for i in P:
            if i < 1.0:
                stringP+="0"+str(int(i * 10))
            else:
                stringP+=str(int(i * 10))
        return stringP
    
    def _labelGen(self, P, X):
        """Generates label for use in plot legends."""
        label = []
        for paramNum,param in enumerate(P):
            name = self.parametersPlot[paramNum] + ':'
            labelEle = name + str(param)
            label.append(labelEle)
        return ' '.join(label) + ' $X$:' + str(X)

    def _sedLabelGen(self, sedtype):
        """Generates an appropriate SED label given a valid SED type."""
        if sedtype == 'mss':
            return 'Kurucz MS'
        elif sedtype == 'qsos':
            return 'Quasars'
        elif sedtype == 'gals':
            return 'Galaxies'
        elif sedtype == 'wds':
            return 'White Dwarfs'
        elif sedtype == 'mlts':
            return 'MLT Dwarfs'
        elif sedtype == 'sns':
            return 'Supernovas'
        return

    def _figNameGen(self, saveFig, figName, P1, X1, P2, X2):
        """Generates a string for figure names: 'X1_P1_X2_P2_figName' """
        if saveFig == True:
            if figName != None:
                figName='X'+str(int(X1*10))+'_'+self._pToString(P1)+'_'+'X'+str(int(X2*10))+'_'+self._pToString(P2)+'_'+figName
            else:
                figName='X'+str(int(X1*10))+'_'+self._pToString(P1)+'_'+'X'+str(int(X2*10))+'_'+self._pToString(P2)
        else:
            figName = None
        return figName

    def _regressionNameGen(self, comp1, comp2, atmo, bins, err, regressionSed, deltaGrey, add='', 
        pickle=False, f=None):
        """Generates a string for pickle files. """
        X_obs = 'X' + str(int(atmo.X*10))
        P_obs = self._pToString(atmo.P)
        comps = comp1 + '_' + comp2
        X_std = 'XSTD' + str(int(STDAIRMASS*10))
        DG = 'DG' + str(int(deltaGrey*10.0))
        ERR = 'E' + str(int(err*1000))  
        REG = ''
        bins = str(bins) + 'b'
        ext = ''

        if f != None:
            REG = regressionSed + '_' + f
        else:
            REG = regressionSed 

        if add != '':
            add = '_' + add 

        if pickle == True:
            ext = add + '.pkl'
        else:
            ext = add

        return '%s_%s_%s_%s_%s_%s_%s_%s%s' % (X_obs, P_obs, comps, X_std, DG, ERR, REG, bins, ext)

    def _componentCheck(self, comp, bins):
        """Returns a range of values of length bins for a given component."""
        if comp == 'H2O':
            return np.linspace(0.2,5.0,bins), 0
        elif comp == 'O2':
            return np.linspace(0.2,5.0,bins), 1
        elif comp == 'O3':
            return np.linspace(0.2,5.0,bins), 2
        elif comp == 'Rayleigh':
            return np.linspace(0.2,5.0,bins), 3
        elif comp == 'Aerosol':
            return np.linspace(0.2,5.0,bins), 4
        elif comp == 'Alpha':
            return np.linspace(0.2,5.0,bins), 5
        else:
            raise ValueError(comp + ' is not a valid component')

    def _parameterCheck(self, P):
        """Checks if parameter array is of appropriate length."""
        if len(P) != 6:
            raise ValueError('Need 6 parameters to build atmosphere!')
        return

    def _airmassCheck(self, X):
        """Checks if airmass is valid (one that has data)."""
        if self._airmassToString(X) not in self.airmasses:
            raise ValueError('Not a valid airmass, check self.airmasses for valid airmasses')
        return

    def _sedTypeCheck(self, sedtype):
        """Checks if SED type is valid."""
        if sedtype not in SEDTYPES:
            raise ValueError(str(sedtype) + ' is not a valid SED type, valid SED types: ' + str(SEDTYPES))
        return
    
    def _sedReadCheck(self, sedtype):
        """Checks if sed model data has been read in."""
        if sedtype == 'mss':
            if self.mss == None:
                raise ValueError('No Kurucz model data found, please run self.readMSs() or self.readAll()')
        elif sedtype == 'qsos':
            if self.qsos == None:
                raise ValueError('No quasar data found, please run self.readQsos() or self.readAll()')
        elif sedtype == 'gals':
            if self.gals == None:
                raise ValueError('No galaxy data found, please run self.readGals() or self.readAll()')
        elif sedtype == 'wds':
            if self.wds == None:
                raise ValueError('No white dwarf data found, please run self.readWDs() or self.readAll()')
        elif sedtype == 'mlts':
            if self.mlts == None:
                raise ValueError('No mlt dwarf data found, please run self.readMLTs() or self.readAll()')
        elif sedtype == 'sns':
            if self.sns == None:
                raise ValueError('No supernova data found, please run self.readSNs() or self.readAll()')
        return

    def _colorCheck(self, color, mags_std):
        """Checks if given color is valid, if valid returns the color given standard magnitudes."""
        if color in COLORS:
            if color == 'g-i':
                return color, self.gi(mags_std)
            if color == 'u-g':
                return color, self.ug(mags_std)
            if color == 'g-r':
                return color, self.gr(mags_std)
            if color == 'r-i':
                return color, self.ri(mags_std)
            if color == 'i-z':
                return color, self.iz(mags_std)
            if color == 'z-y' or color == 'z-y4':
                return 'z-y', self.zy(mags_std)
        else:
            raise ValueError('Please choose a valid color from ' + str(COLORS))
        return

    def _sedFinder(self, sedtype):
        """Returns seds and sedkeylist given an sedtype."""
        seds = []
        sedkeylist = []

        if sedtype == 'mss':
            seds = self.mss
            sedkeylist = self.msList
        elif sedtype == 'qsos':
            seds = self.qsos
            sedkeylist = self.qsoRedshifts
        elif sedtype == 'gals':
            seds = self.gals
            sedkeylist = self.galList
        elif sedtype == 'wds':
            seds = self.wds
            sedkeylist = self.wdList
        elif sedtype == 'mlts':
            seds = self.mlts
            sedkeylist = self.mltList
        elif sedtype == 'sns':
            seds = self.sns
            sedkeylist = self.snList

        return seds, sedkeylist

    def _subplotFinder(self, filters):
        """Returns rows and columns for a given filter list."""
        subplots = len(filters)
        if subplots == 1:
            return 1, 1
        elif subplots % 2 == 0:
            return subplots / 2, 2
        elif subplots == 3:
            return 3, 1
        else:
            return 3, 2
        return

    def _axisLimiter(self, ax, limits):
        """Sets appropriate axis limits given an axis and threshold limits."""
        min_y, max_y = ax.get_ylim()

        if max_y < limits[0]:
            max_y = limits[0]
        elif min_y > limits[1]:
            min_y = limits[1]
        elif min_y < limits[0] and max_y < limits[1]:
            ax.axhline(limits[0], color='black', linestyle='--')
        elif min_y > limits[0] and max_y > limits[1]:
            ax.axhline(limits[1], color='black', linestyle='--')
        elif min_y < limits[0] and max_y > limits[1]:
            ax.axhline(limits[0], color='black', linestyle='--')
            ax.axhline(limits[1], color='black', linestyle='--')
        else:
            min_y = limits[0]
            max_y = limits[1]

        ax.set_ylim(min_y,max_y)
        return
