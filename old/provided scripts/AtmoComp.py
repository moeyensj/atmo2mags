"""
AtmoComp -

$Id$

ljones@astro.washington.edu

A class designed to read modtran atmospheric component files, investigate various components of the atmosphere and how
they scale with airmass, and also build atmospheres with variable components.

Class data:
 atmo_ind (list of the individual atmospheric absorption features)
 atmo_templates (templates of the individual atmospheric absorption features) - dictionary keyed airmass / comp
 wavelength 
 total transmission
 
"""

import os
import copy
import numpy
import pylab


deg2rad = numpy.pi/180.0
rad2deg = 180.0/numpy.pi
colors = ('b', 'm', 'r', 'g', 'y', 'k')
interplimit = 1e-3


class AtmoComp:
    """Class for reading modtran component files and building an atmosphere."""    
    def __init__(self):
        # lists of the atmospheric components
        self.atmo_comp = ('comb', 'H2O', 'O2', 'O3', 'rayleigh', 'aerosol')
        self.atmo_ind = ('H2O', 'O2', 'O3', 'rayleigh', 'aerosol')
        # list of the airmasses we have templates for
        self.seczlist=None
        # wavelength range of templates
        self.wavelen = None
        # save the transmission templates, mostly for the 'comb' total transmission
        self.atmo_trans = None
        # the templates themselves
        self.atmo_templates = None
        # the airmass we are currently interested in
        self.secz=-99
        # the atmosphere absorption templates scaled to this airmass
        self.atmo_abs = None
        # the atmospheric template coefficients (defaults)
        self.setCoefficients()
        # the total atmospheric absorption at this airmass and with these coefficients
        self.trans_total = None
        # why not go ahead and read the modtran files, for now
        self.readModtranFiles()
        return
    
    
    def readModtranFiles(self, modtranDir=".", modtranRoot="Pachon_MODTRAN", modtranSuffix=".7sc"):    
        """Read the modtran output files, keep the atmospheric components in a dictionary."""
        # find which files are appropriate modtran files
        filenames = os.listdir(".")
        atmofiles = []
        for f in filenames:
            if (f.startswith(modtranRoot)) & (f.endswith(modtranSuffix)):
                atmofiles.append(f)
        # wavelength in nanometers
        self.wavelen = numpy.arange(300.0, 1100.5, 0.5, dtype='float')
        self.atmo_templates = {}
        self.atmo_trans = {}
        self.seczlist = []
        for file in atmofiles:
            # read modtran files from David Burke
            # modtran file names are like Pachon_MODTRAN.12.7sc  .. 12=X=1.2
            # files have 2 header lines which include unknown info
            # then 2 lines of info on columns
            #  WAVELENGTH COMBIN    H2O   UMIX     O3  TRACE     N2    H2O MOLEC AER+CLD  HNO3 AER+CLD
            #   (NM)  TRANS  TRANS  TRANS  TRANS  TRANS   CONT   CONT   SCAT  TRANS  TRANS abTRNS
            # then continue until -9999 flag for wavelength at end of file
            # 
            # Translation from MODTRAN to ATMO_fit:
            #   H2O = H2O TRANS * H2O CONT = water
            #   O2  = UMIX * TRACE  = molecular absorption
            #   O3  = O3 = ozone
            #   Rayleigh = MOLEC SCAT = molecular scattering
            #   aerosol = Use a built-in power-law model. (but could use AER+CLD trans* abTrans it looks like)
            #        (although in these sims, =1.0 all the way through)
            # comb = combined version (i.e. absorption at that airmass for these components calculated by modtran)
            #            
            fin = open(file, 'r')   
            rwavelen = []
            rtrans = {}
            for comp in self.atmo_comp:
                rtrans[comp] = []
            i = 0
            # read data from file
            for lines in fin:
                if i<4: # first four lines are header info
                    if i == 1:
                        values = lines.split()
                        secz = 1/numpy.cos((float(values[2]))*deg2rad)
                        secz = round(secz*10)/10.0
                    i = i+1
                    continue # skip over header lines
                values = lines.split()
                if (float(values[0]) < 0):
                    break  # end of file marker (-9999.)
                rwavelen.append(values[0])
                rtrans['comb'].append(values[1])
                rtrans['H2O'].append(float(values[2])*float(values[7]))
                rtrans['O2'].append(float(values[3])*float(values[5]))
                rtrans['O3'].append(values[4])
                rtrans['rayleigh'].append(values[8])
                rtrans['aerosol'].append(float(values[9])*float(values[11]))
            fin.close()
            # translate data into numpy arrays
            rwavelen = numpy.array(rwavelen, dtype='float')
            trans = {}
            templates = {}
            for comp in self.atmo_comp:
                rtrans[comp] = numpy.array(rtrans[comp], dtype=float)
                # resample transmission onto the same wavelength grid
                trans[comp] = numpy.interp(self.wavelen, rwavelen, rtrans[comp], left=0.0, right=0.0)
                templates[comp] = 1.0 - trans[comp]
            # fix aerosols to somewhat realistic values
            trans['aerosol'] = numpy.exp(-(0.039 * secz)*(self.wavelen/675.0)**-1.7)
            templates['aerosol'] = 1.0 - trans['aerosol']
            # finalize individual airmass
            z = self.seczToString(secz)
            self.seczlist.append(z)
            self.atmo_templates[z] = copy.deepcopy(templates)
            self.atmo_trans[z] = copy.deepcopy(trans)
        # done with individual airmass
        return

    def seczToString(self, secz):
        z = float(secz)        
        return "%.3f" %(z)

    def setCoefficients(self, t0=(3.9/100.0), t1=(0.02/100.0), t2=(-0.03/100.0), alpha=-1.7,
                        mol=0.96, BP=782, O3=0.9, H2O=0.9):
        """Set the atmospheric template coefficients."""
        self.C = {'t0':t0, 't1':t1, 't2':t2, 'alpha':alpha, 'mol':mol, 'BP':BP, 'O3':O3, 'H2O':H2O}
        # the atmospheric template coefficients (defaults)
        #C = {'to':3.9/100.0, 't1':0.02/100.0, 't2':-0.03/100.0, 'alpha':-1.7, 
        #     'mol':0.96, 'BP':782, 'O3':0.9, 'H2O':0.9}
        return 
        
    def interpolateSecz(self, secz):
        """Generate atmospheric absorption profiles for this particular secz by linear interpolation."""
        if abs(secz - self.secz) < interplimit:
            print "Already have calculated atmospheric absorption profiles for this airmass (or close enough) so will use that."
            return        
        # else reset current values
        self.secz = secz
        self.atmo_abs = None
        self.total_trans = None
        if self.seczToString(secz) in self.seczlist:
            print "This is the same as one of the modtran profiles, so I'll use that."
            self.atmo_abs = copy.deepcopy(self.atmo_templates[self.seczToString(secz)])
            self.total_trans = None
            return
        # otherwise, go find closest airmasses which are in atmo_templates
        zlist = numpy.array(self.seczlist, dtype='float')
        condition = ((zlist-secz)<=0)
        seczlo = zlist[condition].max()
        zlo = self.seczToString(seczlo)
        condition = ((zlist-secz)>0)
        seczhi = zlist[condition].min()
        zhi = self.seczToString(seczhi)        
        self.atmo_abs = {}
        ratio = (secz - seczlo) / (seczhi - seczlo)
        # and interpolate
        for comp in self.atmo_comp:
            self.atmo_abs[comp] = (1-ratio)*self.atmo_templates[zlo][comp] + ratio*self.atmo_templates[zhi][comp]
        # reset trans_total to None, as not calculated yet
        self.trans_total = None
        return                      

    def buildAtmos(self, secz, xlim=[300, 1100], doPlot=False):
        """Generate the total atmospheric transmission profile at this airmass, using the coefficients C."""
        # Burke paper says atmosphere put together as 
        # Trans_total (alt/az/time) = Tgray * (e^-Z*tau_aerosol(alt/az/t)) * 
        #         * (1 - C_mol * BP(t)/BPo * A_mol(Z))  -- line 2
        #         * (1 - sqrt(C_mol * BP(t)/BPo) * A_mol(Z))  -- 3
        #         * (1 - C_O3 * A_O3(A) )
        #         * (1 - C_H2O(alt/az/time) * A_H2O(Z))
        # Tau_aerosol = trans['aerosol'] ... but replace with power law (because here all 1's)
        #  typical power law index is about tau ~ lambda^-1
        # A_mol = trans['O2']
                
        # secz = secz of this observation
        # wavelen / atmo_templates == building blocks of atmosphere, with seczlist / atmo_ind keys
        # C = coeffsdictionary = to, t1, t2, alpha0 (for aerosol), C_mol, BP, C_O3, C_H2O  values    
        if (abs(secz - self.secz) > interplimit):
            print "Generating interpolated atmospheric absorption profiles for this airmass %f" %(secz)
            self.interpolateSecz(secz)

        BP0 = 782 # mb
        # set aerosol appropriately with these coefficients
        self.atmo_abs['aerosol'] = 1.0 - numpy.exp(-secz * (self.C['t0'] + self.C['t1']*0.0 + self.C['t2']*0.0) 
                                                   * (self.wavelen/675.0)**self.C['alpha'])
        # set total transmission, with appropriate coefficients
        self.trans_total = numpy.ones(len(self.wavelen), dtype='float')
        self.trans_total = self.trans_total * (1.0 - self.C['mol'] * self.C['BP']/BP0 * self.atmo_abs['rayleigh'])  \
                      * ( 1 - numpy.sqrt(self.C['mol'] * self.C['BP']/BP0) * self.atmo_abs['O2']) \
                      * ( 1 - self.C['O3'] * self.atmo_abs['O3']) \
                      * ( 1 - self.C['H2O'] * self.atmo_abs['H2O']) \
                      * ( 1 - self.atmo_abs['aerosol'])
        # now we can plot the atmosphere
        if doPlot:
            pylab.figure()
            pylab.subplot(212)
            colorindex = 0
            for comp in self.atmo_ind:
                pylab.plot(self.wavelen, self.atmo_abs[comp], colors[colorindex], label='%s' %(comp))
                colorindex = self._next_color(colorindex)        
            leg =pylab.legend(loc=(0.88, 0.3), fancybox=True, numpoints=1, shadow=True)
            ltext = leg.get_texts()
            pylab.setp(ltext, fontsize='small')
            coefflabel = ""
            for comp in ('mol', 't0', 'alpha', 'O3', 'H2O'):
                coefflabel = coefflabel + "C[%s]:%.2f  " %(comp, self.C[comp])
                if (comp=='alpha') | (comp=='mol'):
                    coefflabel = coefflabel + "\n"
            pylab.figtext(0.2, 0.35, coefflabel, fontsize='small')
            pylab.xlim(xlim[0], xlim[1])
            pylab.ylim(0, 1.0)
            pylab.xlabel("Wavelength (nm)")
            pylab.subplot(211)
            pylab.plot(self.wavelen, self.atmo_trans[self.seczToString(1.2)]['comb'], 'r-', label='Standard X=1.2 (no aerosols)')
            pylab.plot(self.wavelen, self.trans_total, 'k-', label='Observed')
            leg = pylab.legend(loc=(0.12, 0.05), fancybox=True, numpoints=1, shadow=True)
            ltext = leg.get_texts()
            pylab.setp(ltext, fontsize='small')
            pylab.xlim(xlim[0], xlim[1])
            pylab.ylim(0, 1.0)
            pylab.title("Example Atmosphere at X=%.2f" %(secz))
        return


    def plotTrans(self, secz, xlim=(300, 1100), newfig=True, savefig=False, figroot='atmos'):
        """Plot atmospheric transmission for all components."""
        if (abs(secz - self.secz) > interplimit):
            print "Generating interpolated atmospheric absorption profiles for this airmass %f" %(secz)
            self.interpolateSecz(secz)
        if newfig:
            pylab.figure()
        colorindex = 0
        for comp in self.atmo_ind:
            pylab.plot(self.wavelen, (1.0-self.atmo_abs[comp]), colors[colorindex], label=comp)
            colorindex = self._next_color(colorindex)
        if self.trans_total != None:
            pylab.plot(self.wavelen, self.trans_total, 'k:')
        leg = pylab.legend(loc='lower right', numpoints=1, fancybox=True)
        ltext = leg.get_texts()
        pylab.setp(ltext, fontsize='small')
        pylab.xlim(xlim[0], xlim[1])
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Transmission")
        pylab.title("Airmass %.2f" %(secz))
        return

    def plotTemplates(self, secz, xlim=(300, 1100), newfig=True, savefig=False, figroot='atmos'):
        """Plot atmospheric absorption templates (not scaled with C) for all components."""
        if (abs(secz - self.secz) > interplimit):            
            print "Generating interpolated atmospheric absorption profiles for this airmass %f" %(secz)
            self.interpolateSecz(secz)
        if newfig:
            pylab.figure()
        colorindex = 0
        for comp in self.atmo_ind:
            pylab.plot(self.wavelen, self.atmo_abs[comp], colors[colorindex], label=comp)
            colorindex = self._next_color(colorindex)
        leg = pylab.legend(loc=(0.85, 0.7), numpoints=1, fancybox=True)
        ltext = leg.get_texts()
        pylab.setp(ltext, fontsize='small')
        pylab.xlim(xlim[0], xlim[1])
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Template")
        pylab.title("Airmass %.2f" %(secz))
        return


    def plotAtmosRatio(self, this_seczlist, xlim=(300, 1100),
                       newfig=True, savefig=False, figroot='atmos_ratio'):
        trans_ratio = numpy.zeros(len(self.wavelen), dtype='float')
        zref = self.seczToString(this_seczlist[0])
        eps = 1e-30
        for comp in self.atmo_comp: 
            tmp = numpy.where(self.atm_trans[zref]==0, eps, self.atm_trans[zref])
        pylab.figure()        
        for z in this_seczlist:            
            trans_ratio = self.atm_trans[self.seczToString(z)][comp] / self.atm_trans[zref][comp]
            pylab.plot(self.wavelen, trans_ratio, label='X=%.2f' %(float(z)))
        pylab.legend(loc='lower right', fancybox=True, numpoints=1)
        pylab.xlim(xlim[0], xlim[1])
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Ratio")
        pylab.title("Changes in transmission with airmass due to component %s" %(comp))
        return

    def plotAtm(self, this_seczlist, xlim=(300, 1100),
                newfig=True, savefig=False, figroot='atmos_piece'):
        for comp in self.atmo_ind: 
            pylab.figure()
            for z in this_seczlist:
                pylab.plot(self.wavelen, self.atmo_trans[z][comp], label='X=%.2f' %(float(z)))
            pylab.legend(loc='lower right', fancybox=True, numpoints=1)
            pylab.xlim(xlim[0], xlim[1])
            pylab.ylim(0, 1)
            pylab.xlabel("Wavelength (nm)")
            pylab.ylabel("Transmission")
            pylab.title("Transmission of component %s" %(comp))
        return    
        
    def plotAbs(self, this_seczlist, xlim=(300, 1100),
                newfig=True, savefig=False, figroot='atmos_template'):
        for comp in self.atmo_ind: 
            pylab.figure()
            for z in this_seczlist:
                pylab.plot(self.wavelen, self.atmo_templates[self.seczToString(z)][comp], label='X=%.2f' %(float(z)))
            pylab.legend(loc=(0.86, 0.7), fancybox=True, numpoints=1)
            pylab.xlim(xlim[0], xlim[1])
            pylab.ylim(0, 0.62)
            pylab.xlabel("Wavelength (nm)")
            pylab.ylabel("Absorption")
            pylab.title("Absorption due to component %s" %(comp))
        return

    def writeAtmo(self, filename):
        if self.trans_total == None:
            raise Exception("No trans_total defined yet. Build an atmosphere first.")
        f = open(filename, 'w')
        print >>f, "# Test atmosphere "
        print >>f, "# Airmass = %.3f" %(self.secz)
        for coeff in self.C.keys():
            print >>f, "# Coefficient: C['%s'] = %f" %(coeff, self.C[coeff])
        for i in range(len(self.wavelen)):
            print >>f, "%.3f  %.8f" %(self.wavelen[i], self.trans_total[i])
        f.close()
        return

    def _next_color(self, colorindex):
        colorindex = colorindex + 1
        if colorindex > len(colors):
            colorindex = 0
        return colorindex


