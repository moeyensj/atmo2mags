#### README ####

The AtmoBuilder class is a powerful and continuously growing tool aimed at providing the functionality required to 
constrain the atmosphere and its atmospheric absorption parameters by using known stellar magnitudes. 

Required Software:
  - LSST Software Stack (including prerequisites)
  - AstroML

Installation:
  - Clone this repository: $ git clone https://github.com/moeyensj/atmo2mags.git
  - Setup the LSST enviroment ($ source ~/lsst/loadLSST.bash)
    - $ eups distrib install sims_photUtils -t sims 
  - Change MINWAVELEN, MAXWAVELEN, WAVELENSTEP in Sed.py, Bandpass.py (part of sims_photUtils) to match MODTRAN data: 300, 1100.5, 0.5 respectively. 
  - Run the following as command line or include in the loadLSST.bash file (or equivalent):
    - setup sims_photUtils -t sims
    - setup sims_sed_library -t sims
      
Coming Soon:
 - AtmoBuilderUsage.ipynb
