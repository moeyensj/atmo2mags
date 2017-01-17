
## AtmoBuilder / Atmo Tutorial

Sections:
- [Prerequisites and Installation](#prerequisites-and-installation)
- [Background Science and Equations](#background-science-and-equations)
- [The Basics: Initialize AtmoBuilder](#the-basics-initialize-atmobuilder)
- [The Basics: Atmo Object](#the-basics-atmo-class)
- [The Basics: SED Types](#the-basics-sed-types)
- [The Basics: Simple Operations](#the-basics-simple-operations)
- [Regressions: Two Parameters](#regressions-two-parameters)

### Prerequisites and Installation

If you are viewing this notebook through notebookviewer online and wish to install and use AtmoBuilder then please clone the [github repository](#https://github.com/moeyensj/atmo2mags) to get started.

AtmoBuilder requires the following to installed and appropriately sourced (including all pre-requisite software):
    - LSST Software Stack
        - sims_photUtils
        - sims_sed_library
    - AstroML

AtmoBuilder's wavelength gridding is limited by the availabe MODTRAN data. As a result this requires the PhysicalParameters.py in sims_photUtils to be altered as follows:

Stack Default:  
    - minwavelen = 300
    - maxwavelen = 1150.1
    - wavelenstep = 0.1
    
Required by AtmoBuilder:
    - minwavelen = 300
    - maxwavelen = 1100
    - wavelenstep = 0.5
   
In the future this step may not be necessary, keep an eye out on updates to the github repository!

### Background Science and Equations

Atmospheric extinction can be modeled by 5 individual components: molecular (Rayleigh) scattering, 
aerosol (Mie) scattering, molecular absorption by ozone (O$_3$), water vapor and combined O$_2$ trace species (see Section 5.2 in LSST document LSE-180):

$S^{atm}(\lambda) = \Pi_{k=1}^5 \, {\rm e}^{- t_k \, \tau_k^{std}(\lambda, X)} 
                           =  {\rm e}^{-\sum_{k=1}^5  t_k \, \tau_k^{std}(\lambda, X)}$

Here $X$ is airmass, $t_k$ is the ratio of the column density (or optical depth) of the $k$-th 
component and its value for the standard atmosphere, and $\tau_k^{std}(\lambda, X)$ is the 
wavelength-dependent optical depth of the $k$-th component for the standard 
atmosphere (with $t_k=1$) viewed at airmass $X$. 

Aerosol (Mie) scattering data is not provided with MODTRAN data and has been modeled with a power-law relationship:

$ \tau_{aerosol}(\lambda, X, \alpha) = 0.1 \, X \, \left({550 \, {\rm nm} \over \lambda} \right)^\alpha $

The power-law exponent ($\alpha$) is chosen to be the 6th free parameter when modeling atmospheres. In summation, we can define a single parameter array of $t_k$ values for each of the 5 extincton components and 1 for the power-law exponent:

$ [t_{H_2O}, t_{O_2}, t_{O_3}, t_{Rayleigh}, t_{Aerosol}, t_{\alpha}] $

The standard atmosphere is defined by the following parameter array and airmass:

    - P: [1.0,1.0,1.0,1.0,1.0,1.7]
    - X: 1.2
    
In order to create an atmosphere for use with AtmoBuilder, a parameter array and airmass like the one above must be given. Details on the specific procedure can be found in the following sections.

### The Basics: Initialize AtmoBuilder

Importing the AtmoBuilder class:


```python
from atmo2mags import AtmoBuilder

% matplotlib inline
```

Initializing the class:


```python
ab = AtmoBuilder()
```

    Found 16 MODTRAN files:
    Pachon_MODTRAN.10.7sc
    Pachon_MODTRAN.11.7sc
    Pachon_MODTRAN.12.7sc
    Pachon_MODTRAN.13.7sc
    Pachon_MODTRAN.14.7sc
    Pachon_MODTRAN.15.7sc
    Pachon_MODTRAN.16.7sc
    Pachon_MODTRAN.17.7sc
    Pachon_MODTRAN.18.7sc
    Pachon_MODTRAN.19.7sc
    Pachon_MODTRAN.20.7sc
    Pachon_MODTRAN.21.7sc
    Pachon_MODTRAN.22.7sc
    Pachon_MODTRAN.23.7sc
    Pachon_MODTRAN.24.7sc
    Pachon_MODTRAN.25.7sc
    MODTRAN files have been read.
    
    Read filter data from LSST software stack.
    Filters: ['u', 'g', 'r', 'i', 'z', 'y4']
    Read hardware data from LSST software stack.


When initialized AtmoBuilder will output the code above to notify you that it has found and read MODTRAN data, and that it has found filter and hardware data from the LSST software stack. 

Retrieving MODTRAN data:


```python
# List of read airmasses
ab.airmasses
```




    ['1.000',
     '1.100',
     '1.200',
     '1.300',
     '1.400',
     '1.500',
     '1.600',
     '1.700',
     '1.800',
     '1.900',
     '2.000',
     '2.100',
     '2.200',
     '2.300',
     '2.400',
     '2.500']




```python
# Standard transmission profile keyed to airmass and component
ab.transDict[1.2]['Rayleigh']
```




    array([ 0.    ,  0.3524,  0.3551, ...,  0.9949,  0.9949,  0.9949])



Plotting the filter and hardware response curves (hardware includes filters):


```python
ab.filterPlot()
```


![png](README/output_18_0.png)



```python
ab.hardwarePlot()
```


![png](README/output_19_0.png)


### The Basics: Atmo Class

AtmoBuilder utilizes an independent atmo object which keeps track of several useful elements. Let's start by creating two atmospheres:


```python
atmo = ab.buildAtmo([1.0,1.0,1.0,2.2,1.4,1.7],1.4)
atmo_std = ab.buildAtmo([1.0,1.0,1.0,1.0,1.0,1.7],1.2)
```

The line of code above has two important keyword parameters:
    - P: parameter array containing 6 values, one for each of the atmospheric extinction components
    - X: airmass

Retrieving atmo attributes:


```python
# Parameters and airmass
atmo.P, atmo.X
```




    ([1.0, 1.0, 1.0, 2.2, 1.4, 1.7], 1.4)




```python
# Total transmission [array]
atmo.sb
```




    array([  0.00000000e+00,   5.62775604e-04,   7.80424731e-04, ...,
             7.95777857e-01,   7.87900709e-01,   7.91135642e-01])




```python
# Individual component transmission profiles [component-keyed dictionary]
atmo.sbDict
```




    {'Aerosol': array([ 0.54288017,  0.54343224,  0.54398304, ...,  0.84640984,
             0.84647403,  0.84653816]),
     'H2O': array([ 0.    ,  1.    ,  1.    , ...,  0.9525,  0.943 ,  0.9468]),
     'O2': array([ 0.,  1.,  1., ...,  1.,  1.,  1.]),
     'O3': array([ 0.    ,  0.015 ,  0.0204, ...,  1.    ,  1.    ,  1.    ]),
     'Rayleigh': array([ 0.        ,  0.06903965,  0.07032593, ...,  0.98706593,
             0.98706593,  0.98706593])}



AtmoBuilder has a lot of various plotting functions including one for atmospheres. The following function will take an atmo object and plot its transmission profile:


```python
ab.transPlot(atmo)
```


![png](README/output_29_0.png)


AtmoBuilder can also plot the individual components:


```python
ab.transPlot(atmo, includeComponents=True)
```


![png](README/output_31_0.png)


Want to compare the atmosphere to another atmosphere?


```python
ab.transPlot(atmo, atmo2=atmo_std)
```


![png](README/output_33_0.png)


Of course, keeping track of a standard atmosphere like the one we created can become laborious. To avoid this, internally in AtmoBuilder there are class constants which define the parameters and airmass for a standard atmosphere. For all major functionality you need not worry about keeping track of the the standard atmosphere. For example:


```python
ab.transPlot(atmo,includeStdAtmo=True)
```


![png](README/output_35_0.png)


### The Basics: SED Types

In order to run regressions and to use many of the basic operations and functions in AtmoBuilder we need SED data to be read in. This is not done when the class is initialized because it can take a fair amount of time depending on how many SEDs you want to read in.

By default AtmoBuilder uses six SED types:
    - Kurucz Main Sequence Stars (read from LSST's SED library) ['mss']
    - White Dwarf Stars (read from LSST's SED library) ['wds']
    - Galaxies (generated at different redshifts from LSST's SED library) ['gals']
    - MLT Dwarf Stars (read from LSST's SED library) ['mlts']
    - Quasars (generated at different redshifts from AtmoBuilder repository) ['qsos']
    - Supernovae (generated at different redshifts from AtmoBuilder repository) ['sns']

Reading in Kurucz MS stars:


```python
ab.readMSs()
```

    # Read 988 MS stars from /Users/joachim/lsst/DarwinX86/sims_sed_library/2014.10.06/starSED/kurucz/


You can call similar functions for each of the six SED types to read them. If you want to use all of the available SED types (or a large subset) you can call the following (this will take a little bit of time):


```python
ab.readAll()
```

    # Read 988 MS stars from /Users/joachim/lsst/DarwinX86/sims_sed_library/2014.10.06/starSED/kurucz/
    # Read 849 white dwarfs from /Users/joachim/lsst/DarwinX86/sims_sed_library/2014.10.06/starSED/wDs/
    # Read 74 mlt stars from /Users/joachim/lsst/DarwinX86/sims_sed_library/2014.10.06/starSED/mlt/
    # Generated 2520 galaxies at redshifts between 0.000000 and 3.000000
    # Generated 76 quasars at redshifts between 0.000000 and 7.500000
    # Generated 39 sn's at redshifts between 0.000000 and 1.200000 on days ['0', '20', '40']


If you want to run some simple operations it may be helpful to retrieve the SED data you read in. With AtmoBuilder it is remarkably easy to do so as it is all stored as class attributes. Lets retrieve some quasar data:


```python
qso = ab.qsos
qso_z = ab.qsoRedshifts
```

We can see what the SED attributes are set when SED data is read in two ways:


```python
ab?
```

... or we can look at each of the individual SED read functions:


```python
ab.readMSs?
```

### The Basics: Simple Operations

Once you have some SED data read in you can use them to calculate magnitudes, delta magnitudes, and so on... 

Lets combine hardware data and an atmosphere object to get the combined throughput:


```python
throughput = ab.combineThroughputs(atmo)
throughput_std = ab.combineThroughputs(atmo_std)
```

We can then also plot the combined throughput:


```python
ab.throughputPlot(throughput)
```


![png](README/output_53_0.png)


Like transPlot, we can plot two throughputs or include the standard throughput:


```python
ab.throughputPlot(throughput, bpDict2=throughput_std)
ab.throughputPlot(throughput, includeStdAtmo=True)
```


![png](README/output_55_0.png)



![png](README/output_55_1.png)


A characteristic function we care about is the normalized bandpass response function for a given throughput. The phi function can be retrieved from a filter-keyed bandpass dictionary as follows (example for u filter):


```python
throughput['u'].phi
```




    array([  0.00000000e+00,   2.41980424e-10,   6.86163658e-10, ...,
             4.21530335e-10,   3.53349898e-10,   5.46307213e-19])



Because these throughputs are bandpass objects, you can also access the wavelength range and the overall throughput:


```python
throughput['u'].sb
```




    array([  0.00000000e+00,   1.06043522e-09,   3.01199080e-09, ...,
             6.75594271e-09,   5.66577908e-09,   8.76373126e-18])




```python
throughput['u'].wavelen
```




    array([  300. ,   300.5,   301. , ...,  1099. ,  1099.5,  1100. ])



AtmoBuilder can also plot phi functions:


```python
ab.phiPlot(throughput)
```


![png](README/output_62_0.png)


Of course it might be handy to compare one phi function to another:


```python
ab.phiPlot(throughput,throughput_std)
```


![png](README/output_64_0.png)


AtmoBuilder can also plot the change in normalized bandpass response function:


```python
ab.dphiPlot(throughput, throughput_std)
```


![png](README/output_66_0.png)


We might want to compare one dphi to another:


```python
throughput_temp = ab.combineThroughputs(ab.buildAtmo([2.3,2.1,1.2,0.9,1.3,1.8],1.4))
ab.dphiPlot(throughput, throughput_std, throughput_temp)
```


![png](README/output_68_0.png)


Sometimes its hard to really see the difference between two phi functions as is often the case when plotting regression data. We can use the ddphiPlot to plot the difference between to dphis in the form of a subtraction between them:


```python
ab.ddphiPlot(throughput,throughput_temp,throughput_std)
```


![png](README/output_70_0.png)


Using the throughput we can calculate the magnitudes of stars we read in (defaults to using Kurucz MS SEDs)


```python
mags = ab.mags(throughput)
mags_std = ab.mags(throughput_std)
```

Or we can calculate the magnitudes and change in magnitudes of other SEDs:


```python
qso_mags = ab.mags(throughput, seds=qso, sedkeylist=qso_z)
qso_mags_std = ab.mags(throughput_std, seds=qso, sedkeylist=qso_z)
```

    AtmoBuilder.py:884: FutureWarning: comparison to `None` will result in an elementwise object comparison in the future.
      if sedkeylist == None:



```python
dmags = ab.dmags(mags, mags_std)
qso_dmags = ab.dmags(qso_mags, qso_mags_std)
```

Both mags and dmags are returned as filter-keyed dictionaries of magnitudes:


```python
mags['u']
dmags['u']
```




    array([ -60.23624135,  -60.40850702,  -60.56826619,  -60.74566544,
            -60.91177112,  -61.08369099,  -61.25698358,  -61.99937469,
            -62.66635849,  -63.32696276,  -64.00057694,  -64.66445761,
            -65.33650774,  -66.00255451,  -66.66112885,  -67.33962434,
            -68.01727869,  -68.68508883,  -69.3486042 ,  -69.68886445,
            -69.91999707,  -70.15160589,  -70.3832493 ,  -70.60231546,
            -70.83213889,  -71.05492605,  -71.28677986,  -71.51070901,
            -71.73450645,  -71.95922916,  -72.18802748,  -72.41771595,
            -72.53058869,  -72.72820371,  -72.93014927,  -73.132368  ,
            -73.33655565,  -73.53354453,  -73.73531967,  -73.94081379,
            -74.13836494,  -74.34397373,  -74.53915895,  -74.73827928,
            -74.9433136 ,  -75.04396874,  -59.52402349,  -59.69855054,
            -59.87310733,  -60.0551013 ,  -60.22616858,  -64.55743204,
            -65.23532682,  -65.9216536 ,  -66.59582775,  -67.27782341,
            -67.96125851,  -68.64347645,  -69.32404503,  -69.66348194,
            -69.90049712,  -70.14033457,  -70.36996138,  -70.60783434,
            -70.84684351,  -71.07597787,  -71.32100036,  -71.55846612,
            -71.78867348,  -72.02556848,  -72.25842021,  -72.49352107,
            -72.61787325,  -72.82266811,  -73.02214165,  -73.2315028 ,
            -81.27689532,  -81.62754728,  -81.97118041,  -82.31789381,
            -76.01767967,  -75.39459305,  -74.79016611,  -74.19762842,
            -73.61790584,  -73.05243811,  -72.49593438,  -77.06009157,
            -77.34977956,  -77.63712864,  -77.92769151,  -53.13487875,
            -53.17204579,  -53.20372047,  -53.23868547,  -74.52555866,
            -74.44455135,  -74.36559598,  -60.61545957,  -74.3316647 ,
            -74.32737852,  -74.32483912,  -74.32443684,  -74.32503356,
            -74.33166749,  -74.35387338,  -74.3688597 ,  -74.39606558,
            -74.4169389 ,  -74.44510294,  -74.48307664,  -74.51990728,
            -58.16282678,  -57.91452287,  -57.68063794,  -57.45445311,
            -57.21993644,  -57.00045331,  -56.7810099 ,  -56.57767205,
            -56.35062053,  -56.15995865,  -55.95081383,  -55.75049548,
            -55.56084755,  -74.54147212,  -74.62642796,  -74.70444065,
            -74.7865556 ,  -74.87225916,  -74.96549384,  -75.05287665,
            -75.14097372,  -75.23394256,  -75.31929274,  -75.41676617,
            -75.50843317,  -75.60408753,  -55.4626833 ,  -55.2990263 ,
            -55.14843755,  -54.99466203,  -54.84151896,  -54.69055569,
            -54.54452003,  -54.40485222,  -54.2622259 ,  -54.12204818,
            -53.98483925,  -53.71759728,  -75.65001982,  -75.8503057 ,
            -76.0455661 ,  -76.24235527,  -76.44127804,  -76.63895707,
            -76.83761866,  -77.03730666,  -77.23602519,  -77.43618269,
            -77.63486273,  -77.83232918,  -78.03556837,  -53.64882934,
            -53.51646906,  -53.38853758,  -53.25523235,  -53.12632906,
            -53.00275709,  -52.87204056,  -52.75043919,  -52.62236111,
            -52.49420722,  -78.13389795,  -78.44596891,  -78.75953864,
            -79.07167419,  -79.38301796,  -79.69936356,  -80.00943957,
            -80.31960105,  -80.63038525,  -80.94461882,  -81.25315381,
            -81.56445003,  -81.88067133,  -51.48921194,  -51.4156024 ,
            -51.34699472,  -82.03337031,  -82.3956664 ,  -82.74751966,
            -83.11960392,  -83.4760254 ,  -76.64269017,  -76.43570475,
            -76.23346998,  -76.04249687,  -75.8254595 ,  -75.63051196,
            -75.43058508,  -75.22554744,  -75.03250423,  -74.82962963,
            -74.63609879,  -74.15640921,  -74.1500809 ,  -74.15745049,
            -74.15213616,  -74.16100426,  -74.15392265,  -74.15901407,
            -74.15574484,  -74.16134754,  -74.16677407,  -74.16903447,
            -74.16777365,  -74.17282006,  -54.26554062,  -53.97204235,
            -53.67362206,  -53.38425277,  -52.80573001,  -74.1729936 ,
            -74.378248  ,  -74.57585198,  -74.77799222,  -74.97951527,
            -75.18448016,  -75.38693033,  -75.5853645 ,  -75.78650837,
            -75.98913024,  -76.18580718,  -76.39013381,  -76.59072937,
            -52.66275721,  -52.44724008,  -52.24001703,  -52.02110629,
            -51.81826069,  -51.60413008,  -51.39482373,  -51.18693049,
            -50.97949997,  -50.77277202,  -50.56368064,  -50.3590797 ,
            -50.1522377 ,  -76.69096985,  -77.03416777,  -77.3766119 ,
            -77.71767727,  -78.06402224,  -78.40386794,  -78.74529612,
            -79.08706858,  -79.42860879,  -79.77297255,  -80.10596699,
            -80.4438102 ,  -80.78460329,  -50.05312819,  -49.9736761 ,
            -49.89522255,  -49.81931451,  -49.74024005,  -49.66420725,
            -49.58869353,  -49.51797908,  -49.44146371,  -49.36437778,
            -49.29258898,  -49.21790087,  -49.14344944,  -80.95289563,
            -81.33350045,  -81.71989   ,  -82.10091639,  -82.49489692,
            -82.875586  ,  -83.255307  ,  -49.11316924,  -49.14567356,
            -49.19047843,  -49.23730199,  -49.27729172,  -49.31620662,
            -49.35664158,  -49.39249513,  -49.45027481,  -49.48395335,
            -49.52011389,  -49.56812697,  -73.91716046,  -74.14587426,
            -74.37676676,  -74.60787844,  -74.83258527,  -75.06090161,
            -75.29023098,  -49.94183416,  -49.74481108,  -49.54772904,
            -49.3560214 ,  -49.1599496 ,  -48.96577795,  -48.77231837,
            -48.5806793 ,  -48.38479851,  -75.40659146,  -75.76467891,
            -76.12691081,  -76.48583002,  -76.8416848 ,  -77.20192808,
            -77.55720095,  -77.91652991,  -78.27368369,  -78.62264732,
            -78.99328377,  -79.33325018,  -79.69626575,  -48.29330832,
            -48.23973737,  -48.1884778 ,  -48.13905835,  -48.09008041,
            -48.0413517 ,  -47.99203793,  -47.94394768,  -47.896292  ,
            -47.84964398,  -47.801565  ,  -47.75500522,  -47.70889452,
            -79.88064301,  -80.27521475,  -80.67229624,  -81.08532529,
            -81.48158099,  -81.88294344,  -82.29249778,  -82.68589332,
            -47.6765522 ,  -47.74127043,  -47.79915311,  -47.86402462,
            -47.91470523,  -47.97477663,  -48.02661995,  -48.08196115,
            -48.13537526,  -48.19268223,  -48.24830803,  -48.32173465,
            -94.61236108,  -95.97552566,  -48.95535512,  -49.07186926,
            -49.219737  ,  -49.36751334,  -49.49811506,  -49.64848341,
            -49.77024546,  -49.91501155,  -50.0637811 ,  -50.12878758,
            -50.29936941,  -50.49163195,  -50.68223429,  -50.85628441,
            -71.28272945,  -71.30989179,  -71.55544158,  -71.80454168,
            -72.05557059,  -72.305308  ,  -72.55436808,  -72.80061197,
            -73.05116048,  -73.29858276,  -73.54632854,  -73.79565656,
            -74.04125747,  -74.29007229,  -74.41253938,  -74.7804466 ,
            -75.14883982,  -75.51346074,  -75.87788381,  -76.24440045,
            -76.60627633,  -76.96566361,  -77.34006217,  -77.69323531,
            -78.05961786,  -78.42312211,  -78.77363576,  -46.90649116,
            -46.87957869,  -46.8539378 ,  -46.82723429,  -46.80231829,
            -46.77448331,  -78.95871421,  -79.36352024,  -79.79365441,
            -80.20199797,  -80.6205534 ,  -81.021797  ,  -81.43633807,
            -81.86203814,  -82.25915057,  -46.76542244,  -46.81536593,
            -46.88001155,  -46.94871999,  -47.02169495,  -47.06966611,
            -47.12369617,  -47.20888888,  -47.27800027,  -47.33780205,
            -47.39915977,  -47.45714633,  -47.52185678,  -96.7272796 ,
            -98.11162544,  -99.4809028 , -100.86011999,  -47.55823051,
            -47.69091453,  -47.84353897,  -47.97209662,  -48.12800943,
            -48.27526702,  -48.40286147,  -48.54546014,  -48.69646312,
            -48.82314174,  -48.97385388,  -49.1111279 ,  -49.25660228,
            -49.33621703,  -49.53078762,  -49.7030987 ,  -49.90157738,
            -50.09830605,  -50.29226342,  -69.62490838,  -69.68838243,
            -69.7628741 ,  -69.82435337,  -69.89386727,  -69.95641348,
            -70.02993371,  -70.09553254,  -70.16462503,  -70.23220915,
            -70.30250557,  -70.34033848,  -70.596495  ,  -70.85212951,
            -71.10926683,  -71.36123578,  -71.61946058,  -71.87363352,
            -72.13455658,  -72.38667427,  -72.63780861,  -72.89641277,
            -73.14669998,  -73.40469246,  -73.53007707,  -73.904201  ,
            -74.27294589,  -74.64388462,  -75.0177044 ,  -75.38857563,
            -75.75817567,  -76.12606884,  -76.4966294 ,  -76.86379591,
            -77.23031805,  -77.58939197,  -77.96675591,  -46.15759892,
            -46.13622108,  -46.11706323,  -46.09591372,  -46.07798317,
            -78.15990188,  -78.59069898,  -79.02642436,  -79.4661429 ,
            -79.89963929,  -80.33073715,  -80.76760655,  -81.1979908 ,
            -81.63141691,  -46.06383624,  -46.13773356,  -46.19740472,
            -46.2705788 ,  -46.33959804,  -46.40563024,  -46.46652156,
            -46.55212721,  -46.62090774,  -46.68786244,  -46.74382056,
            -46.8144954 ,  -46.89165134,  -46.92917658,  -47.06696677,
            -47.22834551,  -47.38136189,  -47.51424422,  -47.6651414 ,
            -47.81799677,  -47.95487684,  -48.11872478,  -48.26344531,
            -48.40898249,  -48.56138103,  -48.70004618,  -99.98683805,
           -101.4460491 ,  -48.77859978,  -48.97609087,  -49.17694905,
            -49.3767663 ,  -49.56825322,  -49.78079311,  -68.92943602,
            -69.00212525,  -69.07877384,  -69.14908472,  -69.22716122,
            -69.30237598,  -69.37726837,  -69.44669906,  -69.52439285,
            -69.56329759,  -69.81383937,  -70.06983576,  -70.31950362,
            -70.57317728,  -70.82480081,  -71.07510068,  -71.32798194,
            -71.57483943,  -71.82591559,  -72.07557205,  -72.32439819,
            -72.57786817,  -72.6982858 ,  -73.09219346,  -73.48874615,
            -73.88297617,  -74.2741372 ,  -74.66855895,  -75.06072269,
            -75.45242935,  -75.8469124 ,  -76.22535387,  -76.61905114,
            -77.00764092,  -77.39879926,  -45.69510196,  -45.68415201,
            -45.66272484,  -77.59356899,  -78.02679475,  -78.47353394,
            -78.90448184,  -79.33830329,  -79.77477653,  -80.20900137,
            -80.64476739,  -81.08413845,  -81.51101035,  -45.65794563,
            -45.73431549,  -45.81222844,  -45.88976122,  -45.95860974,
            -46.03862282,  -46.11320229,  -46.18694297,  -46.25827172,
            -46.33683983,  -46.41731943,  -46.49337685,  -46.56956824,
            -46.61210264,  -46.76088765,  -46.92121721,  -47.04787288,
            -47.22431698,  -47.36670126,  -47.52240068,  -47.66897465,
            -47.8236443 ,  -47.98354779,  -48.12804979,  -48.29744772,
            -48.43609552,  -99.79664333, -101.06750075, -102.33222161,
           -103.58747483,  -48.51750619,  -48.71960351,  -48.91566976,
            -49.11157813,  -49.32747985,  -49.53101082,  -68.41014301,
            -68.4990069 ,  -68.58116383,  -68.66454597,  -68.7505224 ,
            -68.8341587 ,  -68.91805281,  -69.00131356,  -69.04502461,
            -69.3163815 ,  -69.58677886,  -69.85455201,  -70.12858571,
            -70.39658177,  -70.66759735,  -70.93452155,  -71.20135835,
            -71.47236795,  -71.73711429,  -72.00577401,  -72.26739217,
            -72.4059273 ,  -72.79751405,  -73.18604892,  -73.57567576,
            -73.96561916,  -74.35521945,  -74.74775173,  -75.12583949,
            -75.50927883,  -75.90532877,  -76.2880743 ,  -76.67100833,
            -77.04646767,  -77.24562136,  -77.70375988,  -78.15053581,
            -78.59708452,  -79.05373178,  -79.49493685,  -79.95217412,
            -80.39579995,  -80.84563792,  -81.29723321,  -45.6169505 ,
            -45.70252952,  -45.78160795,  -45.87384054,  -45.96057154,
            -46.04219094,  -46.13010183,  -46.21548754,  -46.29691175,
            -46.37631936,  -46.46528451,  -46.49960335,  -46.6671565 ,
            -46.82818256,  -46.959673  ,  -47.11565469,  -47.2808072 ,
            -47.43688866,  -47.58420357,  -47.74137153,  -47.90745094,
            -48.05486281,  -48.20261459,  -48.36692556, -102.53220335,
           -104.00580061,  -48.44166174,  -48.6400586 ,  -48.85001235,
            -49.06724558,  -49.27675139,  -49.46812866,  -49.68367022,
           -129.41061696,  -68.195978  ,  -68.29016758,  -68.38118748,
            -68.47388568,  -68.56285886,  -68.65439673,  -68.74806107,
            -68.79413805,  -69.07040386,  -69.34517232,  -69.62150829,
            -69.89482857,  -70.16962406,  -70.4438721 ,  -70.71404078,
            -70.98665941,  -71.26170807,  -71.53279695,  -71.80726856,
            -72.07544159,  -72.21221651,  -72.60724671,  -73.00073077,
            -73.39392013,  -73.78762312,  -74.17968153,  -74.57159471,
            -74.95384399,  -75.35491996,  -75.73831841,  -76.12884345,
            -76.51340536,  -76.90014724,  -77.09663142,  -77.53584375,
            -78.00000924,  -78.4452303 ,  -78.90364584,  -79.34496273,
            -79.79469108,  -80.24096199,  -80.68330963,  -81.12183998,
            -45.87014559,  -45.94578364,  -46.02459358,  -46.10932411,
            -46.18029746,  -46.25490098,  -46.34505423,  -46.4275746 ,
            -46.49696488,  -46.54073965,  -46.69141327,  -46.85891635,
            -46.99617613,  -47.16293236,  -47.31592339,  -47.48119207,
            -47.61706328,  -47.7819168 ,  -47.9264062 ,  -48.09022897,
            -48.2342557 ,  -48.3914402 , -103.87315895,  -48.46416521,
            -48.6702458 ,  -49.08536775,  -49.27446553,  -49.68564785,
           -129.79755795,  -68.1754787 ,  -68.26596138,  -68.36250164,
            -68.45940016,  -68.55533645,  -68.65079636,  -68.69903653,
            -68.97446184,  -69.2519229 ,  -69.53033549,  -69.807889  ,
            -70.085873  ,  -70.36089194,  -70.63763876,  -70.91171555,
            -71.18633597,  -71.46085686,  -71.73116001,  -72.00703944,
            -72.14351923,  -72.54479316,  -72.94695039,  -73.34713041,
            -73.74868486,  -74.1452312 ,  -74.54730469,  -74.93647818,
            -75.33559855,  -75.74136033,  -76.13243224,  -76.51450814,
            -76.92401687,  -77.11331811,  -77.57134053,  -78.0259911 ,
            -78.48711573,  -78.94169344,  -79.39814683,  -79.84481853,
            -80.29931078,  -80.75678578,  -81.20618763,  -45.97342402,
            -46.04511335,  -46.12425942,  -46.21487364,  -46.28846347,
            -46.36853761,  -46.46027481,  -46.54299678,  -46.58667202,
            -46.73476493,  -46.87632102,  -47.0342778 ,  -47.19724339,
            -47.36084606,  -47.49549633,  -47.65429149,  -47.81156536,
            -47.96377786,  -48.11695048,  -48.25787757,  -48.41027367,
           -103.94097576, -105.40200598, -106.86094067, -108.31144278,
            -48.49461597,  -48.69270438,  -48.8928376 ,  -49.0950278 ,
            -49.30695855,  -49.51055829,  -49.70397607, -109.03462707,
           -109.67799879, -130.00191786, -129.63563436, -129.44810534,
           -129.27029705, -129.08643713, -128.90523984,  -61.36021931,
            -61.51988041,  -61.60521744,  -62.2384942 ,  -62.87225518,
            -63.50888666,  -64.14878937,  -64.79049044,  -65.44131841,
            -66.07611353,  -66.71699826,  -67.35776927,  -67.99896874,
            -68.6518172 ,  -69.2952062 ,  -69.61474588,  -69.82537088,
            -70.04131052,  -70.2477687 ,  -70.46049099,  -70.68116784,
            -70.89089848,  -71.09763131,  -71.31353123,  -71.51928242,
            -71.73552853,  -71.94893336,  -72.16215246,  -72.2621171 ,
            -72.45134071,  -72.64032937,  -72.83916589,  -73.02590817,
            -73.21706587,  -73.40337755,  -73.59773279,  -73.78738454,
            -73.97213881,  -74.16417223,  -74.3624843 ,  -74.5527421 ,
            -74.65236621,  -74.79150736,  -74.94040124,  -75.09432419,
            -75.23740113,  -75.38957005,  -72.24743646,  -72.43518364,
            -72.61211043,  -72.80195091,  -72.99076052,  -73.17740898,
            -73.37324751,  -73.55428047,  -73.74336879,  -73.92925578,
            -74.11320106,  -74.30246831,  -74.3913525 ,  -74.54372105,
            -74.68685065,  -74.83443658,  -74.98921607,  -75.13049815,
            -75.27848711,  -75.42524826,  -75.57431955,  -73.28947256,
            -73.47001051,  -73.65820331,  -73.83836146,  -74.02536237,
            -74.117932  ,  -74.25303296,  -74.40042473,  -74.54558563,
            -74.68868334,  -74.8278027 ,  -74.9713895 ,  -75.11261732,
            -75.25924289,  -75.39935994,  -75.54117478,  -75.68303489,
            -75.83309809,  -75.89985858,  -73.97876937,  -74.11533093,
            -74.25464507,  -74.39092125,  -74.5323158 ,  -74.66907296,
            -74.81229009,  -74.95221131,  -75.09408213,  -75.16389834,
            -75.25208406,  -75.3443716 ,  -75.43006959,  -75.51997562,
            -75.60592133,  -75.70305266,  -75.79312316,  -75.88508035,
            -77.02070042,  -77.0119938 ,  -76.99029197,  -76.96617623,
            -76.95171587,  -72.80601591,  -72.90458517,  -72.99909662,
            -74.72528163,  -74.74175145,  -74.74969756,  -74.75230399,
            -74.77152343,  -74.77677517,  -74.8116063 ,  -74.8166487 ,
            -74.81812062,  -74.84723709,  -74.84956201,  -74.81141376,
            -74.77552669,  -74.72827258,  -74.68341388,  -74.65758997,
            -74.61992286,  -74.57277144,  -74.54359607,  -74.50598477,
            -74.47106899,  -74.43022037,  -74.39798119,  -74.38531938])



Suppose we want to see the change in magnitudes between two atmospheres for a single SED type, we can use the following function: 


```python
ab.dmagPlot(throughput,throughput_std,'qsos')
```


![png](README/output_79_0.png)


Or another SED type example:


```python
ab.dmagPlot(throughput,throughput_std,'gals')
```


![png](README/output_81_0.png)


### Regressions: Two Parameters

Regressions (currently limited to two varying parameters) are the most integral part of AtmoBuilder and are run using a single function. 


```python
ab.computeAtmoFit?
```

Making use of the atmospheres we created at the beginning of this tutorial, lets run some regressions. Note that we only changed two parameters relative to the standard atmosphere when we created the atmosphere. These were the parameters for the Rayleigh and Aerosol components of atmospheric extinction, so let us regress those two components. (I have added the pickleString keyword argument which saves the pickle dump files and plots with an additional string component in its file name)


```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 1.94 2.14
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43


    /Users/joachim/lsst/DarwinX86/anaconda/2.2.0/lib/python2.7/site-packages/matplotlib/text.py:52: UnicodeWarning: Unicode equal comparison failed to convert both arguments to Unicode - interpreting them as being unequal
      if rotation in ('horizontal', None):
    /Users/joachim/lsst/DarwinX86/anaconda/2.2.0/lib/python2.7/site-packages/matplotlib/text.py:54: UnicodeWarning: Unicode equal comparison failed to convert both arguments to Unicode - interpreting them as being unequal
      elif rotation == 'vertical':



![png](README/output_86_2.png)



![png](README/output_86_3.png)



![png](README/output_86_4.png)


Sometimes the contour plots do not adequately show the log-likelihood matrix, we can directly plot the logL matrix instead of the regression contours using 'useLogL':


```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, useLogL=True, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 1.94 2.14
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_88_1.png)



![png](README/output_88_2.png)



![png](README/output_88_3.png)


It is also possible to plot both the contours and the log-likelihood matrix on the same plot:


```python
ab.computeAtmoFit('Rayleigh','Aerosol',atmo, plotBoth=True, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 1.94 2.14
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_90_1.png)



![png](README/output_90_2.png)



![png](README/output_90_3.png)


A colorbar for the log-likelihood column might be handy and can be included using a keyword argument. Unfortunately, it does resize the middle column slightly. 


```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, useLogL=True, includeColorBar=True, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_mss_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 1.94 2.14
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_92_1.png)



![png](README/output_92_2.png)



![png](README/output_92_3.png)


By default regressions will run using the Kurucz main sequence star SED data. However, we can also specify a different regression SED: (I have also set generate Dphi to False to condense the notebook a little)


```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, regressionSed='qsos', generateDphi=False, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 2.04 1.84
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_94_1.png)


Like the regression SEDs, we can also set which comparison SEDs to use. By default, computeAtmoFit will always use the remaining SEDs (all of them except the original regression SED).


```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, regressionSed='qsos', comparisonSeds=['gals','sns'], generateDphi=False, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E5_qsos_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 2.04 1.84
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_96_1.png)


The default percent error value for dmags in mmags is set to be 0.5%, however we can also set it to other values. For example, we can set it to be 5%:


```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, err=0.05, generateDphi=False, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E50_mss_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E50_mss_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E50_mss_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E50_mss_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E50_mss_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG0_E50_mss_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 1.94 2.14
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_98_1.png)


One additional feature of interest is the deltaGrey keyword parameter. It is essentially an additional extinction factor implemented to add some offset representative of cloudy and potentially grey skies.

It works in the following way:
    - deltaGrey = 0.0 : default, no effect
    - deltaGrey >= 0.0 : subtract deltaGrey value from dmags (dG should be in mmags, for example: 20)
    - deltaGrey < 0.0 : subtract median dmag value from dmags


```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, deltaGrey=20, generateDphi=False, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG200_E5_mss_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG200_E5_mss_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG200_E5_mss_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG200_E5_mss_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG200_E5_mss_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG200_E5_mss_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 1.94 2.14
    r 2.35 1.22
    i 2.14 1.43
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_100_1.png)



```python
ab.computeAtmoFit('Rayleigh', 'Aerosol', atmo, deltaGrey=-1, generateDphi=False, pickleString='Tutorial')
```

    Computing nonlinear regression for Rayleigh and Aerosol.
    Observed atmosphere parameters: [1.0, 1.0, 1.0, 2.2, 1.4, 1.7]
    Observed atmosphere airmass:    1.4
    Standard atmosphere parameters: [1.0, 1.0, 1.0, 1.0, 1.0, 1.7]
    Standard atmosphere airmass:    1.2
    Observed atmosphere parameter for Rayleigh: 2.2
    Observed atmosphere parameter for Aerosol: 1.4
    
    Calculating best parameters for u filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG-10_E5_mss_u_50b_Tutorial.pkl'
    Saved LogL for u filter.
    Calculating best parameters for g filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG-10_E5_mss_g_50b_Tutorial.pkl'
    Saved LogL for g filter.
    Calculating best parameters for r filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG-10_E5_mss_r_50b_Tutorial.pkl'
    Saved LogL for r filter.
    Calculating best parameters for i filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG-10_E5_mss_i_50b_Tutorial.pkl'
    Saved LogL for i filter.
    Calculating best parameters for z filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG-10_E5_mss_z_50b_Tutorial.pkl'
    Saved LogL for z filter.
    Calculating best parameters for y4 filter...
    @pickle_results: using precomputed results from 'pickles/X14_P101010221417_Rayleigh_Aerosol_XSTD12_DG-10_E5_mss_y4_50b_Tutorial.pkl'
    Saved LogL for y4 filter.
    
    Best fit parameters (Filter, Rayleigh, Aerosol):
    u 2.24 1.12
    g 2.04 1.84
    r 2.35 1.22
    i 3.06 0.82
    z 2.14 1.43
    y4 2.14 1.43



![png](README/output_101_1.png)



```python

```
