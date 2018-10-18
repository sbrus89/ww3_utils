import collections

switch = collections.OrderedDict()

######################
# Mandatory switches
######################

# Machine-dependent code
switch['F90']  = [True,  'Fortran 0 style data and time capturing and program abort']
switch['DUM']  = [False, 'Dummy to be used if WaveWatch II is to be installed on previously untried hardware']

# Hardward model/message passing protocol
switch['SHRD'] = [False, 'Shared memory model (no message passing)']
switch['DIST'] = [True,  'Disrbuted memory model'] 
switch['MPI']  = [True,  'Message Passing Interface']

# Word length to determine record length in direct access files
switch['LRB4'] = [True,  '4 byte words']
switch['LRB8'] = [False, '8 byte words']

# Compliation as subroutine or stand-alone program
switch['NOPA'] = [True,  'Compliation as stand-alone program']
switch['PALM'] = [False, 'Compliation as subroutine']

# Propagation schemes and GSE alleviation method
switch['PR0']  = [False, 'No propagation scheme / GSE alleviation used']
switch['PR1']  = [False, 'First order propagation scheme, no GSE alleviation']
switch['PR2']  = [False, 'Higher-order schemes with Booij and Holthuijsen (1987) dispersion correction']
switch['PR3']  = [True,  'Higher-order schemes with Tolman (2002a) averaging technique']
switch['PRX']  = [False, 'Experimental (user supplied)']
switch['UNO']  = [False, 'Second-order (UNO) propagation scheme']
switch['UQ']   = [True,  'Third-order (UQ) propagation scheme']

# Flux computation
switch['FLX0'] = [True,  'No routine used; flux computation included in source terms']
switch['FLX1'] = [False, 'Friction velocity according to Charnock, 1955; Wu, 1982']
switch['FLX2'] = [False, 'Friction velocity from Tolman and Chalikov input']
switch['FLX3'] = [False, 'Friction velocity from Tolman and Chalikov input with cap']
switch['FLX4'] = [False, 'Friction velocity according to Hwang (2001)']
switch['FLXX'] = [False, 'Experimental (user supplied)']

# Linear input
switch['LN0']  = [False, 'No linear input']
switch['SEED'] = [False, 'Spectral seeding']
switch['LN1']  = [True,  'Cavaleri and Malanotte-Rizzoli with filter']
switch['LNX']  = [False, 'Experimental (user supplied)']

# Input and dissipattion
switch['ST0']  = [False, 'No input and dissipation used']
switch['ST1']  = [False, 'WAM3 source term package']
switch['ST2']  = [False, 'Tolman and Chalikov (1996) source term package']
switch['STAB2']= [False, 'Enable stability correction for ST2']
switch['ST3']  = [False, 'WAM4 and variants source term package']
switch['STAB3']= [False, 'Enable stability correction from Abdalla and Bidlot (2002)']
switch['ST4']  = [True,  'Ardhuin et al. (2010) source term package']
switch['ST6']  = [False, 'BYDRZ source term package']
switch['STX']  = [False, 'Experimental (user supplied)']
switch['STAB0']= [True,  '']

# Nonlinear interactions
switch['NL0']  = [False, 'No nonlinear interactions used']
switch['NL1']  = [True,  'Discrete interaction approximation (DIA)']
switch['NL2']  = [False, 'Exact interaction approximation (WRT)']
switch['NL3']  = [False, 'Generalized Multiple DIA (GMD)']
switch['NL4']  = [False, 'Two-scale approximation (TSA)']
switch['NLX']  = [False, 'Experimental (user supplied)']

# Bottom friction
switch['BT0']  = [False, 'No bottom friction used']
switch['BT1']  = [True,  'JONSWAP bottom friction formulation']
switch['BT4']  = [False, 'SHOWEX bottom friction formulation']
switch['BT8']  = [False, 'Dalrymple and Liu formulation (fluid mud seafloor)']
switch['BT9']  = [False, 'Ng formulation (fluid mud seafloor)']
switch['BTX']  = [False, 'Experimental (user supplied)']

# Damping by sea ice
switch['IC0']  = [True,  'No damping by sea ice']
switch['IC1']  = [False, 'Simple formulation']
switch['IC2']  = [False, 'Liu et al. formulation']
switch['IC3']  = [False, 'Wang and Shen formulation']
switch['IC4']  = [False, 'Frequency-dependent damping by sea ice']

# Scattering by sea ice
switch['IS0']  = [True,  'No scattering by sea ice'] 
switch['IS1']  = [False, 'Diffusive scattering by sea ice (simple)']
switch['IS2']  = [False, 'Floe-size dependent scattering and dissipation']

# Reflection
switch['REF0'] = [True,  'No reflection']
switch['REF1'] = [False, 'Enables reflection of shorelines and icebergs']

# Depth-induced breaking 
switch['DB0']  = [False, 'No depth-induced breaking used']
switch['DB1']  = [True,  'Battjes-Janssen']
switch['DB2']  = [False, 'Experimental (user supplied)']

# Triad interactions
switch['TR0']  = [True,  'No triad interactions used']
switch['TR1']  = [False, 'Lumped Triad Interaction (LTA) method']
switch['TRX']  = [False, 'Experimental (user supplied)']

# Bottom scattering
switch['BS0']  = [True,  'No bottom scattering used']
switch['BS1']  = [False, 'Magne and Ardhuin']
switch['BSX']  = [False, 'Experimental (user supplied)']

# Supplemental source term
switch['XX0']  = [True,  'No supplemental source term used']
switch['XXX']  = [False, 'Experimental (user supplied)']

# Wind interpolation in time
switch['WNT0'] = [False, 'No interpolation']
switch['WNT1'] = [True,  'Linear interpolation']
switch['WNT2'] = [False, 'Approximately quadratic interpolation']

# Wind interpolation in space
switch['WNX0'] = [False, 'Vector interpolation']
switch['WNX1'] = [True,  'Approximately linear speed interpolation']
switch['WNX2'] = [False, 'Approximately quadratic speed interpolation']

# Current interpolation in time
switch['CRT0'] = [False, 'No interpolation']
switch['CRT1'] = [True,  'Linear interpolation']
switch['CRT2'] = [False, 'Approximately quadratic interpolation']

# Current interpolation in space
switch['CRX0'] = [False, 'Vector interpolation']
switch['CRX1'] = [True,  'Approximate linear speed interpolation']
switch['CRX2'] = [False, 'Approximate quadratic speed interpolation']

# User supplied GRIB package
switch['NOGRB']= [True,  'No package included']
switch['NCEP1']= [False, 'NCEP GRIB1 package for IBM SP']
switch['NCEP2']= [False, 'NCEP GRIB2 package for IBM SP']

####################
# Optional switches
####################

# Output control
switch['O0']   = [True,  'Output of namelists in grid preprocessor']
switch['O1']   = [True,  'Output of boundary points in grid preprocessor']
switch['O2']   = [True,  'Output of the grid point status map in grid preprocessor']
switch['O2a']  = [False, 'Generation of land-sea mask file mask.ww3 in grid preprocessor']
switch['O2b']  = [False, 'Output of obstruction map in grid preprocessor']
switch['O2c']  = [False, 'Print status map in grid preprocessor']
switch['O3']   = [True,  'Additional output in loop over fields in field preprocessor']
switch['O4']   = [True,  'Print plot of normalized one-dimensional energy spectrum in initial conditions program']
switch['O5']   = [True,  'Print plot of normalized two-dimensional energy spectrum in initial conditions program']
switch['O6']   = [True,  'Print plot of normalized spatial distribution of wave heights in initial conditions program (not adapted for distributed memory)']
switch['O7']   = [True,  'Echo input data for homogeneous fields in generic shell']
switch['O7a']  = [False, 'Diagnostic output for output points']
switch['O7b']  = [False, 'Diagnostic output for output points in ww3_multi']
switch['O8']   = [False, 'Filter field output for extremely small wave heights in wave model (useful for some propagation tests)']
switch['O9']   = [False, 'Assign a negative wave height to negative energy in wave model (used in testing phase of new propagation schemes)']
switch['O10']  = [False, 'Identify main elements of multi-grid model extensions in standard output']
switch['O11']  = [False, 'Additional log output on management algorithm in log.mww3']
switch['O12']  = [False, 'Identify removed boundary points in overlapping grids (center)']
switch['O13']  = [False, 'Identify removed boundary points in overlapping grids (edge)']
switch['O14']  = [False, 'Generate log file with buoy data buoy log.ww3 for output type ITYPE = 0 in ww3_outp']
switch['O15']  = [False, 'Generate log file with time stamps of input data file times.XXX in ww3_prep']
switch['O16']  = [False, 'Generate GrADS output of grid partitioning in ww3_gspl']

# OpenMP parallelization
switch['OMPG'] = [False, 'General loop parallelization directives used for both exclusive OpenMP parallelization and hybrid MPI-OpenMP parallelization']
switch['OMPX'] = [False, 'Loop parallelization directives used for exclusive OpenMP parallelization']
switch['OMPH'] = [False, 'Loop parallelization directives used for hybrid MPI-OpenMP parallelization']

# Coninuously moving grids
switch['MGP']  = [False, 'Activate propagation correction']
switch['MGW']  = [False, 'Apply wind correction in moving grid approach']
switch['MGG']  = [False, 'Activate GSE alleviation correction']

# Complier dependent switches
switch['C90']  = [False, 'Compiler directives for Cray C90 (vectorization)']
switch['NEC']  = [False, 'Compiler directives for NEC SX6/SX8 (vectorization)']

# Miscellaneous
switch['ARC']  = [False, 'Arctic grid option for SMC grid']
switch['COU']  = [False, 'Activates the calculation of variables required for coupling']
switch['DSS0'] = [False, 'Switch off frequency dispersion in diffusive dispersion correction']
switch['FLD1'] = [False, 'Sea-state dependent tau Reichl et al. (2014)']
switch['FLD2'] = [False, 'Sea-state dependent tau Donelan et al. (2012)']
switch['IG1']  = [False, 'Second-order spectrum and free infragravity waves']
switch['MLIM'] = [True,  'Use Miche-style shallow water limiter']
switch['MPIBDI']=[False, 'Experimental parallelization of multi-grid model initialization']
switch['MPIT'] = [False, 'Test output for MPI initializations']
switch['MPRF'] = [False, 'Profiling of individual models and nesting in ww3 multi']
switch['NC4']  = [True,  'Activates the NetCDF-4 API in the NetCDF pre- and post-processing programs']
switch['NCC']  = [False, 'NCEP Coupler']
switch['NCO']  = [False, 'Code modifications for operational implementation at NCO (NCEP Central Operations)']
switch['NLS']  = [False, 'Activate nonlinear smoother']
switch['NNT']  = [False, 'Generate file test data nnn.ww3 with spectra and nonlinear interactions for training and testing of NNIA']
switch['OASIS']= [False, 'Initializes OASIS Coupler']
switch['OASACM']=[False, 'OASIS atmospheric model coupling fields']
switch['OASOCM']=[False, 'OASIS oceanic model coupling fields']
switch['REFRX']= [False, 'Enables refraction based on spatial gradients in phase velocity']
switch['REFT'] = [False, 'Test output for shoreline reflection (which is activated with ref1)']
switch['RTD']  = [False, 'Rotated grid option']
switch['RWND'] = [False, 'Correct wind speed for current velocity']
switch['S']    = [False, 'Enable subroutine tracing in the main WAVEWATCH III subroutines']
switch['SCRIP']= [False, 'Enable SCRIP remapping routines']
switch['SCRIPNC']=[False,'Enable storage of remapping weights in NetCDF files']
switch['SEC1'] = [False, 'Enable the use of global time steps less than 1 s, but does not allow output at time steps less than 1 s']
switch['SMC']  = [False, 'Activate SMC grid']
switch['T']    = [False, 'Enable test output throughout the program(s)']
switch['Tn']   = [False, 'Enable test output throughout the program(s)']
switch['TDYN'] = [False, 'Dynamic increment of swell age in diffusive dispersion correction (test cases only)']
switch['TIDE'] = [False, 'Enables tidal analysis']
switch['TIDET']= [False, 'test output for tidal analysis']
switch['TRKNC']= [False, 'Activates the NetCDF API in the wave system tracking post-processing program']
switch['XW0']  = [False, 'Swell diffusion only in ULTIMATE QUICKEST scheme']
switch['XW1']  = [False, 'Wave growth diffusion only in ULTIMATE QUICKEST scheme']



  
         
