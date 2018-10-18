grid_name   = "'glo_15m'"    # Grid name

#-------------------------
# Spctral grid parameters
#-------------------------
freq_inc    = '1.07'         # Frequency increment factor
first_freq  = '0.035'        # Minimum frequency 
num_freq    = '50'           # Number of frequency bins
num_dir     = '36'           # Number of directional bins
rel_dir_off = '0.5'          # Relative offset of first direction

#-------------
# Model flags
#-------------
fldry       = 'F'            # Dry run
flcx        = 'T'            # Activate x component of propagation
flcy        = 'T'            # Activate y component of propagattion
flcth       = 'T'            # Activate direction shifts
flck        = 'F'            # Activate wavenumber shifts
flsou       = 'T'            # Activate source terms

#----------------------
# Time step parameters
#----------------------
max_gl_dt   = '900'          # Max global time step
max_geo_dt  = '300'          # Max x-y time step
max_spec_dt = '450'          # Max k-theta time step
min_src_dt  = '30'           # Min source term time step

#---------------------
# Namelist parameters
#---------------------
namelist = {'MISC':[['FLAGTR','2']]}

#----------------------------
# Grid definition parameters
#----------------------------
gstrg       = "'RECT'"       # Type of grid (RECT = retilinear, CURV = curvilinear, UNST = unstructured)
flagll      = 'T'            # Coordinate system flag (T = spherical, F = Cartesian)
cstrg       = "'SMPL'"       # Closure type (NONE = no closure, SMPL = simple closure, TRPL, tripole closure)
nx          = '1440'         # Number of points in x direction
ny          = '720'          # Number of points in y direction

# Rectilinlinear grid
sx          = '0.25'         # Grid spacing in x direction
sy          = '0.25'         # Grid spacing in y direction
grd_sca_fac = '1.0'          # Scaling factor (division)
x11         = '0.0'          # x coordinate of (1,1)
y11         = '-90.0'        # y coordinate of (1,1)
crd_sca_fac = '1.0'          # Scaling factor (division)

# Bottom depth file
lim_bot_dep = '-0.1'         # Limiting bottom depth
min_dep     = '2.5'          # Minimum water depth
dep_unit    = '20'           # Unit number          
dep_sca_fac = '0.001'        # Scaling factor
dep_idla    = '1'            # Layout indicator (1 = Read line-by-line bottom to top, 2 = Like 1, single read statement, 3 = Read line-by-line top to bottom, 4 = Like 3, single read statement)
dep_idfm    = '1'            # Format indicator (1 = Free format, 2 = Fixed format, 3 = Unformatted)
dep_frmt    = "'(....)'"     # Format 
dep_from    = "'NAME'"       # File type parameter
dep_name    = "'glo_15m.bot'"# Name of file

# Obstruction file
obs_unit    = '21'           # Unit number
obs_sca_fac = '0.01'         # Scale factor (multiplication)
obs_idla    = '1'            # Layout indicator
obs_idfm    = '1'            # Format indicator
obs_frmt    = "'(....)'"     # Format 
obs_from    = "'NAME'"       # File type parameter
obs_name    = "'glo_15m.obst'"# Name of file

# Mask file
msk_unit    = '22'           # Unit number
msk_idla    = '1'            # Layour indicator
msk_idfm    = '1'            # Format indicator
msk_frmt    = "'(....)'"     # Format 
msk_from    = "'NAME'"       # File type parameter
msk_name    = "'glo_15m.mask'"# Name of file



