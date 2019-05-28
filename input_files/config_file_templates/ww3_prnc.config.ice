# Field types:
#   IC1 - Ice thickness
#   IC5 - Ice floe mean diameter
#   ICE - Ice concentration
#   ISI - Icebergs and sea-ice
#   LEV - Water levels
#   WND - Winds
#   WNS - Winds, including air-sea temperature difference
#   CUR - Currents
#   DAT - Data for assimilation
field_type  : "'ICE'"

# Format types:
#   AI - Transfer field as is
#   LL - Field defined on regular lon-lat or Cartesian grid
format_type : "'LL'"

# Include time in file
time_flag   : "T"

# Include header in fille
header_flag : "T"

# Name of spatial dimensions in nc file
lon_dimen   : "lon"
lat_dimen   : "lat"

# Name of variables to use in nc file
u_var       : "ICEC_L1"
v_var       : ""

# File contianing forcing data
filename    : "'ice.nc'"
