import subprocess

#--------------------
# Request parameters
#--------------------
start_date = '2005-06-01 00:00'
end_date   = '2005-07-01 00:00'
data_set   = 'CFSR hourly'
parameters = ['uwind','vwind']
products   = ['1 hr forecast','2 hr forecast','3 hr forecast','4 hr forecast','5 hr forecast','6 hr forecast']
grids      = ['0.5x0.5']
levels     = ['specified height above ground: 10m']
region_box = ['-180','180','-90','90']

#------------------------
# Request paramter codes
#------------------------
data_set_codes  = {'CFSR hourly'  :'093.1',
                   'CFSRv2 hourly':'094.1'}
parameter_codes = {'uwind':'3%217-0.2-1:0.2.2',
                   'vwind':'3%217-0.2-1:0.2.3'}
product_codes   = {'6 hr forecast':'3',
                   '3 hr forecast':'119',
                   '1 hr forecast':'486',
                   '2 hr forecast':'488',
                   '4 hr forecast':'490',
                   '5 hr forecast':'492'}
grid_codes      = {'0.312x0.312':'83',
                   '0.5x0.5'    :'57',
                   '1.875x1.904':'3' ,
                   '2.5x2.5'    :'4'}
level_codes     = {'specified height above ground: 10m':'223'}

#----------------
# Log in command
#----------------
command = ['curl','-c','rda_auth_cookies','-d','"email=sbrus@nd.edu&passwd=irishswimming&action=login"','https://rda.ucar.edu/cgi-bin/login']
subprocess.call(command)

#----------------------
# Data request command
#----------------------
command = ['curl','-b','rda_auth_cookies','-d','"dsid=ds'         +data_set_codes[data_set]                          +'&' \
                                              + 'rtype=S&'                                                                \
                                              + 'rinfo=dsnum='    +data_set_codes[data_set]                          +';' \
                                              + 'startdate='      +start_date                                        +';' \
                                              + 'enddate='        +end_date                                          +';' \
                                              + 'product='        +','.join([product_codes[x]   for x in products])  +';' \
                                              + 'parameters='     +','.join([parameter_codes[x] for x in parameters])+';' \
                                              + 'level='          +','.join([level_codes[x]     for x in levels])    +';' \
                                              + 'grid_definition='+','.join([grid_codes[x]      for x in grids])     +';' \
                                              + 'slat='           +region_box[2]                                     +';' \
                                              + 'nlat='           +region_box[3]                                     +';' \
                                              + 'wlon='           +region_box[0]                                     +';' \
                                              + 'elon='           +region_box[1]                                     +';' \
                                              + 'ofmt=netCDF"','https://rda.ucar.edu/php/dsrqst.php']
command =  ' '.join(command)
subprocess.call(command,shell=True)
