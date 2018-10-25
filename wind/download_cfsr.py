import subprocess

#--------------------
# Request parameters
#--------------------
data = 'ssh'
start_date = '2005-07-01 00:00'
end_date   = '2005-08-01 00:00'
data_set   = 'CFSR hourly'
region_box = ['-180','180','-90','90']

if data == 'wind':
  parameters = ['uwind','vwind']
  products   = ['1 hr forecast',
                '2 hr forecast',
                '3 hr forecast',
                '4 hr forecast',
                '5 hr forecast',
                '6 hr forecast']
  levels     = ['specified height above ground: 10m']
  grids      = ['0.5x0.5 (720x361)']
elif data == 'currents':
  parameters = ['ucurrent','vcurrent']
  products   = ['1 hr average (initial+0 to initial+1)',
                '1 hr average (initial+1 to initial+2)',
                '1 hr average (initial+2 to initial+3)',
                '1 hr average (initial+3 to initial+4)',
                '1 hr average (initial+4 to initial+5)',
                '1 hr average (initial+5 to initial+6)']
  levels     = ['depth below sea level: 5m']
  grids      = ['0.5x0.5 (720x360)']
elif data == 'ssh':
  parameters = ['ssh']
  products   = ['1 hr average (initial+0 to initial+1)',
                '1 hr average (initial+1 to initial+2)',
                '1 hr average (initial+2 to initial+3)',
                '1 hr average (initial+3 to initial+4)',
                '1 hr average (initial+4 to initial+5)',
                '1 hr average (initial+5 to initial+6)']
  levels     = ['ground or water surface']
  grids      = ['0.5x0.5 (720x360)']

#------------------------
# Request paramter codes
#------------------------
data_set_codes  = {'CFSR hourly'  :'093.1',
                   'CFSRv2 hourly':'094.1'}
parameter_codes = {'uwind':'3%217-0.2-1:0.2.2',
                   'vwind':'3%217-0.2-1:0.2.3',
                   'ucurrent':'3%217-4.2-1:10.1.2',
                   'vcurrent':'3%217-4.2-1:10.1.3',
                   'ssh':'3%217-4.2-1:10.3.195'}
product_codes   = {'6 hr forecast':'3',
                   '3 hr forecast':'119',
                   '1 hr forecast':'486',
                   '2 hr forecast':'488',
                   '4 hr forecast':'490',
                   '5 hr forecast':'492',
                   '1 hr average (initial+0 to initial+1)':'944',
                   '1 hr average (initial+1 to initial+2)':'945',
                   '1 hr average (initial+2 to initial+3)':'946',
                   '1 hr average (initial+3 to initial+4)':'947',
                   '1 hr average (initial+4 to initial+5)':'948',
                   '1 hr average (initial+5 to initial+6)':'949'}
grid_codes      = {'0.312x0.312'       :'83',
                   '0.5x0.5 (720x361)' :'57',
                   '0.5x0.5 (720x360)' :'63',
                   '1.875x1.904'       :'3' ,
                   '2.5x2.5'           :'4'}
level_codes     = {'specified height above ground: 10m':'223',
                   'depth below sea level: 5m'         :'128',
                   'ground or water surface'           :'107,210'}

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
print command
subprocess.call(command,shell=True)
