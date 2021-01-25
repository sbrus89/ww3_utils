import subprocess

#title = 'T62 grid'
#nlat = 94
#nlon = 192
#outfile = 'T62.SCRIP.nc'

title = 'CFSv2 0.205x0.204 grid'
nlat = 880 
nlon = 1760
outfile = 'CFSv2.SCRIP.nc'

cmd = 'ncremap -G '
cmd = cmd + "ttl='"+title+"'"
cmd = cmd + '#latlon='+str(nlat)+','+str(nlon)
cmd = cmd + '#lat_typ=gss'
cmd = cmd + '#lon_typ=grn_ctr'
cmd = cmd + ' -g '+outfile

print(cmd)
subprocess.call(cmd,shell=True)
