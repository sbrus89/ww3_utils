import os
import define_switches 

pwd = os.getcwd()

switches = define_switches.switch

f = open(pwd+'/switch','w')

for switch in switches:
    if switches[switch][0] == True:
        f.write(switch + ' ')
        print switch + ':' + (8-len(switch))*' ' + switches[switch][1]
