import define_switches 

switches = define_switches.switch

f = open('switch','w')

for switch in switches:
    if switches[switch][0] == True:
        f.write(switch + ' ')
        print switch + ':' + (8-len(switch))*' ' + switches[switch][1]
