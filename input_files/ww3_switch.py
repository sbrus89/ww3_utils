import os
import sys
sys.dont_write_bytecode = True
import ww3_switch_config  

def write_switch():
  
  pwd = os.getcwd()
  
  switches = ww3_switch_config.switch
  
  f = open(pwd+'/switch','w')
  
  for switch in switches:
      if switches[switch][0] == True:
          f.write(switch + ' ')
          print switch + ':' + (8-len(switch))*' ' + switches[switch][1]

if __name__ == '__main__':
  write_switch()


