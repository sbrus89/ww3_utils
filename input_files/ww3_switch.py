import os
import sys
sys.dont_write_bytecode = True
import yaml
import pprint

def write_switch():
  
  pwd = os.getcwd()
  
  f = open(pwd+'/ww3_switch.config')
  switches = oyaml.load(f)

  
  f = open(pwd+'/switch','w')
  
  for switch in switches:
      if switches[switch][0] == True:
          f.write(switch + ' ')
          print switch + ':' + (8-len(switch))*' ' + switches[switch][1]

if __name__ == '__main__':
  write_switch()


