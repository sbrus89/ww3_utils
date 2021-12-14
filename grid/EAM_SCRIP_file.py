import subprocess

ne ='30'
np = '4'

cmds = [
         'GenerateCSMesh --alt --res '+ne+
                       ' --file ne'+ne+'.g', 
         
         'GenerateVolumetricMesh --in ne'+ne+'.g'+
                               ' --out ne'+ne+'pg'+np+'.g'+
                               ' --np '+np+
                               ' --uniform',
         
         'ConvertExodusToSCRIP --in ne'+ne+'pg'+np+'.g'
                                             ' --out ne'+ne+'pg'+np+'.scrip.nc',
       ]

for cmd in cmds:
  print(cmd)
  subprocess.call(cmd,shell=True)
