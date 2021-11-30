import subprocess

ne ='30'
np = '4'


subprocess.call('GenerateCSMesh --alt --res '+ne+
                              ' --file ne'+ne+'.g', 
                              shell=True)

subprocess.call('GenerateVolumetricMesh --in ne'+ne+'.g'+
                                      ' --out ne'+ne+'pg'+np+'.g'+
                                      ' --np '+np+
                                      ' --uniform',
                                      shell=True)

subprocess.call('ConvertExodusToSCRIP --in ne'+ne+'pg'+np+'.g'
                                    ' --out ne'+ne+'pg'+np+'.scrip.nc',
                                    shell=True)
