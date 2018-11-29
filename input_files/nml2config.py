import pprint

string =  r'''&SLN1 CLIN =  80.0, RFPM =  1.00, RFHF =  0.50 
  &SIN3 ZWND = 10.0, ALPHAO = 0.00950, Z0MAX = 0.00000, BETAMAX = 1.25000,
        SINTHP = 2.00000, ZALP = 0.01100, TAUWSHELTER = 0.00000, SWELLFPAR = 0,
        SWELLF = 0.00000, SWELLF2 = 0.00000, SWELLF3 = 0.00000, SWELLF4 = 100000.0, SWELLF5 = 1.20000, Z0RAT = 1.00000 
  &SNL1 LAMBDA =  0.250, NLPROP = 0.278E+08, KDCONV =  0.750, KDMIN =  0.500,
        SNLCS1 =  5.500, SNLCS2 =  0.833, SNLCS3 =  -1.250 
  &SDS3 SDSC1 = -0.2100E+01, SDSC2 =  0.0000E+00, SDSC3 =  0.0000E+00,
        SDSC4 =  0.0000E+00, SDSC5 =  0.0000E+00, SDSC6 =  0.1000E+01, WNMEANP = 0.50,
        SDSDELTA1 = 0.40, SDSDELTA2 = 0.60, SDSLF = 1.00, SDSHF = 1.00, 
        SDSBR =   0.1200E-02, SDSBR2 = 0.80, SDSP = 2.00, SDSISO = 1, SDSCOS =0.0, SDSDTH =180.0, 
        SDSBM0 =  1.00, SDSBM1 = 0.00, SDSBM2 = 0.00, SDSBM3 = 0.00, SDSBM4 = 0.00 
  &SBT1 GAMMA = -0.6700E-01 
  &SDB1 BJALFA =  1.000, BJGAM =  0.730, BJFLAG = .TRUE. 
  &PRO3 CFLTM = 0.70, WDTHCG = 1.50, WDTHTH = 1.50 
  &MISC CICE0 = 0.330, CICEN = 0.670, PMOVE = 0.500, XSEED = 1.000, FLAGTR = 4,
        XP = 0.150, XR = 0.100, XFILT = 0.050, IHM =  100, HSPM = 0.050,
        WSM = 1.700, WSC = 0.333, FLC = .TRUE., FMICHE = 0.750 '''
string = string.replace('\n',' ')
namelists = string.split('&')
pprint.pprint(namelists)
output = '{'
for nml in namelists:
  if nml:
    name = nml.split()[0]
    print name
    output = output + "'"+name+"':["
    parameters = nml.replace(name,'').split(',')
    for i,param in enumerate(parameters):
      if i > 0:
        output = output + ','
      var = param.split('=')[0].strip()
      val = param.split('=')[1].strip()
      output = output + "['"+var+"','"+val+"']"
    output = output + '],\n'
output = output.rstrip(',\n') + '}'

print output
