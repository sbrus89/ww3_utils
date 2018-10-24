import glob
import collections
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')


runs = collections.OrderedDict()
runs['1/4 degree']='/users/sbrus/scratch4/WW3_timing/glo_15m/'
runs['1/2 degree']='/users/sbrus/scratch4/WW3_timing/glo_30m/'
runs['1 degree'] ='/users/sbrus/scratch4/WW3_timing/glo_1d/'


##################################################################
##################################################################

def get_timing_information(ofile):
  
  f = open(ofile).read().splitlines()

  init_time = None
  elapsed_time = None
  for line in f:
    if line.find('Initialization time :') >= 0:
      init_time = float(line.split()[3])
    if line.find('Elapsed time        :') >= 0:
      elapsed_time = float(line.split()[3])
  if init_time != None and elapsed_time != None:
    return init_time,elapsed_time
  else:
    print 'Timing information not found in screen output file'
    raise SystemExit(0)

##################################################################
##################################################################

def plot_data(xdata,ydata,plot_type,xlabel,ylabel,filename):

  plt.figure()
  lines = []
  labels = []
  for run in xdata: 
    if plot_type == 'loglog':
      l, = plt.loglog(xdata[run],ydata[run],'-o')
    elif plot_type == 'semilogx':
      l, = plt.semilogx(xdata[run],ydata[run],'-o')
    lines.append(l)
    labels.append(run)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.legend(lines,labels)
  plt.tight_layout()
  plt.savefig(filename,bbox_inches='tight')
  plt.close()
      

##################################################################
##################################################################

if __name__ == '__main__':

  # Get scaling data for each run
  data = {}
  for run in runs:
   
      # Initialize the data structure
      data[run] = {'np':[],
                   'wct':[]}

      # Find directories for each core count
      direc = runs[run]
      nps = glob.glob(direc+'np*')

      # Get scaling data for each core count
      for np_dir in nps:
          print np_dir
         
          # Find screen output files
          ofiles = glob.glob(np_dir+'/*.o*')

          # Average run time over all jobs ran
          wct_avg = 0.0      
          for f in ofiles:
              init, elapsed = get_timing_information(f)        
              wct = elapsed - init
              print wct
              wct_avg  = wct_avg  + wct 
          wct_avg  = wct_avg /float(len(ofiles))
          print "average: ", wct_avg

          # Add result to the data structure
          data[run]['wct'].append(wct_avg)
          data[run]['np'].append(float(np_dir.split('np')[1]))
      
      # Sort the data by core counts
      data[run]['np'] = np.array(data[run]['np'])
      data[run]['wct']= np.array(data[run]['wct'])
      ind = data[run]['np'].argsort()
      data[run]['np'] = data[run]['np'][ind]
      data[run]['wct']= data[run]['wct'][ind]   

  # Plot scaling curve for each run
  x = collections.OrderedDict()
  y = collections.OrderedDict()
  for run in runs:
    x[run] = data[run]['np']
    y[run] = data[run]['wct']
  plot_data(x,y,'loglog','number of processors','wall clock time (s)','wct.png')

  # Plot speedup for each run
  x = collections.OrderedDict()
  y = collections.OrderedDict()
  for run in runs:
    n = len(data[run]['np'])-1
    x[run] = np.zeros(n)
    y[run] = np.zeros(n)
    for i in range(n):
      x[run][i] = data[run]['np'][i+1]
      y[run][i] = data[run]['wct'][i]/data[run]['wct'][i+1]
  plot_data(x,y,'semilogx','number of processors','parallel speedpup','parallel_speedup.png') 

  # Plot speedup over high res case
  x = collections.OrderedDict()
  y = collections.OrderedDict()
  high_res = '1/4 degree'
  for run in runs:
    if run != high_res:
      x[run] = data[run]['np']
      y[run] = np.divide(data[high_res]['wct'],data[run]['wct'])
  plot_data(x,y,'semilogx','number of processors','speedup over '+high_res,'resolution_speedup.png')

  # Plot simulated years/day
  x = collections.OrderedDict()
  y = collections.OrderedDict()
  for run in runs:
    x[run] = data[run]['np']
    y[run] = (3600.0*24.0)/(data[run]['wct']*365.0)
  plot_data(x,y,'semilogx','number of processors','simulated years per day (estimated)','sim_yrs_per_day.png')
