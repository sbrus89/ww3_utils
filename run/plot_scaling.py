import glob
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')


runs ={'glo_15m':'/users/sbrus/scratch4/WW3_timing/glo_15m/',
       'glo_30m':'/users/sbrus/scratch4/WW3_timing/glo_30m/',
       'glo_1d' :'/users/sbrus/scratch4/WW3_timing/glo_1d/'}


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

  # Plot data for each run
  plt.figure()
  lines = []
  labels = []
  for run in runs:
    l, = plt.loglog(data[run]['np'],data[run]['wct'],'-o')
    lines.append(l)
    labels.append(run)
  plt.xlabel('number of processors')
  plt.ylabel('time (s)')
  plt.legend(lines,labels)
  plt.tight_layout()
  plt.savefig('wct.png',bbox_inches='tight')
  plt.close()
      
    
