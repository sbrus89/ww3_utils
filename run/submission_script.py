import datetime

def write_submission_script(machine,ncores,job_name,queue,exe=None,time=None,filename='run.sub',account=None,email='sbrus@anl.gov',pre_cmds=[],post_cmds=[],compiler='intel'):

  # Machine specific settings
  if machine == "grizzly":
    cores_per_node = 36
    if compiler == 'gnu':
      module_load = ['gcc/5.3.0',
                     'openmpi/1.10.5',
                     'netcdf/4.4.1 ',
                     'parallel-netcdf/1.5.0 ',
                     'pio/1.7.2']
    elif compiler == 'intel':
      module_load = ['intel/17.0.1',
                     'mvapich2/2.2',
                     'netcdf/4.4.1 ',
                     'parallel-netcdf/1.5.0 ',
                     'pio/1.7.2']
    else:
      print('compiler not supported')
      raise SystemExit(0)
    module_purge = True
    module_use = ['/usr/projects/climate/SHARED_CLIMATE/modulefiles/all/']
    queues = {'standard'   :{'np_min':1   ,'np_max':53640,'t_lim':16.0},
              'interactive':{'np_min':1   ,'np_max':2520 ,'t_lim':4.0 }}
    if account == None:
      account = 'e3sm'
  elif machine == "badger":
    cores_per_node = 36
    if compiler == 'gnu':
      module_load = ['gcc/7.4.0',
                     'openmpi/2.1.2',
                     'hdf5-serial/1.8.16',
                     'netcdf-serial/4.4.0']
    else:
      print('compiler not supported')
      raise SystemExit(0)
    module_purge = True
    module_use = []
    queues = {'standard'   :{'np_min':1,'np_max':3600,'t_lim':16.0},
              'interactive':{'np_min':1,'np_max':72  ,'t_lim':2.0 }}
    if account == None:
      account = 'e3sm'
  elif machine == "blues":
    cores_per_node = 36
    if compiler == 'intel':
      module_load = ['intel/17.0.0-pwabdn2',
                     'mvapich2/2.2-verbs-qwuab3b',
                     'netcdf/4.4.1-tckdgwl',
                     'netcdf-fortran/4.4.4-urmb6ss']
    else:
      print('compiler not supported')
      raise SystemExit(0)
    module_purge = True
    module_use = []
    queues = {'acme-small'   :{'np_min':1,'np_max':180 ,'t_lim':48.0},
              'acme-medium'  :{'np_min':181,'np_max':2160,'t_lim':24.0 },
              'acme-large'   :{'np_min':2161,'np_max':8640,'t_lim':12.0 }}
    partition = queue
    if account == None:
      account = 'condo'
  
  
  # Checks
  if queue not in queues:
    print("queue not supported on " + machine)
    raise SystemExit(0)

  if ncores > 1 and ncores % cores_per_node != 0:
    print("number of cores should be multple of the number of cores per node")
    raise SystemExit(0)

  if ncores > queues[queue]['np_max'] or ncores < queues[queue]['np_min']:
    print("number of cores requested is not in allowed range for queue")
    raise SystemExit(0)
  
  if time:
    if time > queues[queue]['t_lim']:
      print("time exceeds time limit of queue")
      raise SystemExit(0)
  else:
    time = queues[queue]['t_lim']  


  # Calculate some parameters
  nnodes = int(ncores/cores_per_node)
  if nnodes < 1:
    nnodes = 1
  time = str(datetime.timedelta(hours=time))
  
  if machine == 'blues':
    queue = None

  # SBATCH options  
  sbatch = {'nodes'    : str(nnodes),
            'time'     : time,
            'account'  : account,
            'job-name' : job_name,
            'output'   : job_name+'.o%j',
            'error'    : job_name+'.e%j',
            'qos'      : queue,
            'partition': partition,
            'mail-user': email,
            'mail-type': 'all'}

  # Write submission script header
  f = open(filename,'w')
  f.write('#!/bin/bash\n')
  for opt in sbatch:
    if sbatch[opt] != None:
      f.write('#SBATCH --'+opt+'='+sbatch[opt]+'\n')
  f.write('\n')

  # Write submission script module loads
  f.write('# Modules for '+machine+'\n')
  if module_purge:
    f.write('module purge\n')
  for mod in module_use:
    f.write('module use '+mod+'\n')
  for mod in module_load:
    f.write('module load '+mod+'\n')
  f.write('\n')

  # Write submission script commands
  for cmd in pre_cmds:
    f.write(cmd+'\n')
  if exe:
    f.write('srun -n '+str(ncores)+' ./'+exe+'\n')
  for cmd in post_cmds:
    f.write(cmd+'\n')
  f.close()



#########################################################################################
#########################################################################################


if __name__ == '__main__':
  import argparse
  import socket

  # Command line arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("--machine"   ,type=str,   help="machine name"          , default=None)
  parser.add_argument("--cores"     ,type=int,   help="number of cores"       , required=True)
  parser.add_argument("--job_name"  ,type=str,   help="job name"              , default='run')
  parser.add_argument("--queue"     ,type=str,   help="queue name"            , default='interactive')
  parser.add_argument("--exe"       ,type=str,   help="main executable name"  , default=None)
  parser.add_argument("--time"      ,type=float, help="time limit (hours)"    , default=None)
  parser.add_argument("--file_name" ,type=str,   help="submission script name", default='run.sub')
  parser.add_argument("--account"   ,type=str,   help="name of allocation    ", default=None)
  parser.add_argument("--email"     ,type=str,   help="email to send run information", default='sbrus@anl.gov')
  parser.add_argument("--pre_cmds"  ,type=str,   help="commands to run before main executable", nargs='*',default=[])
  parser.add_argument("--post_cmds" ,type=str,   help="commands to run after main executable" , nargs='*',default=[])
  parser.add_argument("--compiler"   ,type=str,   help="compiler name"          , default='intel')
  args = parser.parse_args()

  # Determine machine
  if args.machine:
    machine = args.machine
  else:
    host = socket.gethostname()
    if host[0:2] == 'gr':
      machine = 'grizzly'
    elif host[0:2] == 'wf':
      machine = 'wolf'
    elif host[0:2] == 'bl':
      machine = 'blues'

  # Write the submission script
  write_submission_script(machine,
                          args.cores,
                          args.job_name,
                          args.queue,
                          args.exe,
                          args.time,
                          args.file_name,
                          args.account,
                          args.email,
                          args.pre_cmds,
                          args.post_cmds,
                          args.compiler)
 
 
