
import subprocess
import submission_script
import os

pes = [36,72,144,288,576,1152,1800]
time = 16.0
queue = 'standard'
#queue = 'interactive'

if __name__ == '__main__':

  pwd = os.getcwd()


  for i,np in enumerate(pes):
     
    # Create run directory
    np_dir = 'np'+str(np)
    subprocess.call('mkdir '+np_dir,shell=True)
    os.chdir(np_dir)
    ln_files = ['wind.ww3','ww3_shel.inp','mod_def.ww3','*.in','restart.ww3','ww3_shel']
    for f in ln_files:
      subprocess.call('ln -sf ../'+f+' .',shell=True)

    # Create information for submission script
    job_name = 'ww3_np'+str(np)
    sub_file = job_name+'.sub'
    pre_cmds = ['cd '+np_dir]

    # Write the submission script
    submission_script.write_submission_script(machine='grizzly',
                                              time=time,
                                              ncores=np,
                                              job_name=job_name,
                                              queue=queue,
                                              exe='ww3_shel',
                                              filename=sub_file,
                                              pre_cmds=pre_cmds)
    if i > 0:
      time = time/2.0

    run_cmd = ['sbatch',sub_file]
    print(' '.join(run_cmd))
    output = subprocess.Popen(run_cmd, stdout=subprocess.PIPE).communicate()[0]
    print(output)

    os.chdir('..')


    
