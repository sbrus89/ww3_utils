
import subprocess
import submission_script
import os

pes = [36,72,144,288,576,1152]

if __name__ == '__main__':

  pwd = os.getcwd()


  for i,np in enumerate(pes):

    # Create information for submission script
    job_name = 'ww3_np'+str(np)
    sub_file = job_name+'.sub'
    np_dir = 'np'+str(np)
    post_cmds = ['mkdir '+np_dir,
                 'mv '+job_name+'.o* '+np_dir,
                 'mv '+job_name+'.e* '+np_dir]

    # Write the submission script
    submission_script.write_submission_script(machine='grizzly',
                                              ncores=np,
                                              job_name=job_name,
                                              queue='interactive',
                                              exe='ww3_shel',
                                              filename=sub_file,
                                              post_cmds=post_cmds)

    # Submit the runs with depenencies
    if os.path.isfile(pwd+'/ww3_shel'):
      if i > 0:
        run_cmd = ['sbatch','--dependency=afterok:'+job_id,sub_file]
      else:
        run_cmd = ['sbatch',sub_file]

      print ' '.join(run_cmd)
      output = subprocess.Popen(run_cmd, stdout=subprocess.PIPE).communicate()[0]
      print output

      job_id = output.split()[-1]

    
