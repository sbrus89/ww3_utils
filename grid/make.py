import pprint

targets = {'rotate':['rotate.F90',  'rotate_mod.F90',  'read_write_gmsh.F90',  'globals.F90',  'grid_file_mod.F90'],
           'rotate_test':['rotate_test.F90', 'rotate_mod.F90'],
           'gmsh2vtk':['gmsh2vtk.F90', 'globals.F90', 'grid_file_mod.F90', 'read_write_gmsh.F90', 'write_vtk.F90'],
           'cull_waves_mesh':['cull_waves_mesh.F90', 'kdtree2.F90', 'in_cell_mod.F90','write_vtk.F90'],
           'clean_mesh':['globals.F90', 'grid_file_mod.F90', 'edge_connectivity_mod.F90', 'fix_elements.F90', 'clean_mesh.F90']}

FC = 'ifort'
FFLAGS = '-C -g -traceback'
LIBRARY_LINKS = ['netcdf','netcdff']
LIBRARY_PATHS = ['/blues/gpfs/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-17.0.0/netcdf-fortran-4.4.4-urmb6ss/libs']
INCLUDE_PATH = ['/blues/gpfs/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-17.0.0/netcdf-fortran-4.4.4-urmb6ss/include']


modules = {}
for target in targets:
  for fname in targets[target]:

      if fname not in modules:
        modules[fname] = {'defined':[],'depends':[]}

        with open(fname,'r') as f:
          source = f.read().lower().splitlines()          # make all lowercase to aviod case sensitivity

        for line in source:

          # Find module dependencies by finding USE statements
          if line.strip().find('use ') == 0:
            mod_name = line.split()[1].strip(',')         # strip the , to account for: USE module, ONLY: ....
            if mod_name not in modules[fname]['depends']:
              modules[fname]['depends'].append(mod_name)

          # Find module definitions by finding MODULE statements
          if line.strip().find('module ') == 0:  
            modules[fname]['defined'].append(line.split()[1])

pprint.pprint(modules)
        



makefile = []
makefile.append('FC := '+FC)
makefile.append('FFLAGS := '+FFLAGS)
lib = 'LIB :='
for path in LIBRARY_PATHS:
   lib = lib + ' -L'+path
for link in LIBRARY_LINKS:
  lib = lib + ' -l'+link
makefile.append(lib)
inc = 'INC :='
for path in INCLUDE_PATH:
  inc = inc + ' -I'+path
makefile.append(inc)
makefile.append('')
makefile.append('ODIR = odir/')
makefile.append('INC += -I$(ODIR)')
makefile.append('')
makefile.append('########################################################################')
makefile.append('# Object lists')
makefile.append('########################################################################')
makefile.append('')
for target in targets:
  line = target+'_objects = '
  for f in targets[target]:
    line = line + ' '+f.split('.')[0]+'.o'
  makefile.append(line)
makefile.append('')
for target in targets:
  line = target+'_obj = $(patsubst %.o, $(ODIR)%.o, $('+target+'_objects))'
  makefile.append(line)
makefile.append('')
makefile.append('########################################################################')
makefile.append('# Main Executable Targets')
makefile.append('########################################################################')
makefile.append('')
for target in targets:
  makefile.append('.PHONY : '+target)
makefile.append('')
makefile.append('$(ODIR) :')
makefile.append('\tmkdir -p $@')
makefile.append('')
for target in targets:
  makefile.append(target+' : $(ODIR) $('+target+'_obj)')
  makefile.append('\t$(FC) $(FFLAGS) -o $@ $('+target+'_obj) $(LIB) $(INC)')
  makefile.append('')
makefile.append('')
makefile.append('########################################################################')
makefile.append('# Object File Targets')
makefile.append('########################################################################')
makefile.append('')
files_written = []
for fname in modules:
  if fname not in files_written:
    line = '$(ODIR)'+fname.split('.')[0]+'.o : '+fname
    for dep in modules[fname]['depends']:
      for f in modules:
        for mod in modules[f]['defined']:
          if mod == dep and fname != f:
            line = line + ' $(ODIR)'+f.split('.')[0]+'.o'
    makefile.append(line)
    makefile.append('\t$(FC) $(FFLAGS) $(INC) -c $< -o $@')
    for mod in modules[fname]['defined']:
      makefile.append('\tmv '+mod+'.mod $(ODIR)')
    makefile.append('')
    files_written.append(fname)
  

with open('Makefile','w') as f:
  f.write('\n'.join(makefile))





 
