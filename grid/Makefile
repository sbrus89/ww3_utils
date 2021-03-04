FC := ifort
FFLAGS := -C -g -traceback
LIB := -L/blues/gpfs/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-17.0.0/netcdf-fortran-4.4.4-urmb6ss/libs -lnetcdf -lnetcdff
INC := -I/blues/gpfs/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-17.0.0/netcdf-fortran-4.4.4-urmb6ss/include

ODIR = odir/
FFLAGS += -I$(ODIR)

########################################################################
# Object lists
########################################################################

rotate_objects =  rotate.o rotate_mod.o read_write_gmsh.o globals.o grid_file_mod.o
rotate_test_objects =  rotate_test.o rotate_mod.o

rotate_obj = $(patsubst %.o, $(ODIR)%.o, $(rotate_objects))
rotate_test_obj = $(patsubst %.o, $(ODIR)%.o, $(rotate_test_objects))

########################################################################
# Main Executable Targets
########################################################################

.PHONY : rotate
.PHONY : rotate_test

$(ODIR) :
	mkdir -p $@

rotate : $(ODIR) $(rotate_obj)
	$(FC) $(FFLAGS) -o $@ $(rotate_obj) $(LIB) $(INC)

rotate_test : $(ODIR) $(rotate_test_obj)
	$(FC) $(FFLAGS) -o $@ $(rotate_test_obj) $(LIB) $(INC)


########################################################################
# Object File Targets
########################################################################

$(ODIR)rotate.o : rotate.F90 $(ODIR)rotate_mod.o $(ODIR)read_write_gmsh.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)rotate_mod.o : rotate_mod.F90 $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv rotate_mod.mod $(ODIR)

$(ODIR)read_write_gmsh.o : read_write_gmsh.F90 $(ODIR)grid_file_mod.o $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv read_write_gmsh.mod $(ODIR)

$(ODIR)globals.o : globals.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)grid_file_mod.o : grid_file_mod.F90 $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)

$(ODIR)rotate_test.o : rotate_test.F90 $(ODIR)rotate_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
