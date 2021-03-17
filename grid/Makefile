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
gmsh2vtk_objects =  gmsh2vtk.o globals.o grid_file_mod.o read_write_gmsh.o write_vtk.o
clean_mesh_objects =  globals.o grid_file_mod.o edge_connectivity_mod.o fix_elements.o clean_mesh.o

rotate_obj = $(patsubst %.o, $(ODIR)%.o, $(rotate_objects))
rotate_test_obj = $(patsubst %.o, $(ODIR)%.o, $(rotate_test_objects))
gmsh2vtk_obj = $(patsubst %.o, $(ODIR)%.o, $(gmsh2vtk_objects))
clean_mesh_obj = $(patsubst %.o, $(ODIR)%.o, $(clean_mesh_objects))

########################################################################
# Main Executable Targets
########################################################################

.PHONY : rotate
.PHONY : rotate_test
.PHONY : gmsh2vtk
.PHONY : clean_mesh

$(ODIR) :
	mkdir -p $@

rotate : $(ODIR) $(rotate_obj)
	$(FC) $(FFLAGS) -o $@ $(rotate_obj) $(LIB) $(INC)

rotate_test : $(ODIR) $(rotate_test_obj)
	$(FC) $(FFLAGS) -o $@ $(rotate_test_obj) $(LIB) $(INC)

gmsh2vtk : $(ODIR) $(gmsh2vtk_obj)
	$(FC) $(FFLAGS) -o $@ $(gmsh2vtk_obj) $(LIB) $(INC)

clean_mesh : $(ODIR) $(clean_mesh_obj)
	$(FC) $(FFLAGS) -o $@ $(clean_mesh_obj) $(LIB) $(INC)


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

$(ODIR)gmsh2vtk.o : gmsh2vtk.F90 $(ODIR)globals.o $(ODIR)write_vtk.o $(ODIR)grid_file_mod.o $(ODIR)read_write_gmsh.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)write_vtk.o : write_vtk.F90 $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv write_vtk.mod $(ODIR)

$(ODIR)edge_connectivity_mod.o : edge_connectivity_mod.F90 $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_connectivity_mod.mod $(ODIR)

$(ODIR)fix_elements.o : fix_elements.F90 $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv fix_elements.mod $(ODIR)

$(ODIR)clean_mesh.o : clean_mesh.F90 $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)edge_connectivity_mod.o $(ODIR)fix_elements.o
	$(FC) $(FFLAGS) -c $< -o $@
