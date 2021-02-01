      PROGRAM gmsh2vtk

      USE globals, ONLY: rp
      USE write_vtk
      USE grid_file_mod, ONLY: grid_type
      USE read_gmsh, ONLY: read_gmsh_file

      TYPE(grid_type) :: mesh
      INTEGER :: i,ne 
      CHARACTER(20) :: filename,vtkfile
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect

      vtkfile = 'new_orleans.vtk'
      CALL read_gmsh_file('new_orleans.msh',mesh)

      ALLOCATE(ect(3,mesh%ne))
      ne = 0
      DO i = 1,mesh%ne
        IF (mesh%el_type(i) == 2) THEN
          ne = ne + 1
          ect(:,ne) = mesh%ect(:,i)
        ENDIF
      ENDDO
          

      PRINT*, 'Writing '//TRIM(ADJUSTL(vtkfile))
      CALL write_vtk_file(TRIM(ADJUSTL(vtkfile)),mesh%nn,mesh%xy,ne,ect,mesh%depth)

  

      END PROGRAM gmsh2vtk
