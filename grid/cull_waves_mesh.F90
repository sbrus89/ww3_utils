      PROGRAM cull_waves_mesh

      USE netcdf
      USE in_cell_mod, ONLY: in_cell_init, pt_in_cell, check, pi, lonlat2xyz, R
      USE globals, ONLY: rp
      USE write_vtk
      USE read_write_gmsh

      IMPLICIT NONE

      INTEGER :: waves_ncid
      INTEGER :: ocean_ncid
      INTEGER :: nCells_dimid, nVertices_dimid
      INTEGER :: cellsOnVertex_varid, lonCell_varid, latCell_varid
      INTEGER :: bottomDepth_varid

      CHARACTER(100) :: waves_mesh_file
      CHARACTER(100) :: waves_mesh_culled_vtk
      CHARACTER(100) :: waves_mesh_culled_gmsh
      CHARACTER(100) :: ocean_mesh_file

      INTEGER :: ne_waves, nn_waves
      INTEGER :: ne_new, nn_new
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect_waves
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect_new
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x_waves, y_waves
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x_new, y_new
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_waves
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_new
      REAL(rp) :: xy(2)
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xyz
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xy_new

      INTEGER :: ncells_ocean
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_ocean

      INTEGER :: el, nd
      INTEGER :: i
      INTEGER :: in_cell
      INTEGER, DIMENSION(:), ALLOCATABLE :: nepn
      INTEGER, DIMENSION(:), ALLOCATABLE :: new_node_numbers

      INTEGER, DIMENSION(:), ALLOCATABLE :: wave_node_in_ocean
      INTEGER, DIMENSION(:), ALLOCATABLE :: wave_elements_keep

      NAMELIST / inputs / waves_mesh_file, ocean_mesh_file
      NAMELIST / output / waves_mesh_culled_vtk, waves_mesh_culled_gmsh
 
      OPEN(UNIT=15, FILE='cull_waves_mesh.nml')
      READ(UNIT=15, NML=inputs)
      READ(UNIT=15, NML=output)
      CLOSE(15)
 
      PRINT*, 'waves_mesh_file = ', TRIM(ADJUSTL(waves_mesh_file))
      PRINT*, 'ocean_mesh_file = ', TRIM(ADJUSTL(ocean_mesh_file))
      PRINT*, 'waves_mesh_file_vtk = ', TRIM(ADJUSTL(waves_mesh_culled_vtk))
      PRINT*, 'waves_mesh_file_gmsh = ', TRIM(ADJUSTL(waves_mesh_culled_gmsh))
      CALL check(NF90_OPEN(TRIM(ADJUSTL(waves_mesh_file)), NF90_NOWRITE, waves_ncid)) 
      CALL check(NF90_INQ_DIMID(waves_ncid, 'nCells', nCells_dimid))
      CALL check(NF90_INQ_DIMID(waves_ncid, 'nVertices', nVertices_dimid))

      CALL check(NF90_INQUIRE_DIMENSION(waves_ncid, nCells_dimid, len=nn_waves))
      CALL check(NF90_INQUIRE_DIMENSION(waves_ncid, nVertices_dimid, len=ne_waves))

      CALL check(NF90_INQ_VARID(waves_ncid, 'cellsOnVertex', cellsOnVertex_varid))
      CALL check(NF90_INQ_VARID(waves_ncid, 'lonCell', lonCell_varid))
      CALL check(NF90_INQ_VARID(waves_ncid, 'latCell', latCell_varid))

      ALLOCATE(x_waves(nn_waves), y_waves(nn_waves))
      CALL check(NF90_GET_VAR(waves_ncid, lonCell_varid, x_waves))
      CALL check(NF90_GET_VAR(waves_ncid, latCell_varid, y_waves))
      ALLOCATE(ect_waves(3, ne_waves))
      CALL check(NF90_GET_VAR(waves_ncid, cellsOnVertex_varid, ect_waves))
      CALL check(NF90_CLOSE(waves_ncid))

      CALL check(NF90_OPEN(ocean_mesh_file, NF90_NOWRITE, ocean_ncid))
      CALL check(NF90_INQ_DIMID(ocean_ncid, 'nCells', nCells_dimid))
      CALL check(NF90_INQUIRE_DIMENSION(ocean_ncid, nCells_dimid, len=ncells_ocean))
      CALL check(NF90_INQ_VARID(ocean_ncid, 'bottomDepth', bottomDepth_varid))
      ALLOCATE(depth_ocean(ncells_ocean))
      CALL check(NF90_GET_VAR(ocean_ncid, bottomDepth_varid, depth_ocean))
      CALL check(NF90_CLOSE(ocean_ncid))


      
      CALL in_cell_init(TRIM(ADJUSTL(ocean_mesh_file)))


      ALLOCATE(wave_node_in_ocean(nn_waves))
      wave_node_in_ocean = 0
      ALLOCATE(wave_elements_keep(ne_waves))
      wave_elements_keep = 0
      ALLOCATE(depth_waves(nn_waves))
      depth_waves = -10d0

      ne_new = 0
elems:DO el = 1,ne_waves

         IF (mod(el,1000) == 0) THEN
           PRINT*, el
         ENDIF

   nds:  DO i = 1,3
         
          nd = ect_waves(i,el)

          IF (nd == 0) THEN
            wave_elements_keep(el) = 0
            EXIT nds
          ENDIF

          IF (wave_node_in_ocean(nd) == 1) THEN
            wave_elements_keep(el) = 1
            ne_new = ne_new + 1
            CYCLE nds
          ENDIF

          xy(1) = x_waves(nd)
          xy(2) = y_waves(nd)

          CALL pt_in_cell(xy,in_cell)

          IF (in_cell > 0) THEN
            wave_node_in_ocean(nd) = 1
            depth_waves(nd) = depth_ocean(in_cell)
            IF (wave_elements_keep(el) == 0) THEN
              wave_elements_keep(el) = 1
              ne_new = ne_new + 1
            ENDIF
          ENDIF    

        ENDDO nds
      ENDDO elems


      ! Create new element connectivity table
      ALLOCATE(ect_new(3,ne_new))
      ne_new = 0
      DO el = 1,ne_waves
        IF (wave_elements_keep(el) == 1) THEN
          ne_new = ne_new + 1
          DO i = 1,3
            ect_new(i,ne_new) = ect_waves(i,el)
          ENDDO
        ENDIF
      ENDDO
      PRINT*, ne_waves, ne_new

      ! Find number of elements per node
      ALLOCATE(nepn(nn_waves))
      nepn = 0
      DO el = 1,ne_new
        DO i = 1,3
          nd = ect_new(i,el)
          nepn(nd) = nepn(nd) + 1
        ENDDO
      ENDDO

      ! Build new coordinate list
      ALLOCATE(x_new(nn_waves),y_new(nn_waves))
      ALLOCATE(depth_new(nn_waves))
      ALLOCATE(new_node_numbers(nn_waves))
      nn_new = 0
      DO nd = 1,nn_waves 
        IF (nepn(nd) /= 0) THEN
          nn_new = nn_new + 1
          x_new(nn_new) = x_waves(nd)
          y_new(nn_new) = y_waves(nd)
          depth_new(nn_new) = depth_waves(nd)
          new_node_numbers(nd) = nn_new
        ENDIF
      ENDDO

      ! Adjust element connectivity node numbers
      DO el = 1,ne_new
        DO i = 1,3
          ect_new(i,el) = new_node_numbers(ect_new(i,el))
        ENDDO
      ENDDO

      ! Convert to Cartesian coordinates for vtk output
      ALLOCATE(xyz(3,nn_new))
      ALLOCATE(xy_new(2,nn_new))
      DO i = 1,nn_new
        CALL lonlat2xyz(x_new(i),y_new(i),R,xyz(1,i),xyz(2,i),xyz(3,i))
        xy_new(1,i) = x_new(i)
        xy_new(2,i) = y_new(i)
      ENDDO
     
      CALL write_vtk_file(TRIM(ADJUSTL(waves_mesh_culled_vtk)),nn_new,xyz,ne_new,ect_new,depth_new)

      CALL write_header(TRIM(ADJUSTL(waves_mesh_culled_gmsh)))
      CALL write_nodes(nn_new,xy_new,depth_new)
      CALL write_elements(ne_new,ect_new)

      END PROGRAM cull_waves_mesh
