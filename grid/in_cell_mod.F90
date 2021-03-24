      MODULE in_cell_mod

      USE globals, ONLY: rp
      USE kdtree2_module
      USE netcdf

      IMPLICIT NONE

      INTEGER :: ocean_ncid
      INTEGER :: nCells_dimid, nVertices_dimid, maxEdges_dimid
      INTEGER :: cellsOnVertex_varid, lonCell_varid, latCell_varid
      INTEGER :: verticesOnCell_varid, nEdgesOnCell_varid
      INTEGER :: lonVertex_varid, latVertex_varid

      INTEGER :: nCells, nVertices
      INTEGER :: maxEdges
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: verticesOnCell
      INTEGER, DIMENSION(:), ALLOCATABLE :: nEdgesOnCell
      REAL(rp), DIMENSION(:), ALLOCATABLE :: lonCell, latCell
      REAL(rp), DIMENSION(:), ALLOCATABLE :: lonVertex, latVertex

      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: lonlatCells
      REAL(rp), DIMENSION(:), ALLOCATABLE :: area

      TYPE(kdtree2), POINTER :: tree_xy
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest
      INTEGER :: srchdp
      REAL(rp), PARAMETER :: tol = 1d-7

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE in_cell_init(ocean_mesh_file)

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: ocean_mesh_file

      INTEGER :: i

      PRINT*, "Initializing in_cell"

      ! Open ocean mesh
      CALL check(NF90_OPEN(ocean_mesh_file, NF90_NOWRITE, ocean_ncid))

      ! Get dimension IDs
      CALL check(NF90_INQ_DIMID(ocean_ncid, 'nCells', nCells_dimid))
      CALL check(NF90_INQ_DIMID(ocean_ncid, 'nVertices', nVertices_dimid))
      CALL check(NF90_INQ_DIMID(ocean_ncid, 'maxEdges', maxEdges_dimid))

      ! Get dimension values for ocean mesh
      CALL check(NF90_INQUIRE_DIMENSION(ocean_ncid, nCells_dimid, len=nCells))
      CALL check(NF90_INQUIRE_DIMENSION(ocean_ncid, nVertices_dimid, len=nVertices))
      CALL check(NF90_INQUIRE_DIMENSION(ocean_ncid, maxEdges_dimid, len=maxEdges))

      ! Get variables IDs
      CALL check(NF90_INQ_VARID(ocean_ncid, 'verticesOnCell', verticesOnCell_varid))
      CALL check(NF90_INQ_VARID(ocean_ncid, 'nEdgesOnCell', nEdgesOnCell_varid))
      CALL check(NF90_INQ_VARID(ocean_ncid, 'lonCell', lonCell_varid))
      CALL check(NF90_INQ_VARID(ocean_ncid, 'latCell', latCell_varid))
      CALL check(NF90_INQ_VARID(ocean_ncid, 'lonVertex', lonVertex_varid))
      CALL check(NF90_INQ_VARID(ocean_ncid, 'latVertex', latVertex_varid))

      ! Get variable values for oean mesh
      ALLOCATE(verticesOnCell(maxEdges, nCells))
      CALL check(NF90_GET_VAR(ocean_ncid, verticesOnCell_varid, verticesOnCell))
      ALLOCATE(nEdgesOnCell(nCells))
      CALL check(NF90_GET_VAR(ocean_ncid, nEdgesOnCell_varid, nEdgesOnCell))
      ALLOCATE(lonCell(nCells),latCell(nCells))
      CALL check(NF90_GET_VAR(ocean_ncid, lonCell_varid, lonCell))
      CALL check(NF90_GET_VAR(ocean_ncid, latCell_varid, latCell))
      ALLOCATE(lonVertex(nVertices),latVertex(nVertices))
      CALL check(NF90_GET_VAR(ocean_ncid, lonVertex_varid, lonVertex))
      CALL check(NF90_GET_VAR(ocean_ncid, latVertex_varid, latVertex))
      CALL check(NF90_CLOSE(ocean_ncid))

      ! Initialize kd-tree based on ocean cell centers
      srchdp = 20
      ALLOCATE(closest(srchdp))
      ALLOCATE(lonlatCells(2,nCells))
      DO i = 1,nCells
        lonlatCells(1,i) = lonCell(i)
        lonlatCells(2,i) = latCell(i)
      ENDDO
      tree_xy => kdtree2_create(lonlatCells(1:2,1:nCells), rearrange=.true., sort=.true.)

      ! Compute ocean cell areas
      ALLOCATE(area(nCells))
      DO i = 1,nCells
        CALL compute_total_area(i,lonCell(i),latCell(i),area(i))
      ENDDO 

      PRINT*, "Done"

      RETURN
      END SUBROUTINE in_cell_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE in_cell(xy,in_flag)

      IMPLICIT NONE

      REAL(rp), INTENT(IN) :: xy(2)
      INTEGER, INTENT(OUT) :: in_flag

      INTEGER :: i
      INTEGER :: cell
      REAL(rp) :: area_sum

      CALL kdtree2_n_nearest(tp=tree_xy,qv=xy,nn=srchdp,results=closest)

      in_flag = 0
      DO i = 1,srchdp

        cell = closest(i)%idx
        CALL compute_total_area(cell,xy(1),xy(2),area_sum)

        IF (abs(area_sum-area(cell)) < tol) THEN
          in_flag = 1
        ENDIF

      ENDDO

      RETURN
      END SUBROUTINE in_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE compute_total_area(cell,x3,y3,total_area)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: cell
      REAL(rp), INTENT(IN) :: x3
      REAL(rp), INTENT(IN) :: y3
      REAL(rp), INTENT(OUT) :: total_area

      INTEGER :: j
      INTEGER :: n1,n2
      REAL(rp) :: x1,x2
      REAL(rp) :: y1,y2

      total_area = 0d0
      DO j = 1,nEdgesOnCell(cell)
        n1 = mod(j+0,nEdgesOnCell(cell))+1
        n2 = mod(j+1,nEdgesOnCell(cell))+1

        x1 = lonVertex(verticesOnCell(n1,cell))
        y1 = latVertex(verticesOnCell(n1,cell))

        x2 = lonVertex(verticesOnCell(n2,cell))
        y2 = latVertex(verticesOnCell(n2,cell))
        
        total_area = total_area + 0.5d0*abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
      ENDDO

      RETURN
      END SUBROUTINE compute_total_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check(status)

      USE netcdf

      IMPLICIT NONE
      INTEGER :: status

      IF (status /= NF90_NOERR) THEN
        PRINT("(A, A)"),  "fatal error from ",  TRIM(NF90_STRERROR(status))
      ENDIF

      RETURN
      END SUBROUTINE check

      END MODULE in_cell_mod
