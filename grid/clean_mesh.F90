      PROGRAM clean_mesh

      USE globals, ONLY: rp,nverts
      USE grid_file_mod, ONLY: read_header,read_coords,read_connectivity,&
                               write_header,write_coords,write_connectivity, &
                               write_open_boundaries, write_flow_boundaries
      USE edge_connectivity_mod, ONLY: elements_per_node,find_edge_pairs,find_interior_edges, &
                                       find_element_edges,find_neighbor_elements,find_boundary_nodes
      USE fix_elements, ONLY: flag_problem_elements,flag_isolated_element_patches, &
                              flag_single_element_passages, &
                              fix_single_node_connections_across_islands, &
                              create_new_ect,flag_single_element_islands, & 
                              add_elements_to_ect,remove_unconnected_nodes

      IMPLICIT NONE

      CHARACTER(100) :: filename
      CHARACTER(100) :: grid_name
      CHARACTER(3) :: frmt 
      INTEGER :: i
      INTEGER :: nn,nn_new
      INTEGER :: ne,ne_new
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xy,xy_new
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect,ect_new
      INTEGER :: nbeds
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: bednds
      INTEGER :: nbed 
      INTEGER, DIMENSION(:), ALLOCATABLE :: bedn
      INTEGER :: nbnd 
      INTEGER :: nbou
      INTEGER :: nvel
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbseg
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbnds
      INTEGER :: nope
      INTEGER :: neta
      INTEGER, DIMENSION(:), ALLOCATABLE ::obseg
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: obnds
      INTEGER, DIMENSION(:), ALLOCATABLE :: bndn
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_type
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth
      INTEGER, DIMENSION(:), ALLOCATABLE :: nepn
      INTEGER :: mnepn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: epn
      INTEGER :: ned
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ged2el
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ged2nn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ged2led
      INTEGER :: nied
      INTEGER, DIMENSION(:), ALLOCATABLE :: iedn
      INTEGER, DIMENSION(:), ALLOCATABLE :: ed_type
      INTEGER, DIMENSION(:), ALLOCATABLE :: recv_edge
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el2ged
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el2el
      INTEGER :: nfill_elements
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fill_elements
      INTEGER, DIMENSION(:), ALLOCATABLE :: keep_element
      
      
      ! Read in .msh file
      filename = './Global_50km.14'
      CALL read_header(0,filename,grid_name,ne,nn)
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)

      ALLOCATE(keep_element(ne))

      DO i = 1,5
        ! Compute element connectivity tables
        PRINT*, "Comuting element/edge connectivity"
        CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn) 
        CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)
        CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge,nbed,bedn)
        CALL find_element_edges(ne,ned,ged2el,ged2led,el2ged)
        CALL find_neighbor_elements(ne,ned,ged2el,ged2led,el2el)
        CALL find_boundary_nodes(nn,nbed,bedn,ged2nn,nbnd,bndn)

        ! Remove elements to clean mesh
        keep_element = 1
        PRINT*, "Fixing bad element connectivity"
        CALL flag_problem_elements(ne,ect,nepn,el2el,keep_element)
        CALL flag_isolated_element_patches(ne,el2el,keep_element) 
        CALL fix_single_node_connections_across_islands(nn,nbnd,bndn,nbed,bedn,nepn,epn,ged2nn,keep_element)
        CALL flag_single_element_passages(ne,nn,ect,nbnd,bndn,keep_element)
        CALL create_new_ect(ne,ect,keep_element,ne_new,ect_new)
        
        ne = ne_new
        ect = ect_new 
      ENDDO

      ! Fill in any single element islands
      ALLOCATE(fill_elements(3,ne))
      PRINT*, "Filling single element islands"
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn) 
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge,nbed,bedn)
      CALL flag_single_element_islands(ned,nbed,bedn,ged2nn,xy,nfill_elements,fill_elements)
      CALL add_elements_to_ect(ne,ect,nfill_elements,fill_elements)

      ! To Do: remove X elements: interior node shares four elements,
      ! combine into two elements

      DEALLOCATE(el_type)
      ALLOCATE(el_type(ne))
      el_type = 1

      ! Remove unconnected nodes
      PRINT*, "Removing unconnected nodes"
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn) 
      CALL remove_unconnected_nodes(nn,nepn,xy,ne,ect,nn_new,xy_new,ect_new)

      ! Write new .msh file    
      PRINT*, "Writing .msh file"
      nbeds = 0
      CALL write_header('./clean.14',grid_name,ne,nn)
      CALL write_coords(nn,xy,depth)
      CALL write_connectivity(ne,ect,el_type,nverts)
      nope = 0
      neta = 0
      CALL write_open_boundaries(nope,neta,obseg,obnds)
      nbou = 0
      nvel = 0
      CALL write_flow_boundaries(nbou,nvel,fbseg,fbnds)

      END PROGRAM clean_mesh
