      PROGRAM rotate_test

      USE rotate_mod

      IMPLICIT NONE

      REAL(rp) :: LON_POLE,LAT_POLE 
      REAL(rp) :: RADPHI, RADTHETA
      INTEGER :: NP
      REAL(rp), DIMENSION(:), ALLOCATABLE :: LON,LAT
      REAL(rp), DIMENSION(:), ALLOCATABLE :: LONR,LATR
      REAL(rp), DIMENSION(:), ALLOCATABLE :: LONR2,LATR2
      REAL(rp), DIMENSION(:), ALLOCATABLE :: UVEC,VVEC 
      REAL(rp), DIMENSION(:), ALLOCATABLE :: UVECR,VVECR 
      REAL(rp), DIMENSION(:), ALLOCATABLE :: UVECUR,VVECUR 
      REAL(rp), DIMENSION(:), ALLOCATABLE :: UVECR2,VVECR2 
      REAL(rp), DIMENSION(:), ALLOCATABLE :: UVECUR2,VVECUR2 
      REAL(rp), DIMENSION(:), ALLOCATABLE :: ANGLED
      INTEGER :: ncid,lon_dimid,lat_dimid,t_dimid
      INTEGER :: lon_varid,lat_varid
      INTEGER :: u_varid,v_varid
      INTEGER :: nlon,nlat,ntime
      REAL(rp), DIMENSION(:), ALLOCATABLE :: lon_nc,lat_nc
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: u_nc,v_nc
      INTEGER :: i,j,n
      REAL (rp):: ROTMAT(3,3) = reshape( (/ 1.0, 0.0, 0.0, &
                                  0.0, 1.0, 0.0, &
                                  0.0, 0.0, 1.0 /), (/3, 3/)) ;
      REAL (rp), ALLOCATABLE, DIMENSION(:,:,:) :: RVELF, RVELI
      INTEGER :: IERR
      CHARACTER (:), ALLOCATABLE :: nc_file

      LON_POLE = -42.8906d0 
      LAT_POLE = 72.3200d0
      nc_file = 'wnd10mx0.5.gdas.200505.grb2.nc'
    

      CALL check(NF90_OPEN(nc_file,NF90_NOWRITE,ncid))      

      CALL check(NF90_INQ_DIMID(ncid,'lon',lon_dimid))
      CALL check(NF90_INQ_DIMID(ncid,'lat',lat_dimid))
      CALL check(NF90_INQ_DIMID(ncid,'time',t_dimid))

      CALL check(NF90_INQUIRE_DIMENSION(ncid,lon_dimid,len=nlon))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,lat_dimid,len=nlat))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,t_dimid,len=ntime))
      PRINT*, nlon,nlat,ntime

      CALL check(NF90_INQ_VARID(ncid,'lon',lon_varid))
      CALL check(NF90_INQ_VARID(ncid,'lat',lat_varid))
      CALL check(NF90_INQ_VARID(ncid,'U_GRD_L103',u_varid))
      CALL check(NF90_INQ_VARID(ncid,'V_GRD_L103',v_varid))

      ALLOCATE(lon_nc(nlon),lat_nc(nlat))
      CALL check(NF90_GET_VAR(ncid,lon_varid,lon_nc))
      CALL check(NF90_GET_VAR(ncid,lat_varid,lat_nc))
      ALLOCATE(u_nc(nlon,nlat,ntime),v_nc(nlon,nlat,ntime))
      CALL check(NF90_GET_VAR(ncid,u_varid,u_nc))
      CALL check(NF90_GET_VAR(ncid,v_varid,v_nc))

      NP = nlon*nlat
      ALLOCATE(LON(NP), LAT(NP))
      ALLOCATE(LONR(NP), LATR(NP))
      ALLOCATE(RVELF(2,2,NP), RVELI(2,2,NP))
      ALLOCATE(LONR2(NP), LATR2(NP))
      ALLOCATE(ANGLED(NP))
      ALLOCATE(UVEC(NP),VVEC(NP))
      ALLOCATE(UVECR(NP),VVECR(NP))
      ALLOCATE(UVECUR(NP),VVECUR(NP))
      ALLOCATE(UVECR2(NP),VVECR2(NP))
      ALLOCATE(UVECUR2(NP),VVECUR2(NP))

      n = 1
      DO i = 1,nlat
        DO j = 1,nlon

          IF (lon_nc(j) > 180d0) THEN
            LON(n) = lon_nc(j) - 360d0
          ELSE
            LON(n) = lon_nc(j)
          ENDIF

          LAT(n) = lat_nc(i)

          UVEC(n) = u_nc(j,i,1)
          VVEC(n) = v_nc(j,i,1)

          n = n+1
        ENDDO
      ENDDO

      LONR = 0d0
      LATR = 0d0
      LONR2 = 0d0
      LATR2 = 0d0
      UVECR = 0d0
      VVECR = 0d0
      UVECUR = 0d0
      VVECUR = 0d0
      RVELF = 0d0
      RVELI = 0d0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ADCIRC
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      RADPHI = LON_POLE*DEG2RAD ;
      RADTHETA = LAT_POLE*DEG2RAD ;

      CALL GET_ROTMAT_ZNVEC( ROTMAT, RADPHI, RADTHETA )   ;

      CALL CHECK_RTOTMAT( ROTMAT, IERR )

      DO i = 1,NP
        LON(i) = LON(i)*DEG2RAD
        LAT(i) = LAT(i)*DEG2RAD
      ENDDO
      CALL SPCOORSROTS1( ROTMAT, LONR, LATR, LON, LAT, NP ) ;

      CALL SPVECROTSMAT( RVELF, RVELI, ROTMAT, &
                          LONR, LATR, LON, LAT, NP ) ;

      CALL MAP2DSPVECS( UVECR, VVECR, UVEC, VVEC, RVELF, NP ) 
      CALL MAP2DSPVECS( UVECUR, VVECUR, UVECR, VVECR, RVELI, NP ) 

      DO i = 1,NP
        LON(i) = LON(i)*RAD2DEG
        LAT(i) = LAT(i)*RAD2DEG
      ENDDO
 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! WAVEWATCH III
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL W3LLTOEQ ( LAT , LON, LATR2, LONR2,     &              
                      ANGLED, LAT_POLE, LON_POLE, NP )             

      UVECR2 = UVEC
      VVECR2 = VVEC
      CALL W3XYRTN ( NP, UVECR2, VVECR2, -ANGLED)
      UVECUR2 = UVECR
      VVECUR2 = VVECR
      CALL W3XYRTN ( NP, UVECUR2, VVECUR2, ANGLED)

      DO i = 1,np
         PRINT("(3(F18.12,F18.12,15x),F18.12)"), LON(i), LAT(i), RAD2DEG*LONR(i), RAD2DEG*LATR(i), LONR2(i), LATR2(i), ANGLED(i)
      !   PRINT("(4(F24.12))"), UVEC(i), VVEC(i), UVECR2(i), VVECR2(i)
      !   PRINT("(4(F24.12))"), UVEC(i), VVEC(i), UVECUR(i), VVECUR(i)
      !   PRINT("(2(F8.3,F8.3,15x),6(F18.12,F18.12,15x))"),LON(i), LAT(i), RAD2DEG*LONR(i), RAD2DEG*LATR(i), UVEC(i), VVEC(i), UVECR(i), VVECR(i), UVECUR(i), VVECUR(i), UVECUR2(i), VVECUR2(i), UVECR2(i), VVECR2(i)
      ENDDO

      END PROGRAM rotate_test      

