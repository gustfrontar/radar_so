MODULE cs
!=======================================================================
!
! [PURPOSE:] Thinning of radar data
!
! [HISTORY:] This version produce a superobing of the observations but
! the data is stored in azimuth , range , elevation. Conversion to the 
! lat , lon , z is performed by the observation operator.
!
!=======================================================================
!$USE OMP_LIB

!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)


CONTAINS

SUBROUTINE write_radar(nlon,nlat,nlev,nvar,data_in,ndata_in,grid_az,grid_el,grid_ra,error,ido,lambdar,  &
                       filename,radar_lon,radar_lat,radar_z)
IMPLICIT NONE
INTEGER, INTENT(IN)     :: nlon ,nlat,nlev,nvar
REAL(r_size),INTENT(IN) :: data_in( nlon,nlat,nlev,nvar ) 
REAL(r_size),INTENT(IN) :: grid_az( nlon,nlat,nlev,nvar )
REAL(r_size),INTENT(IN) :: grid_el( nlon,nlat,nlev,nvar )
REAL(r_size),INTENT(IN) :: grid_ra( nlon,nlat,nlev,nvar )
REAL(r_size),INTENT(IN) :: error( nvar ) 
INTEGER     ,INTENT(IN) :: ido(nvar) 
REAL(r_size),INTENT(IN) :: lambdar
INTEGER     ,INTENT(IN) :: ndata_in(  nlon,nlat,nlev,nvar)
CHARACTER(*),INTENT(IN) :: filename
REAL(r_size),INTENT(IN) :: radar_lon , radar_lat , radar_z

REAL(r_size)            :: max_obs(nvar) , min_obs(nvar)

INTEGER  :: ii,jj,kk,iv,nobs
REAL(r_sngl) :: wk(7)

max_obs(nvar) = -999.0d10
min_obs(nvar) = 999.0d10

nobs=0

   DO iv=1,nvar
    DO ii=1,nlon
     DO jj=1,nlat
      DO kk=1,nlev

       !correspond to the location where the stronger echoes are located.
       IF( ndata_in(ii,jj,kk,iv) > 0 )THEN
           wk(2)=REAL(grid_az(ii,jj,kk,iv),r_sngl)
           wk(3)=REAL(grid_el(ii,jj,kk,iv),r_sngl)
           wk(4)=REAL(grid_ra(ii,jj,kk,iv),r_sngl)
           wk(1)=REAL(ido(iv),r_sngl)
           wk(6)=REAL(error(iv),r_sngl)
           wk(5)=REAL(data_in(ii,jj,kk,iv),r_sngl)
           wk(7)=REAL(lambdar,r_sngl)
           WRITE(99)wk
           nobs = nobs + 1
           if( data_in(ii,jj,kk,iv) > max_obs(iv) )max_obs(iv)=data_in(ii,jj,kk,iv)
           if( data_in(ii,jj,kk,iv) < min_obs(iv) )min_obs(iv)=data_in(ii,jj,kk,iv)
           !WRITE(*,*)wk
       ENDIF
      ENDDO
     ENDDO
    ENDDO
    WRITE(*,*)'Max obs :',max_obs(iv),' Min obs:',min_obs(iv),' var = ',ido(iv)
   ENDDO
   WRITE(*,*)'Total obs :',nobs

END SUBROUTINE write_radar


!2D interpolation using box average. Destination grid is assumed to be regular.
SUBROUTINE com_interp_boxavereg(xini,dx,nx,yini,dy,ny,zini,dz,nz,nvar,xin,yin,zin,datain,nin    &
               &                ,data_ave,data_max,data_min,data_std,data_n,data_w,undef,weigth,weigth_ind,is_angle)
  IMPLICIT NONE
  INTEGER , INTENT(IN)           :: nx , ny , nz , nvar , nin
  REAL(r_sngl),INTENT(IN)        :: dx , dy , dz , xini , yini , zini
  REAL(r_sngl),INTENT(IN)        :: undef
  REAL(r_sngl),INTENT(IN)        :: xin(nin),yin(nin),zin(nin),datain(nin,nvar)
  LOGICAL     ,INTENT(IN)        :: weigth(nvar)
  LOGICAL     ,INTENT(IN)        :: is_angle(nvar)
  INTEGER     ,INTENT(IN)        :: weigth_ind(nvar)
  REAL(r_sngl),INTENT(OUT)       :: data_ave(nx,ny,nz,nvar)
  REAL(r_sngl),INTENT(OUT)       :: data_max(nx,ny,nz,nvar)
  REAL(r_sngl),INTENT(OUT)       :: data_min(nx,ny,nz,nvar)
  REAL(r_sngl),INTENT(OUT)       :: data_std(nx,ny,nz,nvar)
  REAL(r_sngl),INTENT(OUT)       :: data_w(nx,ny,nz,nvar)
  INTEGER,INTENT(OUT)            :: data_n(nx,ny,nz,nvar)
  REAL(r_sngl)                   :: w(nin)
  REAL(r_sngl)                   :: tmp_data

  INTEGER                        :: ii , ix , iy , iz , iv


data_n=0

!private(tmp_data,w,ii,ix,iy,iz,iv)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE( tmp_data,w,ii,ix,iy,iz,iv )

DO iv = 1,nvar !Loop over the variables (we can perform OMP over this loop)

  IF ( weigth(iv) )THEN
     w = datain(:,weigth_ind(iv))
  ELSE
     w = 1.0e0
  ENDIF 

  DO ii = 1,nin  !Loop over the input data 

    !Compute the location of the current point with in grid coordinates (rx,ry)
    ix = int( ( xin(ii) - xini ) / dx ) + 1
    iy = int( ( yin(ii) - yini ) / dy ) + 1
    iz = int( ( zin(ii) - zini ) / dz ) + 1

    !Check is the data is within the grid.
    IF( ix <= nx .and. ix >= 1 .and. iy <= ny .and. iy >= 1 .and.   &
        iz <= nz .and. iz >= 1 .and. datain(ii,iv) /= undef .and.   & 
        w(ii) /= undef )THEN

      IF(  data_n(ix,iy,iz,iv) == 0 )THEN
        data_max(ix,iy,iz,iv) = datain(ii,iv)
        data_min(ix,iy,iz,iv) = datain(ii,iv)
        data_ave(ix,iy,iz,iv) = datain(ii,iv) * w(ii)
        data_std(ix,iy,iz,iv) = ( datain(ii,iv) ** 2 )*w(ii)
        data_w  (ix,iy,iz,iv) = w(ii)
        data_n  (ix,iy,iz,iv) = 1

      ELSE
        data_w(ix,iy,iz,iv) = data_w(ix,iy,iz,iv) + w(ii)
        data_n(ix,iy,iz,iv) = data_n(ix,iy,iz,iv) + 1
        IF ( .not. is_angle(iv) )THEN
            data_ave(ix,iy,iz,iv) = data_ave(ix,iy,iz,iv) + datain(ii,iv) * w(ii)
            data_std(ix,iy,iz,iv) = data_std(ix,iy,iz,iv) + ( datain(ii,iv) ** 2 )*w(ii)
        ELSE
            CALL min_angle_distance( data_ave(ix,iy,iz,iv)/data_n(ix,iy,iz,iv) , datain(ii,iv) , tmp_data )
            data_ave(ix,iy,iz,iv) = data_ave(ix,iy,iz,iv) + datain(ii,iv) * w(ii)
            data_std(ix,iy,iz,iv) = data_std(ix,iy,iz,iv) + ( datain(ii,iv) ** 2 )*w(ii)
        ENDIF


        IF( datain(ii,iv) > data_max(ix,iy,iz,iv) )THEN
          data_max(ix,iy,iz,iv) = datain(ii,iv)
        ENDIF
        IF( datain(ii,iv) < data_min(ix,iy,iz,iv) )THEN
          data_min(ix,iy,iz,iv) = datain(ii,iv)
        ENDIF

      ENDIF


    ENDIF

  ENDDO

  WHERE( data_n(:,:,:,iv) > 0)
       data_ave(:,:,:,iv) = data_ave(:,:,:,iv) / REAL( data_n(:,:,:,iv) , r_sngl )
       data_std(:,:,:,iv) = SQRT( data_std(:,:,:,iv)/REAL( data_n(:,:,:,iv) , r_sngl ) - data_ave(:,:,:,iv) ** 2 )
  ENDWHERE

  !If this is an angle check that the value is between 0-360.
  IF( is_angle(iv) )THEN
    WHERE( data_ave(:,:,:,iv) > 360.0e0 )
       data_ave(:,:,:,iv) = data_ave(:,:,:,iv) - 360.0e0
    ENDWHERE
    WHERE( data_ave(:,:,:,iv) < 0.0e0 )
       data_ave(:,:,:,iv) = data_ave(:,:,:,iv) + 360.0e0
    ENDWHERE
  ENDIF


ENDDO

!$OMP END PARALLEL DO

END SUBROUTINE com_interp_boxavereg

SUBROUTINE min_angle_distance( ref_val , val , corrected_val )
IMPLICIT NONE
REAL(r_sngl),INTENT(IN)  :: ref_val , val 
REAL(r_sngl),INTENT(OUT) :: corrected_val
REAL(r_sngl)             :: diff
!Given two angles ref_val and val, find an angle equivalent to val that minimizes the distance between ref_val and val.
!ref_val and val are [0-360], but corrected_val [-180,540]

 diff = ref_val - val 
 IF( diff > 180.0e0 )THEN
   corrected_val = val + 360.0e0
 ELSEIF( diff < -180.0e0 )THEN
   corrected_val = val -360.0e0
 ELSE
   corrected_val = val
 ENDIF

END SUBROUTINE min_angle_distance


END MODULE cs
