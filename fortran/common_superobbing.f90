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


END MODULE cs
