!=================================================================================================================
 module mpas_chemistry_opp
 use mpas_kind_types
 use mpas_log,only          : mpas_log_write
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_dimension,mpas_pool_get_array
 use mpas_atmphys_constants,only: degrad
 use mpas_atmphys_utilities,only: physics_error_fatal

!from wrf physics:
 use module_ra_cam_support,only: r8,getfactors

 use o3_gfs,only: latsozp,levozp,oz_coeff,timeoz,oz_lat,oz_pres,ozplin,read_o3data


 implicit none
 private
 public:: init_opp,opp_from_MPAS


!initialization and interpolation to MPAS mesh of input data needed to run the CHEM2D Ozone Photochemistry
!Parameterization (CHEM2D-OPP) of McCormack et al. (2006).


 contains


!=================================================================================================================
 subroutine init_opp(configs,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh


!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input


!local variables:
 logical,pointer:: o3climatology,o3prognostic
 integer:: i,iCell,i1,i2,k,m,n
 integer,pointer:: nCells,nMonths,nO3Coeffs,nO3Levels

 real(kind=RKIND):: dlat,dlatCell,lat,lon
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: O3pres
 real(kind=RKIND),dimension(:,:,:,:),pointer:: O3data

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_opp:')

 call mpas_pool_get_config(configs,'config_o3climatology',o3climatology)
 call mpas_pool_get_config(configs,'config_o3prognostic' ,o3prognostic )

 call mpas_pool_get_dimension(mesh,'nCells'    ,nCells   )
 call mpas_pool_get_dimension(mesh,'nMonths'   ,nMonths  )
 call mpas_pool_get_dimension(mesh,'nOPPLevels',nO3Levels)
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs',nO3Coeffs)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)

 call mpas_pool_get_array(chem_input,'OPPpres',O3pres)
 call mpas_pool_get_array(chem_input,'OPPdata',O3data)


!--- read monthly-mean CHEM2-OPP linearized coefficients:
!if(o3prognostic) then
    call read_o3data()
!elseif(o3climatology) then
!   return
!endif


!--- make sure that the dimensions in the input O3CHem2D-OPP file are the same as the ones set in Registry.xml:
 if(nO3Levels /= levozp) then
    call mpas_log_write('--- nO3Levels = $i whereas levozp = $i',intArgs=(/nO3Levels,levozp/))
    call physics_error_fatal('--- nO3Levels is different than levozp from input data file:')
 elseif(nO3Coeffs /= oz_coeff) then
    call mpas_log_write('--- nO3Coeffs = $i whereas oz_coeff = $i',intArgs=(/nO3Coeffs,oz_coeff/))
 elseif(nMonths /= timeoz) then
    call mpas_log_write('--- nMonths = $i whereas timeiz = $i',intArgs=(/nMonths,timeoz/))
 endif


!--- interpolation of input CHEM2D-OPP coefficients to the MPAS grid: this loop is the same as in subroutine
!    init_o3climatology (see module mpas_atmphys_o3climatology, lines 149-184):
 do k = 1,nO3Levels
    O3pres(k) = oz_pres(k)
 enddo

 do iCell = 1,nCells
    lat = latCell(iCell)/degrad
    lon = lonCell(iCell)/degrad
    if(lat .gt. oz_lat(latsozp)) then
     i1 = latsozp
     i2 = latsozp
    elseif(lat .lt. oz_lat(1)) then
       i1 = 1
       i2 = 1
    else
       do i = 1, latsozp
          if(lat.ge.oz_lat(i) .and. lat.lt.oz_lat(i+1)) exit
       enddo
       i1 = i
       i2 = i+1
    endif

    do m = 1,nMonths
       do n = 1,nO3Coeffs
          do k = 1,nO3Levels
             dlat     = oz_lat(i2)-oz_lat(i1)
             dlatCell = lat-oz_lat(i1)

             if(dlat == 0.) then
                O3data(m,n,k,iCell) = ozplin(i1,k,n,m)
             else
                O3data(m,n,k,iCell) = ozplin(i1,k,n,m) &
                     + (ozplin(i2,k,n,m)-ozplin(i1,k,n,m))*dlatCell/dlat
             endif
          enddo
       enddo
    enddo
 enddo

!do m = 1,nMonths
!   do iCell = 1,nCells
!      do k = 1,nO3Levels
!         if(m.eq.1 .or. m.eq.nMonths) &
!            call mpas_log_write('$i $i $i $r $r $r $r $r $r $r',intArgs=(/m,iCell,k/),realArgs=(/O3pres(k), &
!                               (O3data(m,n,k,iCell),n=1,nO3Coeffs)/))
!      enddo
!      call mpas_log_write(' ')
!   enddo
!   call mpas_log_write(' ')
!enddo

 call mpas_log_write('--- end subroutine init_opp:')

 end subroutine init_opp

!=================================================================================================================
 subroutine opp_from_MPAS(julian,curr_julday,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 integer,intent(in):: julian
 real(kind=RKIND):: curr_julday

!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input

!local pointers:
 integer,pointer:: nCellsSolve,nO3Coeffs,nO3Levels,nMonths
 integer:: iCell,imonth,k,n

 real(kind=RKIND),dimension(:,:,:),pointer  :: O3coeffs
 real(kind=RKIND),dimension(:,:,:,:),pointer:: O3data

!local variables:
 logical:: finddate
 logical:: ozncyc

 integer:: ijul,m,nm,np,np1
 integer, dimension(12) ::  date_oz
 data date_oz/16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350/

 real(kind=r8):: intjulian_r8
 real(kind=r8):: cdayozp_r8, cdayozm_r8
 real(kind=r8):: fact1_r8, fact2_r8
 real(kind=RKIND):: fact1,fact2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine opp_from_MPAS:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nMonths'    ,nMonths    )
 call mpas_pool_get_dimension(mesh,'nOPPLevels' ,nO3Levels  )
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs' ,nO3Coeffs  )

 call mpas_pool_get_array(chem_input,'OPPdata'  ,O3data  )
 call mpas_pool_get_array(chem_input,'OPPcoeffs',O3coeffs)

!--- for now, fix the array GFSO3coeffs to GFSO3data for July:
!imonth = 7
!O3coeffs(:,:,:) = O3data(imonth,:,:,:)

 ozncyc = .true.

!julian starts from 0.0 at 0Z on 1 Jan.
 intjulian_r8 = real(julian + 1.0_RKIND, r8)    ! offset by one day

!jan 1st 00z is julian=1.0 here
 ijul=int(intjulian_r8)
!note that following will drift. need to use actual month/day info to compute julian.
 intjulian_r8 = intjulian_r8 - real(ijul, r8)
 ijul = mod(ijul,365)
 if(ijul .eq. 0) ijul=365
 intjulian_r8 = intjulian_r8 + real(ijul, r8)
 np1=1

 finddate=.false.
 do m = 1, nMonths
    if(date_oz(m).gt.intjulian_r8 .and. .not.finddate) then
       np1 = m
       finddate = .true.
    endif
 enddo
 cdayozp_r8=date_oz(np1)

 if(np1 .gt. 1) then
    cdayozm_r8 = date_oz(np1-1)
    np = np1
    nm = np-1
 else
    cdayozm_r8 = date_oz(12)
    np = np1
    nm = nMonths
 endif
 call getfactors(ozncyc,np1,cdayozm_r8,cdayozp_r8,intjulian_r8,fact1_r8,fact2_r8)
 fact1 = real(fact1_r8, RKIND)
 fact2 = real(fact2_r8, RKIND)

 call mpas_log_write('--- end subroutine getfactors:')
 call mpas_log_write('nm    = $i',intArgs =(/nm/))
 call mpas_log_write('np    = $i',intArgs =(/np/))
 call mpas_log_write('fact1 = $r',realArgs=(/fact1/))
 call mpas_log_write('fact2 = $r',realArgs=(/fact2/))

!time interpolation:
 do iCell = 1,nCellsSolve
    do k = 1,nO3Levels
       do n = 1,nO3Coeffs
          O3Coeffs(n,k,iCell) = fact1*O3data(nm,n,k,iCell) + fact2*O3data(np,n,k,iCell)
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine opp_from_MPAS:')

 end subroutine opp_from_MPAS

!=================================================================================================================
 end module mpas_chemistry_opp
!=================================================================================================================
