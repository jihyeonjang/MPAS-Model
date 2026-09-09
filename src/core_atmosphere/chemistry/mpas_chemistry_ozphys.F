!=================================================================================================================
 module mpas_chemistry_ozphys
 use mpas_kind_types
 use mpas_log,only          : mpas_log_write
 use mpas_io_units,only     : mpas_new_unit,mpas_release_unit
 use mpas_derived_types,only: mpas_pool_type,MPAS_LOG_ERR
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_dimension,mpas_pool_get_array
 use mpas_atmphys_constants,only: degrad
 use mpas_atmphys_utilities,only: physics_error_fatal

!from the CCPP physics repository (see src/core_atmosphere/Externals.cfg):
 use machine,only      : kind_phys
 use module_ozphys,only: ty_ozphys


 implicit none
 private
 public:: ozphys,init_ozphys,update_ozphys


!MPAS interface to the ozone photochemistry parameterization of the CCPP physics library
!(physics/photochem/module_ozphys.F90, McCormack et al. 2006).
!
!The scheme itself is used exactly as it is distributed with CCPP. This module only does what the
!CCPP host model (fv3atm) does in GFS_phys_time_vary.fv3.F90, using MPAS data structures:
!
!  init_ozphys   -> ty_ozphys%load_o3prog    reads global_o3prdlos.f77 (once, at initialization)
!                -> ty_ozphys%setup_o3prog   computes the data-grid to MPAS-mesh interpolation indices
!  update_ozphys -> ty_ozphys%update_o3prog  interpolates the forcing coefficients in time
!
!and the run() procedure is called from mpas_chemistry_driver.
!
!Unlike the earlier MPAS implementation, the monthly-mean coefficients are not interpolated to the MPAS
!mesh and stored for all twelve months: the raw data stay inside the ty_ozphys object (a few MB per MPI
!task), and only the coefficients valid at the current time are kept in the chem_input pool.
!
!All fields exchanged with the scheme are declared with the CCPP kind_phys, and converted explicitly to
!and from RKIND, so that the interface is also correct in single-precision MPAS builds.


!ty_ozphys object holding the ozone forcing data and the scheme itself. It is filled once in
!init_ozphys, and is shared by all blocks owned by an MPI task:
 type(ty_ozphys),save:: ozphys

 logical,save:: ozphys_loaded = .false.


 contains


!=================================================================================================================
 subroutine init_ozphys(configs,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input

!local pointers:
 character(len=StrKIND),pointer:: o3prdlos_file
 character(len=StrKIND),pointer:: o3prognostic_scheme

 integer,pointer:: nCells,nMonths,nO3Coeffs,nO3Levels
 integer,dimension(:),pointer:: O3jindx1,O3jindx2

 real(kind=RKIND),dimension(:),pointer:: latCell
 real(kind=RKIND),dimension(:),pointer:: O3ddy,O3pres

!local variables and arrays:
 character(len=128):: err_message

 integer:: iCell,k
 integer:: o3unit

 real(kind=kind_phys),dimension(:),allocatable:: lat_d,ddy_d

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_o3prdlos_file'      ,o3prdlos_file      )
 call mpas_pool_get_config(configs,'config_o3prognostic_scheme',o3prognostic_scheme)

!--- the run procedure is selected in mpas_chemistry_driver, but the namelist option is checked here so
!    that an invalid value stops the model at initialization:
 if(trim(o3prognostic_scheme) /= 'oz_phys_2015' .and. &
    trim(o3prognostic_scheme) /= 'oz_phys_2006') then
    call mpas_log_write('--- config_o3prognostic_scheme = '//trim(o3prognostic_scheme), &
                        messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_ozphys: config_o3prognostic_scheme must be one of ' // &
                             '''oz_phys_2015'' or ''oz_phys_2006''.')
 endif

 call mpas_pool_get_dimension(mesh,'nCells'    ,nCells   )
 call mpas_pool_get_dimension(mesh,'nMonths'   ,nMonths  )
 call mpas_pool_get_dimension(mesh,'nOPPLevels',nO3Levels)
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs',nO3Coeffs)

 call mpas_pool_get_array(mesh,'latCell',latCell)

 call mpas_pool_get_array(chem_input,'OPPpres'  ,O3pres  )
 call mpas_pool_get_array(chem_input,'OPPjindx1',O3jindx1)
 call mpas_pool_get_array(chem_input,'OPPjindx2',O3jindx2)
 call mpas_pool_get_array(chem_input,'OPPddy'   ,O3ddy   )


!--- read the monthly-mean ozone production and loss coefficients. The data are read by every MPI task,
!    as they are in the CCPP host model:
 if(.not. ozphys_loaded) then
    call mpas_new_unit(o3unit,unformatted=.true.)
    if(o3unit < 0) &
       call physics_error_fatal('subroutine init_ozphys: all file units are taken. Change maxUnits in ' // &
                                'mpas_io_units.F')

    call mpas_log_write(' ')
    call mpas_log_write('--- enter subroutine init_ozphys:')
    call mpas_log_write('    reading ozone photochemistry data from '//trim(o3prdlos_file))

    err_message = ozphys%load_o3prog(trim(o3prdlos_file),o3unit)
    call mpas_release_unit(o3unit)

    if(len_trim(err_message) > 0) then
       call mpas_log_write('--- error while reading '//trim(o3prdlos_file)//':',messageType=MPAS_LOG_ERR)
       call mpas_log_write('--- '//trim(err_message),messageType=MPAS_LOG_ERR)
       call physics_error_fatal('subroutine init_ozphys: ty_ozphys%load_o3prog failed.')
    endif

    call mpas_log_write('    ncf   (number of coefficients) = $i',intArgs=(/ozphys%ncf/)  )
    call mpas_log_write('    nlat  (number of latitudes)    = $i',intArgs=(/ozphys%nlat/) )
    call mpas_log_write('    nlev  (number of levels)       = $i',intArgs=(/ozphys%nlev/) )
    call mpas_log_write('    ntime (number of times)        = $i',intArgs=(/ozphys%ntime/))

    ozphys_loaded = .true.
 endif


!--- the dimensions of the input data file must match those defined in Registry_chems.xml:
 if(ozphys%nlev /= nO3Levels) then
    call mpas_log_write('--- nOPPLevels = $i whereas nlev = $i in the input data file', &
                        intArgs=(/nO3Levels,ozphys%nlev/),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_ozphys: update the dimension nOPPLevels in Registry_chems.xml.')
 elseif(ozphys%ncf /= nO3Coeffs) then
    call mpas_log_write('--- nOPPCoeffs = $i whereas ncf = $i in the input data file', &
                        intArgs=(/nO3Coeffs,ozphys%ncf/),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_ozphys: update the dimension nOPPCoeffs in Registry_chems.xml.')
 elseif(ozphys%ntime /= nMonths) then
    call mpas_log_write('--- nMonths = $i whereas ntime = $i in the input data file', &
                        intArgs=(/nMonths,ozphys%ntime/),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_ozphys: the input data file is not monthly-mean data.')
 endif


!--- save the natural log of the pressure levels of the input data, for output only. The scheme uses its
!    own copy (ty_ozphys%po3):
 do k = 1,nO3Levels
    O3pres(k) = ozphys%po3(k)
 enddo


!--- computation of the indices and weights used to interpolate the input data, given on a uniform gaussian
!    grid, to the MPAS mesh. ty_ozphys%setup_o3prog expects latitudes in degrees:
 allocate(lat_d(nCells))
 allocate(ddy_d(nCells))

 do iCell = 1,nCells
    lat_d(iCell) = real(latCell(iCell)/degrad,kind_phys)
 enddo

 call ozphys%setup_o3prog(lat_d,O3jindx1,O3jindx2,ddy_d)

 do iCell = 1,nCells
    O3ddy(iCell) = real(ddy_d(iCell),RKIND)
 enddo

 deallocate(lat_d)
 deallocate(ddy_d)

! call mpas_log_write('--- end subroutine init_ozphys:')

 end subroutine init_ozphys

!=================================================================================================================
 subroutine update_ozphys(julian,curr_julday,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 integer,intent(in):: julian
 real(kind=RKIND),intent(in):: curr_julday

!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input

!local pointers:
 integer,pointer:: nCellsSolve,nO3Coeffs,nO3Levels
 integer,dimension(:),pointer:: O3jindx1,O3jindx2

 real(kind=RKIND),dimension(:),pointer    :: O3ddy
 real(kind=RKIND),dimension(:,:,:),pointer:: O3coeffs

!local variables and arrays:
 integer,parameter:: chunk_size = 256

 integer:: iCell,iCell1,iCell2,k,n,nChunk,n1,n2

 real(kind=kind_phys):: rjday
 real(kind=kind_phys),dimension(:),allocatable    :: ddy_d
 real(kind=kind_phys),dimension(:,:,:),allocatable:: ozpl

!-----------------------------------------------------------------------------------------------------------------
! call mpas_log_write(' ')
! call mpas_log_write('--- enter subroutine update_ozphys:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nOPPLevels' ,nO3Levels  )
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs' ,nO3Coeffs  )

 call mpas_pool_get_array(chem_input,'OPPjindx1',O3jindx1)
 call mpas_pool_get_array(chem_input,'OPPjindx2',O3jindx2)
 call mpas_pool_get_array(chem_input,'OPPddy'   ,O3ddy   )
 call mpas_pool_get_array(chem_input,'OPPcoeffs',O3coeffs)


!--- current time of the year, in days, counted as in the CCPP host model: rjday is the day of the year
!    (1.0 at 00Z on 1 January) plus the fraction of the current day. curr_julday counts from 0.0 at 00Z
!    on 1 January, hence the offset by one day:
 rjday = real(curr_julday,kind_phys) + 1._kind_phys
 if(rjday < ozphys%time(1)) rjday = rjday + 365._kind_phys

 call find_photochem_time_index(ozphys%ntime,ozphys%time,rjday,n1,n2)

 call mpas_log_write('--- time interpolation of the ozone photochemistry coefficients:')
 call mpas_log_write('    julian day = $i',intArgs=(/julian/))
 call mpas_log_write('    rjday      = $r',realArgs=(/real(rjday,RKIND)/))
 call mpas_log_write('    n1         = $i',intArgs=(/n1/))
 call mpas_log_write('    n2         = $i',intArgs=(/n2/))


!--- time interpolation of the coefficients, done in chunks of cells to limit the size of the temporary
!    array needed to reorder the output of ty_ozphys%update_o3prog, which is dimensioned (cells,levels,
!    coefficients), into the MPAS array OPPcoeffs, which is dimensioned (coefficients,levels,cells):
 allocate(ozpl(chunk_size,nO3Levels,nO3Coeffs))
 allocate(ddy_d(chunk_size))

 do iCell1 = 1,nCellsSolve,chunk_size
    iCell2 = min(iCell1+chunk_size-1,nCellsSolve)
    nChunk = iCell2-iCell1+1

    do iCell = 1,nChunk
       ddy_d(iCell) = real(O3ddy(iCell1+iCell-1),kind_phys)
    enddo

    call ozphys%update_o3prog(O3jindx1(iCell1:iCell2),O3jindx2(iCell1:iCell2),ddy_d(1:nChunk), &
                              rjday,n1,n2,ozpl(1:nChunk,:,:))

    do iCell = 1,nChunk
       do k = 1,nO3Levels
          do n = 1,nO3Coeffs
             O3coeffs(n,k,iCell1+iCell-1) = real(ozpl(iCell,k,n),RKIND)
          enddo
       enddo
    enddo
 enddo

 deallocate(ozpl)
 deallocate(ddy_d)

! call mpas_log_write('--- end subroutine update_ozphys:')

 end subroutine update_ozphys

!=================================================================================================================
 subroutine find_photochem_time_index(ntime,time,rjday,n1,n2)
!=================================================================================================================

!This is the routine of the same name contained in GFS_phys_time_vary.fv3.F90 of the CCPP physics library,
!repeated here with the correct kind for the array of times, so that the coefficients are interpolated in
!time exactly as they are in the GFS.

!input arguments:
 integer,intent(in):: ntime
 real(kind=kind_phys),intent(in),dimension(ntime+1):: time
 real(kind=kind_phys),intent(in):: rjday

!output arguments:
 integer,intent(out):: n1,n2

!local variables:
 integer:: j

!-----------------------------------------------------------------------------------------------------------------

 n2 = ntime + 1
 do j = 2,ntime
    if(rjday < time(j)) then
       n2 = j
       exit
    endif
 enddo
 n1 = n2 - 1
 if(n2 > ntime) n2 = n2 - ntime

 end subroutine find_photochem_time_index

!=================================================================================================================
 end module mpas_chemistry_ozphys
!=================================================================================================================
