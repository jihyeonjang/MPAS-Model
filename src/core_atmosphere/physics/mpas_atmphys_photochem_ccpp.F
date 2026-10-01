!=================================================================================================================
 module mpas_atmphys_photochem_ccpp
!=================================================================================================================
!
! This is the only file in MPAS that imports the CCPP photochemistry modules. Everything that MPAS needs
! from the CCPP is re-exported here, so that no other MPAS source has to be compiled with the CCPP module
! directory on its include path. That matters because the CCPP and NoahMP both define a Fortran module
! named "machine" (physics_noahmp/utility/Machine.F90), and only one of the two may be visible to any
! given compilation.
!
! The CCPP sources are used exactly as they come from the ufs-community/ccpp-physics checkout made by
! manage_externals; nothing in physics_photochem/ccpp-physics is modified. The derived types ty_ozphys and
! ty_h2ophys hold the forcing data (latitude axis, pressure axis, time axis and the production and loss
! coefficients) internally, so MPAS does not have to store them in the Registry.
!
! The two instances below are module variables rather than Registry fields because the forcing data are
! global: they do not depend on the MPAS block decomposition, and every block interpolates from the same
! arrays. They are read once, from subroutine init_ozphys and init_h2ophys.
!
! To add another CCPP photochemistry scheme:
!   1. add its source file to physics_photochem/ccpp_photochem.sparse, to CCPP_OBJS in
!      physics_photochem/Makefile and to ATMOSPHERE_CORE_PHOTOCHEM_CCPP_SOURCES in ../CMakeLists.txt,
!   2. import it here, and
!   3. call it from mpas_atmphys_driver_photochem.F.
! No other file has to be touched.
!
!=================================================================================================================

 use mpas_kind_types
 use mpas_derived_types,only: MPAS_LOG_ERR
 use mpas_log,only: mpas_log_write
 use mpas_io_units,only: mpas_new_unit,mpas_release_unit

 use mpas_atmphys_utilities,only: physics_error_fatal

 use module_ozphys ,only: ty_ozphys
 use module_h2ophys,only: ty_h2ophys

 implicit none
 private

 public:: ozphys,h2ophys
 public:: ozphys_loaded,h2ophys_loaded
 public:: init_ozphys,init_h2ophys
 public:: update_ozphys_forcing,update_h2ophys_forcing
 public:: find_photochem_time_index

!--- the CCPP objects that hold the forcing data:
 type(ty_ozphys) ,save:: ozphys
 type(ty_h2ophys),save:: h2ophys

 logical,save:: ozphys_loaded  = .false.
 logical,save:: h2ophys_loaded = .false.

 logical,save:: log_o3_time_once  = .true.
 logical,save:: log_h2o_time_once = .true.

 contains

!=================================================================================================================
 subroutine init_ozphys(o3_scheme,o3file)
!=================================================================================================================

!Reads the ozone production and loss forcing file into the CCPP ty_ozphys object, and checks that the
!number of coefficients in the file matches what the selected scheme expects.

!input arguments:
 character(len=*),intent(in):: o3_scheme
 character(len=*),intent(in):: o3file

!local variables:
 character(len=128):: errmsg
 logical:: lexist,lopened
 integer:: photochem_fileID

!-----------------------------------------------------------------------------------------------------------------

 if(trim(o3_scheme) == 'off') return
 if(ozphys_loaded) return

 call mpas_log_write('')
 call mpas_log_write('--- enter subroutine init_ozphys:')

 inquire(file=trim(o3file),exist=lexist)
 if(.not.lexist) then
    call mpas_log_write('--- config_o3prdlos_file = '//trim(o3file),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_ozphys: the ozone production and loss file was not found.')
 endif

!--- load_o3prog opens the file with convert='big_endian', reads nlat, nlev, ntime and ncf from its
!    header, and allocates the data arrays accordingly. It returns an empty string on success:
 call mpas_new_unit(photochem_fileID,unformatted=.true.)
 if(photochem_fileID < 0) &
    call physics_error_fatal('subroutine init_ozphys: no free Fortran I/O unit is available.')
 errmsg = ozphys%load_o3prog(trim(o3file),photochem_fileID)
 inquire(unit=photochem_fileID,opened=lopened)
 if(lopened) close(photochem_fileID)
 call mpas_release_unit(photochem_fileID)
 if(len_trim(errmsg) > 0) then
    call mpas_log_write('--- '//trim(errmsg),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_ozphys: failed to read '//trim(o3file))
 endif
 call validate_photochem_time_axis('ozone',trim(o3file),ozphys%ntime,ozphys%time)

!
!    run_o3prog_2015 indexes prod(:,5) and prod(:,6), the climatological temperature and the
!    climatological ozone mixing ratio, unconditionally: a file with fewer coefficients would be read
!    out of bounds.
!
!    run_o3prog_2006 updates ozone only inside "if (this%ncf == 2)" and "if (this%ncf == 4)" blocks.
!    There is no branch for any other value, so with the six-coefficient GFS file it reads the file,
!    interpolates the coefficients in latitude, time and height, and then leaves ozone completely
!    unchanged, without any error. Its ozib temporary is also only assigned inside those two blocks,
!    so the optional do3_dt_ozmx diagnostic would read an undefined value.
!
 select case(trim(o3_scheme))

    case('oz_phys_2015')
       if(ozphys%ncf /= 6) then
          call mpas_log_write('--- number of coefficients in '//trim(o3file)//' = $i', &
                              intArgs=(/ozphys%ncf/),messageType=MPAS_LOG_ERR)
          call physics_error_fatal('subroutine init_ozphys: oz_phys_2015 needs an ozone forcing file ' // &
                                   'with 6 coefficients.')
       endif

    case('oz_phys_2006')
       if(ozphys%ncf /= 2 .and. ozphys%ncf /= 4) then
          call mpas_log_write('--- number of coefficients in '//trim(o3file)//' = $i', &
                              intArgs=(/ozphys%ncf/),messageType=MPAS_LOG_ERR)
          call mpas_log_write('--- run_o3prog_2006 supports only ozone forcing files with 2 or 4 coefficients.', &
                              messageType=MPAS_LOG_ERR)
          call mpas_log_write('--- the current GFS-v17 ozone forcing file has 6 coefficients ' // &
                              'and is intended for config_photochem_o3_scheme = ''oz_phys_2015''.', &
                              messageType=MPAS_LOG_ERR)
          call physics_error_fatal('subroutine init_ozphys: oz_phys_2006 requires an ozone forcing file ' // &
                                   'with 2 or 4 coefficients.')
       endif

 end select

 ozphys_loaded = .true.

 call mpas_log_write('--- ozone photochemistry: read '//trim(o3file))
 call mpas_log_write('        nlat  = $i, nlev = $i',intArgs=(/ozphys%nlat,ozphys%nlev/))
 call mpas_log_write('        ntime = $i, ncf  = $i',intArgs=(/ozphys%ntime,ozphys%ncf/))
 call mpas_log_write('        forcing period = $r days',realArgs=(/ozphys%time(ozphys%ntime+1)-ozphys%time(1)/))

 end subroutine init_ozphys

!=================================================================================================================
 subroutine init_h2ophys(h2o_scheme,h2ofile)
!=================================================================================================================

!Reads the water vapor production and loss forcing file into the CCPP ty_h2ophys object.

!input arguments:
 character(len=*),intent(in):: h2o_scheme
 character(len=*),intent(in):: h2ofile

!local variables:
 character(len=128):: errmsg
 logical:: lexist,lopened
 integer:: photochem_fileID

!-----------------------------------------------------------------------------------------------------------------

 if(trim(h2o_scheme) == 'off') return
 if(h2ophys_loaded) return

 call mpas_log_write('')
 call mpas_log_write('--- enter subroutine init_h2ophys:')

 inquire(file=trim(h2ofile),exist=lexist)
 if(.not.lexist) then
    call mpas_log_write('--- config_h2oprdlos_file = '//trim(h2ofile),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_h2ophys: the water vapor production and loss file was not found.')
 endif

 call mpas_new_unit(photochem_fileID,unformatted=.true.)
 if(photochem_fileID < 0) &
    call physics_error_fatal('subroutine init_h2ophys: no free Fortran I/O unit is available.')
 errmsg = h2ophys%load(trim(h2ofile),photochem_fileID)
 inquire(unit=photochem_fileID,opened=lopened)
 if(lopened) close(photochem_fileID)
 call mpas_release_unit(photochem_fileID)
 if(len_trim(errmsg) > 0) then
    call mpas_log_write('--- '//trim(errmsg),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_h2ophys: failed to read '//trim(h2ofile))
 endif
 call validate_photochem_time_axis('water-vapor',trim(h2ofile),h2ophys%ntime,h2ophys%time)

!--- the run procedure uses pltc(:,1), pltc(:,2) and pltc(:,3):
 if(h2ophys%ncf /= 3) then
    call mpas_log_write('--- number of coefficients in '//trim(h2ofile)//' = $i',intArgs=(/h2ophys%ncf/), &
                        messageType=MPAS_LOG_ERR)
    call physics_error_fatal('subroutine init_h2ophys: the water vapor forcing file must have 3 coefficients.')
 endif

 h2ophys_loaded = .true.

 call mpas_log_write('--- water vapor photochemistry: read '//trim(h2ofile))
 call mpas_log_write('        nlat  = $i, nlev = $i',intArgs=(/h2ophys%nlat,h2ophys%nlev/))
 call mpas_log_write('        ntime = $i, ncf  = $i',intArgs=(/h2ophys%ntime,h2ophys%ncf/))
 call mpas_log_write('        forcing period = $r days',realArgs=(/h2ophys%time(h2ophys%ntime+1)-h2ophys%time(1)/))

 end subroutine init_h2ophys

!=================================================================================================================
 subroutine validate_photochem_time_axis(species,file_name,ntime,time)
!=================================================================================================================

!The CCPP file format stores ntime coefficient records but ntime+1 time coordinates. The extra coordinate
!is the first forcing record wrapped forward by one climatological year. MPAS therefore requires a strictly
!increasing 365-day annual axis; accepting another period would make the Gregorian day-of-year mapping and
!the final-to-first cyclic interpolation ambiguous.

 character(len=*),intent(in):: species,file_name
 integer,intent(in):: ntime
 real(kind=RKIND),intent(in),dimension(ntime+1):: time

 integer:: j
 real(kind=RKIND):: annual_period,tol

 if(ntime < 2) then
    call mpas_log_write('--- '//trim(species)//' photochemistry file = '//trim(file_name),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('validate_photochem_time_axis: at least two forcing times are required.')
 endif

 do j=1,ntime
    if(time(j+1) <= time(j)) then
       call mpas_log_write('--- '//trim(species)//' photochemistry file = '//trim(file_name),messageType=MPAS_LOG_ERR)
       call mpas_log_write('--- non-increasing forcing-time index = $i',intArgs=(/j/),messageType=MPAS_LOG_ERR)
       call physics_error_fatal('validate_photochem_time_axis: forcing times must be strictly increasing.')
    endif
 enddo

 annual_period = time(ntime+1)-time(1)
 tol = 1.e-3_RKIND
 if(abs(annual_period-365._RKIND) > tol) then
    call mpas_log_write('--- '//trim(species)//' photochemistry file = '//trim(file_name),messageType=MPAS_LOG_ERR)
    call mpas_log_write('--- forcing annual period = $r days',realArgs=(/annual_period/),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('validate_photochem_time_axis: expected a 365-day annual climatology.')
 endif

 end subroutine validate_photochem_time_axis

!=================================================================================================================
 subroutine update_ozphys_forcing(idx1,idx2,idxh,rjday,ozpl)
!=================================================================================================================

!Use the untouched CCPP interpolation routine for all ordinary monthly intervals. The only special case is
!the cyclic interval from the final forcing record to the first record of the next year. CCPP stores the
!next-year time coordinate in time(ntime+1), but its data array has only ntime records; passing idxt2=1 to
!the stock update_o3prog routine would combine data record 1 with time(1), rather than time(ntime+1), and
!therefore gives incorrect weights across December/January. Handle that single cyclic interval in the MPAS
!adapter without modifying the imported CCPP source.

 integer,intent(in),dimension(:):: idx1,idx2
 real(kind=RKIND),intent(in),dimension(:):: idxh
 real(kind=RKIND),intent(in):: rjday
 real(kind=RKIND),intent(out),dimension(:,:,:):: ozpl

 integer:: n1,n2,nc,l,j,j1,j2
 real(kind=RKIND):: rjday_work,tem,tx1,tx2,t1,t2,cycle_days,tol
 logical:: cyclic

 call mpas_log_write('--- enter subroutine update_ozphys_forcing:')

 cycle_days = ozphys%time(ozphys%ntime+1)-ozphys%time(1)
 rjday_work = rjday
 if(rjday_work < ozphys%time(1)) rjday_work = rjday_work + cycle_days
 call find_photochem_time_index(ozphys%ntime,ozphys%time,rjday_work,n1,n2,cyclic)

!
!--- Temporary one-time diagnostic of forcing-time interpolation
!
!$OMP CRITICAL(photochem_time_log)
if(log_o3_time_once) then

   if(cyclic) then
      t1 = ozphys%time(ozphys%ntime)
      t2 = ozphys%time(ozphys%ntime+1)
   else
      t1 = ozphys%time(n1)
      t2 = ozphys%time(n2)
   endif

   tx1 = (t2-rjday_work)/(t2-t1)
   tx2 = 1._RKIND-tx1

   call mpas_log_write('--- O3 forcing-time interpolation check:')
   call mpas_log_write('    rjday = $r, rjday_work = $r', &
                       realArgs=(/rjday,rjday_work/))
   call mpas_log_write('    data indices: n1 = $i, n2 = $i', &
                       intArgs=(/n1,n2/))
   call mpas_log_write('    time coordinates: t1 = $r, t2 = $r', &
                       realArgs=(/t1,t2/))
   call mpas_log_write('    weights: tx1 = $r, tx2 = $r', &
                       realArgs=(/tx1,tx2/))

   if(cyclic) then
      call mpas_log_write('    cyclic interval = TRUE')
   else
      call mpas_log_write('    cyclic interval = FALSE')
   endif

   log_o3_time_once = .false.
endif
!$OMP END CRITICAL(photochem_time_log)

 if(.not.cyclic) then
    call ozphys%update_o3prog(idx1,idx2,idxh,rjday_work,n1,n2,ozpl)
    return
 endif

 t1 = ozphys%time(ozphys%ntime)
 t2 = ozphys%time(ozphys%ntime+1)
 tol = 100._RKIND*epsilon(1._RKIND)*max(1._RKIND,abs(t2))
 if(rjday_work < t1-tol .or. rjday_work > t2+tol) &
    call physics_error_fatal('update_ozphys_forcing: model time falls outside the cyclic forcing interval.')
 tx1 = (t2-rjday_work)/(t2-t1)
 tx2 = 1._RKIND-tx1

 do nc=1,ozphys%ncf
 do l=1,ozphys%nlev
 do j=1,size(idx1)
    j1 = idx1(j)
    j2 = idx2(j)
    tem = 1._RKIND-idxh(j)
    ozpl(j,l,nc) = tx1*(tem*ozphys%data(j1,l,nc,ozphys%ntime) + idxh(j)*ozphys%data(j2,l,nc,ozphys%ntime)) &
                  +tx2*(tem*ozphys%data(j1,l,nc,1)            + idxh(j)*ozphys%data(j2,l,nc,1))
 enddo
 enddo
 enddo

 call mpas_log_write('--- end subroutine update_ozphys_forcing:')

 end subroutine update_ozphys_forcing

!=================================================================================================================
 subroutine update_h2ophys_forcing(idx1,idx2,idxh,rjday,h2opl)
!=================================================================================================================

!Same cyclic-time correction as update_ozphys_forcing, applied to the H2O forcing table.

 integer,intent(in),dimension(:):: idx1,idx2
 real(kind=RKIND),intent(in),dimension(:):: idxh
 real(kind=RKIND),intent(in):: rjday
 real(kind=RKIND),intent(out),dimension(:,:,:):: h2opl

 integer:: n1,n2,nc,l,j,j1,j2
 real(kind=RKIND):: rjday_work,tem,tx1,tx2,t1,t2,cycle_days,tol
 logical:: cyclic

 call mpas_log_write('--- enter subroutine update_h2ophys_forcing:')

 cycle_days = h2ophys%time(h2ophys%ntime+1)-h2ophys%time(1)
 rjday_work = rjday
 if(rjday_work < h2ophys%time(1)) rjday_work = rjday_work + cycle_days
 call find_photochem_time_index(h2ophys%ntime,h2ophys%time,rjday_work,n1,n2,cyclic)

!
!--- Temporary one-time diagnostic of forcing-time interpolation
!
!$OMP CRITICAL(photochem_time_log)
if(log_h2o_time_once) then

   if(cyclic) then
      t1 = h2ophys%time(h2ophys%ntime)
      t2 = h2ophys%time(h2ophys%ntime+1)
   else
      t1 = h2ophys%time(n1)
      t2 = h2ophys%time(n2)
   endif

   tx1 = (t2-rjday_work)/(t2-t1)
   tx2 = 1._RKIND-tx1

   call mpas_log_write('--- H2O forcing-time interpolation check:')
   call mpas_log_write('    rjday = $r, rjday_work = $r', &
                       realArgs=(/rjday,rjday_work/))
   call mpas_log_write('    data indices: n1 = $i, n2 = $i', &
                       intArgs=(/n1,n2/))
   call mpas_log_write('    time coordinates: t1 = $r, t2 = $r', &
                       realArgs=(/t1,t2/))
   call mpas_log_write('    weights: tx1 = $r, tx2 = $r', &
                       realArgs=(/tx1,tx2/))

   if(cyclic) then
      call mpas_log_write('    cyclic interval = TRUE')
   else
      call mpas_log_write('    cyclic interval = FALSE')
   endif

   log_h2o_time_once = .false.
endif
!$OMP END CRITICAL(photochem_time_log)

 if(.not.cyclic) then
    call h2ophys%update(idx1,idx2,idxh,rjday_work,n1,n2,h2opl)
    return
 endif

 t1 = h2ophys%time(h2ophys%ntime)
 t2 = h2ophys%time(h2ophys%ntime+1)
 tol = 100._RKIND*epsilon(1._RKIND)*max(1._RKIND,abs(t2))
 if(rjday_work < t1-tol .or. rjday_work > t2+tol) &
    call physics_error_fatal('update_h2ophys_forcing: model time falls outside the cyclic forcing interval.')
 tx1 = (t2-rjday_work)/(t2-t1)
 tx2 = 1._RKIND-tx1

 do nc=1,h2ophys%ncf
 do l=1,h2ophys%nlev
 do j=1,size(idx1)
    j1 = idx1(j)
    j2 = idx2(j)
    tem = 1._RKIND-idxh(j)
    h2opl(j,l,nc) = tx1*(tem*h2ophys%data(j1,l,nc,h2ophys%ntime) + idxh(j)*h2ophys%data(j2,l,nc,h2ophys%ntime)) &
                   +tx2*(tem*h2ophys%data(j1,l,nc,1)            + idxh(j)*h2ophys%data(j2,l,nc,1))
 enddo
 enddo
 enddo

 call mpas_log_write('--- end subroutine update_h2ophys_forcing:')

 end subroutine update_h2ophys_forcing

!=================================================================================================================
 subroutine find_photochem_time_index(ntime,time,rjday,n1,n2,cyclic)
!=================================================================================================================

!Brackets rjday between two entries of the forcing time axis. This is a copy of the routine of the same
!name contained in GFS_phys_time_vary.fv3.F90, repeated here for two reasons: it is a contained procedure
!there and so cannot be imported, and its dummy argument for the time axis is declared with the default
!real kind, which only matches the CCPP kind_phys because the UFS compiles ccpp-physics with -r8. MPAS
!does not, so the kind is written out explicitly below.

!input arguments:
 integer,intent(in):: ntime
 real(kind=RKIND),intent(in),dimension(ntime+1):: time
 real(kind=RKIND),intent(in):: rjday

!output arguments:
 integer,intent(out):: n1,n2
 logical,intent(out),optional:: cyclic

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
 if(present(cyclic)) cyclic = n2 > ntime
 if(n2 > ntime) n2 = 1

 end subroutine find_photochem_time_index

!=================================================================================================================
 end module mpas_atmphys_photochem_ccpp
!=================================================================================================================
