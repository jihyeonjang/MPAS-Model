!=================================================================================================================
 module mpas_atmphys_driver_photochem
!=================================================================================================================
!
! Driver for the CCPP ozone and stratospheric water-vapor photochemistry schemes.
!
! The CCPP schemes use model temperature, pressure, latitude, time, and tracer
! mixing ratios to calculate photochemical changes in ozone and water vapor.
! MPAS converts these mixing-ratio changes into physics tendencies and passes
! them to the dynamics.
!
! The main MPAS fields used here are:
!   t_p           temperature                                      [K]
!   pres_hyd_p    hydrostatic layer pressure                       [Pa]
!   pres2_hydd_p  dry hydrostatic interface pressure               [Pa]
!   qv_p          water-vapor mixing ratio                         [kg/kg]
!   qo3_p         ozone mixing ratio                               [kg/kg]
!
! pres_hyd_p is used as the pressure coordinate for interpolating the
! photochemical coefficients. The dry pressure thickness from pres2_hydd_p
! is used for the ozone-column calculation.
!
!=================================================================================================================

 use mpas_kind_types
 use mpas_pool_routines
 use mpas_log

 use mpas_atmphys_constants,only: gravity,degrad
 use mpas_atmphys_manager,only: curr_julday
 use mpas_atmphys_utilities,only: physics_error_fatal
 use mpas_atmphys_vars

 use mpas_atmphys_photochem_ccpp

 implicit none
 private

 public:: init_photochem,       &
          driver_photochem

 contains

!=================================================================================================================
 subroutine init_photochem(configs,mesh,diag_physics)
!=================================================================================================================

!Reads the forcing files and computes, once, the latitude interpolation indices and weights of every cell.
!Called from physics_init.

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: o3_scheme,h2o_scheme
 character(len=StrKIND),pointer:: o3file,h2ofile
 integer,pointer:: nCellsSolve
 real(kind=RKIND),dimension(:),pointer:: latCell
 integer,dimension(:),pointer:: jindx1_o3,jindx2_o3,jindx1_h2o,jindx2_h2o
 real(kind=RKIND),dimension(:),pointer:: ddy_o3,ddy_h2o

!local variables:
 real(kind=RKIND),dimension(:),allocatable:: lat_d

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_photochem_o3_scheme' ,o3_scheme )
 call mpas_pool_get_config(configs,'config_photochem_h2o_scheme',h2o_scheme)
 call mpas_pool_get_config(configs,'config_o3prdlos_file'       ,o3file    )
 call mpas_pool_get_config(configs,'config_h2oprdlos_file'      ,h2ofile   )

 if(trim(o3_scheme) == 'off' .and. trim(h2o_scheme) == 'off') return

 call mpas_log_write('')
 call mpas_log_write('--- enter subroutine init_photochem:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_array(mesh,'latCell',latCell)

!--- the CCPP setup procedures take the latitude in degrees, as the GFS passes xlat_d:
 allocate(lat_d(1:nCellsSolve))
 lat_d(1:nCellsSolve) = latCell(1:nCellsSolve) / degrad

 if(trim(o3_scheme) /= 'off') then
    call init_ozphys(o3_scheme,o3file)
    call mpas_pool_get_array(diag_physics,'jindx1_o3',jindx1_o3)
    call mpas_pool_get_array(diag_physics,'jindx2_o3',jindx2_o3)
    call mpas_pool_get_array(diag_physics,'ddy_o3'   ,ddy_o3   )
    call ozphys%setup_o3prog(lat_d,jindx1_o3(1:nCellsSolve),jindx2_o3(1:nCellsSolve), &
                             ddy_o3(1:nCellsSolve))
 endif

 if(trim(h2o_scheme) /= 'off') then
    call init_h2ophys(h2o_scheme,h2ofile)
    call mpas_pool_get_array(diag_physics,'jindx1_h2o',jindx1_h2o)
    call mpas_pool_get_array(diag_physics,'jindx2_h2o',jindx2_h2o)
    call mpas_pool_get_array(diag_physics,'ddy_h2o'   ,ddy_h2o   )
    call h2ophys%setup(lat_d,jindx1_h2o(1:nCellsSolve),jindx2_h2o(1:nCellsSolve), &
                       ddy_h2o(1:nCellsSolve))
 endif

 deallocate(lat_d)

 call mpas_log_write('--- end subroutine init_photochem.')

 end subroutine init_photochem

!=================================================================================================================
 subroutine driver_photochem(configs,diag_physics,tend_physics,its,ite,run_o3,run_h2o)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
!--- its and ite shadow the module variables of the same name, so that this driver can be called from
!    inside the threaded loop of physics_driver, one chunk of cells at a time:
 integer,intent(in):: its,ite
 logical,intent(in):: run_o3,run_h2o

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 character(len=StrKIND),pointer:: o3_scheme,h2o_scheme
 logical,pointer:: photochem_diagnostics
 integer,dimension(:),pointer:: jindx1_o3,jindx2_o3,jindx1_h2o,jindx2_h2o
 real(kind=RKIND),dimension(:),pointer:: ddy_o3,ddy_h2o
 real(kind=RKIND),dimension(:,:),pointer:: rqo3phten,rqvphten
 real(kind=RKIND),dimension(:,:),pointer:: do3_dt_prd,do3_dt_ozmx,do3_dt_temp,do3_dt_ohoz
 real(kind=RKIND),dimension(:,:),pointer:: dqv_dt_prd,dqv_dt_qvmx

!local variables:
 logical:: o3_on,h2o_on,diag_o3,diag_h2o
 integer:: i,k,nCols,nLevs
 real(kind=RKIND):: rjday,con_1ovg,rdt

 real(kind=RKIND),dimension(:,:),allocatable:: p_c,t_c,dp_c,oz_c,h2o_c
 real(kind=RKIND),dimension(:,:),allocatable:: d1,d2,d3,d4
 real(kind=RKIND),dimension(:,:,:),allocatable:: ozpl,h2opl

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_photochem_o3_scheme'    ,o3_scheme            )
 call mpas_pool_get_config(configs,'config_photochem_h2o_scheme'   ,h2o_scheme           )
 call mpas_pool_get_config(configs,'config_photochem_diagnostics'  ,photochem_diagnostics)

!--- flags for o3, h2o, and photochem_diagnostics
 o3_on    = run_o3  .and. trim(o3_scheme)  /= 'off'
 h2o_on   = run_h2o .and. trim(h2o_scheme) /= 'off'
 diag_o3  = photochem_diagnostics .and. o3_on
 diag_h2o = photochem_diagnostics .and. h2o_on

 if(.not.o3_on .and. .not.h2o_on) return

 nCols = ite-its+1
 nLevs = kte-kts+1
 con_1ovg = 1._RKIND / gravity
 rdt = 1._RKIND / dt_dyn

!--- Match the UFS host convention: rjday is the Gregorian day of year plus the fractional day,
!    counted from 1.0 at 00Z on 1 January. No leap-year rescaling or leap-day hold is applied.
 rjday = curr_julday + 1._RKIND

!=================================================================================================================
!--- ozone photochemistry
!=================================================================================================================

 if(o3_on) then

    if(.not.ozphys_loaded) &
       call physics_error_fatal('subroutine driver_photochem: the ozone forcing data were not loaded.')

    call mpas_pool_get_array(diag_physics,'jindx1_o3',jindx1_o3)
    call mpas_pool_get_array(diag_physics,'jindx2_o3',jindx2_o3)
    call mpas_pool_get_array(diag_physics,'ddy_o3'   ,ddy_o3   )
    call mpas_pool_get_array(tend_physics,'rqo3phten',rqo3phten)

    if(diag_o3) then
       call mpas_pool_get_array(diag_physics,'do3_dt_prd' ,do3_dt_prd )
       call mpas_pool_get_array(diag_physics,'do3_dt_ozmx',do3_dt_ozmx)
       call mpas_pool_get_array(diag_physics,'do3_dt_temp',do3_dt_temp)
       call mpas_pool_get_array(diag_physics,'do3_dt_ohoz',do3_dt_ohoz)
    endif

    allocate(p_c(nCols,nLevs),t_c(nCols,nLevs),dp_c(nCols,nLevs),oz_c(nCols,nLevs))
    call photochem_o3_from_MPAS(its,ite,nCols,nLevs,p_c,t_c,dp_c,oz_c)
    allocate(ozpl(nCols,ozphys%nlev,ozphys%ncf))

!--- Interpolate the forcing coefficients in latitude and time. The latitude indices and weights were
!    computed once, at initialization; the time interpolation is redone every time step.
    call update_ozphys_forcing(jindx1_o3(its:ite),jindx2_o3(its:ite),ddy_o3(its:ite),rjday,ozpl)

    if(diag_o3) then
       allocate(d1(nCols,nLevs),d2(nCols,nLevs),d3(nCols,nLevs),d4(nCols,nLevs))
       d1 = 0._RKIND
       d2 = 0._RKIND
       d3 = 0._RKIND
       d4 = 0._RKIND
    endif

    select case(trim(o3_scheme))

       case('oz_phys_2015')
!--- the 2015 scheme needs all six coefficients and returns all four diagnostic terms.
          if(diag_o3) then
             call ozphys%run_o3prog_2015(con_1ovg,dt_dyn,p_c,t_c,dp_c,ozpl,oz_c,d1,d2,d3,d4)
          else
             call ozphys%run_o3prog_2015(con_1ovg,dt_dyn,p_c,t_c,dp_c,ozpl,oz_c)
          endif

       case('oz_phys_2006')
!--- the 2006 scheme runs with either two or four coefficients, and init_ozphys rejects anything else.
!    With two coefficients the temperature and overhead ozone terms do not exist: the CCPP code would
!    index prod(:,3) and read an uninitialized colo3, so those two optional diagnostics must only be
!    passed when there are four.
          if(diag_o3) then
             if(ozphys%ncf == 4) then
                call ozphys%run_o3prog_2006(con_1ovg,dt_dyn,p_c,t_c,dp_c,ozpl,oz_c,d1,d2,d3,d4)
             else
                call ozphys%run_o3prog_2006(con_1ovg,dt_dyn,p_c,t_c,dp_c,ozpl,oz_c,d1,d2)
             endif
          else
             call ozphys%run_o3prog_2006(con_1ovg,dt_dyn,p_c,t_c,dp_c,ozpl,oz_c)
          endif

       case default
          call physics_error_fatal('subroutine driver_photochem: unknown config_photochem_o3_scheme = ' &
                                   //trim(o3_scheme))

    end select

!--- The CCPP scheme returns an updated mixing ratio. Write the corresponding tendency directly into the
!    Registry-owned MPAS array instead of staging it in a second full-size module array and copying it later.
    do i = its,ite
    do k = kts,kte
       rqo3phten(k,i) = (oz_c(i-its+1,k-kts+1) - qo3_p(i,k,jts)) * rdt
       if(diag_o3) then
          do3_dt_prd(k,i)  = d1(i-its+1,k-kts+1) * rdt
          do3_dt_ozmx(k,i) = d2(i-its+1,k-kts+1) * rdt
          do3_dt_temp(k,i) = d3(i-its+1,k-kts+1) * rdt
          do3_dt_ohoz(k,i) = d4(i-its+1,k-kts+1) * rdt
       endif
    enddo
    enddo

    if(diag_o3) deallocate(d1,d2,d3,d4)
    deallocate(ozpl,oz_c,p_c,t_c,dp_c)

 endif

!=================================================================================================================
!--- stratospheric water vapor photochemistry
!=================================================================================================================

 if(h2o_on) then

    if(.not.h2ophys_loaded) &
       call physics_error_fatal('subroutine driver_photochem: the water vapor forcing data were not loaded.')

    call mpas_pool_get_array(diag_physics,'jindx1_h2o',jindx1_h2o)
    call mpas_pool_get_array(diag_physics,'jindx2_h2o',jindx2_h2o)
    call mpas_pool_get_array(diag_physics,'ddy_h2o'   ,ddy_h2o   )
    call mpas_pool_get_array(tend_physics,'rqvphten',rqvphten)

    if(diag_h2o) then
       call mpas_pool_get_array(diag_physics,'dqv_dt_prd' ,dqv_dt_prd )
       call mpas_pool_get_array(diag_physics,'dqv_dt_qvmx',dqv_dt_qvmx)
    endif

    allocate(p_c(nCols,nLevs),h2o_c(nCols,nLevs))
    call photochem_h2o_from_MPAS(its,ite,nCols,nLevs,p_c,h2o_c)
    allocate(h2opl(nCols,h2ophys%nlev,h2ophys%ncf))

    call update_h2ophys_forcing(jindx1_h2o(its:ite),jindx2_h2o(its:ite),ddy_h2o(its:ite),rjday,h2opl)

    if(diag_h2o) then
       allocate(d1(nCols,nLevs),d2(nCols,nLevs))
       d1 = 0._RKIND
       d2 = 0._RKIND
       call h2ophys%run(dt_dyn,p_c,h2opl,h2o_c,d1,d2)
    else
       call h2ophys%run(dt_dyn,p_c,h2opl,h2o_c)
    endif

    do i = its,ite
    do k = kts,kte
       rqvphten(k,i) = (h2o_c(i-its+1,k-kts+1) - qv_p(i,k,jts)) * rdt
       if(diag_h2o) then
          dqv_dt_prd(k,i)  = d1(i-its+1,k-kts+1) * rdt
          dqv_dt_qvmx(k,i) = d2(i-its+1,k-kts+1) * rdt
       endif
    enddo
    enddo

    if(diag_h2o) deallocate(d1,d2)
    deallocate(h2opl,h2o_c,p_c)

 endif

 end subroutine driver_photochem

!=================================================================================================================
 subroutine photochem_h2o_from_MPAS(its,ite,nCols,nLevs,p_c,h2o_c)
!=================================================================================================================

!Packs the MPAS fields required by the stratospheric water-vapor photochemistry scheme: pressure and
!water-vapor mixing ratio. The CCPP arrays have lower bounds of 1, so the MPAS indices are re-based here.

!input arguments:
 integer,intent(in):: its,ite,nCols,nLevs

!output arguments:
 real(kind=RKIND),intent(out),dimension(nCols,nLevs):: p_c,h2o_c

!local variables:
 integer:: i,k

!-----------------------------------------------------------------------------------------------------------------

 do k = kts,kte
 do i = its,ite
    p_c  (i-its+1,k-kts+1) = pres_hyd_p(i,k,jts)
    h2o_c(i-its+1,k-kts+1) = qv_p(i,k,jts)
 enddo
 enddo

 end subroutine photochem_h2o_from_MPAS

!=================================================================================================================
 subroutine photochem_o3_from_MPAS(its,ite,nCols,nLevs,p_c,t_c,dp_c,oz_c)
!=================================================================================================================

!Packs the MPAS fields required by the ozone photochemistry scheme: pressure, temperature, dry layer
!pressure thickness, and transported ozone mixing ratio. The CCPP arrays have lower bounds of 1, so the
!MPAS indices are re-based here.

!input arguments:
 integer,intent(in):: its,ite,nCols,nLevs

!output arguments:
 real(kind=RKIND),intent(out),dimension(nCols,nLevs):: p_c,t_c,dp_c,oz_c

!local variables:
 integer:: i,k

!-----------------------------------------------------------------------------------------------------------------

 do k = kts,kte
 do i = its,ite
    p_c (i-its+1,k-kts+1) = pres_hyd_p(i,k,jts)
    t_c (i-its+1,k-kts+1) = t_p(i,k,jts)
!--- dry layer thickness, positive: the "dry" hydrostatic pressure decreases upward.
    dp_c(i-its+1,k-kts+1) = pres2_hydd_p(i,k,jts) - pres2_hydd_p(i,k+1,jts)
    oz_c(i-its+1,k-kts+1) = qo3_p(i,k,jts)
 enddo
 enddo

 end subroutine photochem_o3_from_MPAS

!=================================================================================================================
 end module mpas_atmphys_driver_photochem
!=================================================================================================================
