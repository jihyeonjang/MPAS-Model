!=================================================================================================================
 module mpas_chemistry_driver
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_subpool
 use mpas_timer,only: mpas_timer_start,mpas_timer_stop

 use mpas_chemistry_interface
 use mpas_chemistry_ozphys,only: ozphys

!from the CCPP physics repository (see src/core_atmosphere/Externals.cfg):
 use machine,only: kind_phys


 implicit none
 private
 public:: chemistry_driver


!MPAS top chemistry driver.
!Laura D. Fowler (laura@ucar.edu) / (2022-02-08).
!
!subroutine chemistry_driver is the top chemistry driver. It currently calls the ozone photochemistry
!parameterization of the CCPP physics library (physics/photochem/module_ozphys.F90).

!add-ons and modifications to sourcecode:
!----------------------------------------
! * in subroutine chemistry_driver, hit return if config_o3prognostic is set to false.
!   Laura D. Fowler (laura@ucar.edu) / 2021-04-08.
! * replaced the local copy of the GFS ozone photochemistry with the CCPP ty_ozphys implementation. The
!   run procedure is selected with config_o3prognostic_scheme.


 contains


!=================================================================================================================
 subroutine chemistry_driver(domain)
!=================================================================================================================

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables and pointers:
 type(mpas_pool_type),pointer:: configs
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: chem_input
 type(mpas_pool_type),pointer:: state
 type(mpas_pool_type),pointer:: diag
 type(mpas_pool_type),pointer:: diag_chemistry
 type(mpas_pool_type),pointer:: tend_chemistry
 type(block_type),pointer:: block

 type(chems_vars):: mpas_chems


 character(len=StrKIND),pointer:: o3prognostic_scheme

 logical,pointer:: o3prognostic

 integer:: i,k,n
 integer:: its,ite,im,levs,kts,kte,no3,no3levs
 integer:: time_lev

 real(kind=kind_phys):: con_1ovg

! real(kind=RKIND):: dt,grav
! real(kind=RKIND),dimension(:,:),pointer:: qoz,temp

!-----------------------------------------------------------------------------------------------------------------

!--- hit return if config_o3prognostic is set to false:
 call mpas_pool_get_config(domain%configs,'config_o3prognostic',o3prognostic)

 if(.not. o3prognostic) return

 call mpas_pool_get_config(domain%configs,'config_o3prognostic_scheme',o3prognostic_scheme)

! call mpas_log_write(' ')
! call mpas_log_write('--- enter subroutine mpas_chemistry_driver:')

 time_lev = 1

 call mpas_timer_start('chemistry driver')
 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'          ,mesh          )
    call mpas_pool_get_subpool(block%structs,'chem_input'    ,chem_input    )
    call mpas_pool_get_subpool(block%structs,'state'         ,state         )
    call mpas_pool_get_subpool(block%structs,'diag'          ,diag          )
    call mpas_pool_get_subpool(block%structs,'diag_chemistry',diag_chemistry)
    call mpas_pool_get_subpool(block%structs,'tend_chemistry',tend_chemistry)


    !--- defines dimensions used to run the ozone photochemistry parameterization:
    call mpas_chems%chemistry_dims(mesh)


    !--- allocates local arrays needed to run the ozone photochemistry parameterization:
    call mpas_chems%chemistry_allocate()


    !--- initializes local arrays needed to run the ozone photochemistry parameterization:
    call mpas_chems%from_MPAS(block%configs,mesh,chem_input,diag,state,time_lev)


    !--- call to the CCPP ozone photochemistry parameterization. con_1ovg is the inverse of the
    !    gravitational acceleration, and the four optional arguments are the contributions to the ozone
    !    tendency from production and loss, from the ozone mixing ratio, from temperature, and from the
    !    overhead ozone column:
    call mpas_log_write('mpas_chems%im      = $i',intArgs=(/mpas_chems%im/)     )
    call mpas_log_write('mpas_chems%levs    = $i',intArgs=(/mpas_chems%levs/)   )
    call mpas_log_write('mpas_chems%no3     = $i',intArgs=(/mpas_chems%no3/)    )
    call mpas_log_write('mpas_chems%no3levs = $i',intArgs=(/mpas_chems%no3levs/))
    call mpas_log_write('mpas_chems%dt_c    = $r',realArgs=(/real(mpas_chems%dt_c,RKIND)/)  )
    call mpas_log_write('mpas_chems%grav_c  = $r',realArgs=(/real(mpas_chems%grav_c,RKIND)/))

    its     = mpas_chems%its
    ite     = mpas_chems%ite
    kts     = mpas_chems%kts
    kte     = mpas_chems%kte
    no3     = mpas_chems%no3
    no3levs = mpas_chems%no3levs

    con_1ovg = 1._kind_phys / mpas_chems%grav_c

    select case(trim(o3prognostic_scheme))

       case('oz_phys_2015')
          call ozphys%run_o3prog_2015(con_1ovg,           &
                                      mpas_chems%dt_c,    &
                                      mpas_chems%pres_c,  &
                                      mpas_chems%temp_c,  &
                                      mpas_chems%dpres_c, &
                                      mpas_chems%qozcoeffs, &
                                      mpas_chems%qoz,     &
                                      mpas_chems%qozd1,   &
                                      mpas_chems%qozd2,   &
                                      mpas_chems%qozd3,   &
                                      mpas_chems%qozd4    &
                                     )

       case('oz_phys_2006')
          call ozphys%run_o3prog_2006(con_1ovg,           &
                                      mpas_chems%dt_c,    &
                                      mpas_chems%pres_c,  &
                                      mpas_chems%temp_c,  &
                                      mpas_chems%dpres_c, &
                                      mpas_chems%qozcoeffs, &
                                      mpas_chems%qoz,     &
                                      mpas_chems%qozd1,   &
                                      mpas_chems%qozd2,   &
                                      mpas_chems%qozd3,   &
                                      mpas_chems%qozd4    &
                                     )

       case default
          call mpas_log_write('--- config_o3prognostic_scheme = '//trim(o3prognostic_scheme), &
                              messageType=MPAS_LOG_CRIT)

    end select


    !--- returns outputs from the ozone photochemistry parameterization:
    call mpas_chems%to_MPAS(mesh,diag_chemistry,tend_chemistry)


    !--- deallocates local arrays used to run the ozone photochemistry parameterization:
    call mpas_chems%chemistry_deallocate()

    block => block % next
 end do
 call mpas_timer_stop('chemistry driver')


! call mpas_log_write('--- end subroutine mpas_chemistry_driver:')

 end subroutine chemistry_driver

!=================================================================================================================
 end module mpas_chemistry_driver
!=================================================================================================================

