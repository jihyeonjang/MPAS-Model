!=================================================================================================================
 module mpas_chemistry_driver
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_subpool
 use mpas_timer,only: mpas_timer_start,mpas_timer_stop

 use mpas_chemistry_interface
 use o3_gfs


 implicit none
 private
 public:: chemistry_driver


!MPAS top gocart chemistry driver.
!Laura D. Fowler (laura@ucar.edu) / (2022-02-08).
!
!subroutine chemistry_driver is the top chemistry driver that calls the chemistry drivers for the treatment of
!gocart aerosols,gocart seasalts,and gocart dusts.

!add-ons and modifications to sourcecode:
!----------------------------------------
! * in subroutine chemistry_driver, hit return if config_o3prognostic is set to false.
!   Laura D. Fowler (laura@ucar.edu) / 2021-04-08.


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


 character(len=StrKIND):: errmsg

 logical,pointer:: o3prognostic
 logical:: ldiag3d,qdiag3d

 integer:: i,k,n
 integer:: its,ite,im,levs,kts,kte,no3,no3levs
 integer:: errflg
 integer:: time_lev

 real(kind=RKIND):: dt,grav
 real(kind=RKIND),dimension(:,:),pointer:: qoz,temp

!-----------------------------------------------------------------------------------------------------------------

!--- hit return if config_o3prognostic is set to false:
 call mpas_pool_get_config(domain%configs,'config_o3prognostic',o3prognostic)
 if(.not. o3prognostic) return


 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_driver:')

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


    !--- defines dimensions used to run CHEM2D-OPP:
    call mpas_chems%chemistry_dims(mesh)


    !--- allocates local arrays needed to run CHEM2D-OPP:
    call mpas_chems%chemistry_allocate()


    !--- initializes local arrays needed to run CHEM2D-OPP:
    call mpas_chems%from_MPAS(block%configs,mesh,chem_input,diag,state,time_lev)


    !--- call to OPP:
    errflg = 0
    errmsg = ' '

    ldiag3d = .true.
    qdiag3d = .true.

    call mpas_log_write('mpas_chems%im      = $i',intArgs=(/mpas_chems%im/)     )
    call mpas_log_write('mpas_chems%levs    = $i',intArgs=(/mpas_chems%levs/)   )
    call mpas_log_write('mpas_chems%no3     = $i',intArgs=(/mpas_chems%no3/)    )
    call mpas_log_write('mpas_chems%no3levs = $i',intArgs=(/mpas_chems%no3levs/))
    call mpas_log_write('mpas_chems%dt_c    = $r',realArgs=(/mpas_chems%dt_c/)  )
    call mpas_log_write('mpas_chems%grav_c  = $r',realArgs=(/mpas_chems%grav_c/))

    its     = mpas_chems%its
    ite     = mpas_chems%ite
    kts     = mpas_chems%kts
    kte     = mpas_chems%kte
    no3     = mpas_chems%no3
    no3levs = mpas_chems%no3levs

    call o3_gfs_run(mpas_chems%im,        &
                    mpas_chems%levs,      &
                    mpas_chems%no3levs,   &
                    mpas_chems%dt_c,      &
                    mpas_chems%qoz,       &
                    mpas_chems%temp_c,    &
                    mpas_chems%qozp,      &
                    mpas_chems%pres_c,    &
                    mpas_chems%qozcoeffs, &
                    mpas_chems%no3,       &
                    mpas_chems%dpres_c,   &
                    ldiag3d,              &
                    qdiag3d,              &
                    mpas_chems%qozd1,     &
                    mpas_chems%qozd2,     &
                    mpas_chems%qozd3,     &
                    mpas_chems%qozd4,     &
                    mpas_chems%grav_c,    &
                    errmsg,               &
                    errflg                &
                   )


    !--- returns outputs from CHEM2D-OPP:
    call mpas_chems%to_MPAS(mesh,diag_chemistry,tend_chemistry)


    !--- deallocates local arrays used to run CHEM2D-OPP:
    call mpas_chems%chemistry_deallocate()

    block => block % next
 end do
 call mpas_timer_stop('chemistry driver')


 call mpas_log_write('--- end subroutine mpas_chemistry_driver:')
 call mpas_log_write(' ')

 end subroutine chemistry_driver

!=================================================================================================================
 end module mpas_chemistry_driver
!=================================================================================================================

