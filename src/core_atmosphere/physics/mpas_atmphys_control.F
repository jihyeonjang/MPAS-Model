! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_control
 use mpas_dmpar
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_utilities
 use mpas_atmphys_vars, only: l_mp_tables

 implicit none
 private
 public:: physics_namelist_check, &
          physics_registry_init,  &
          physics_tables_init,    &
          physics_compatibility_check

 logical,public:: moist_physics


!MPAS control and initialization routines.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines called in mpas_atmphys_control:
! -------------------------------------------
! * physics_namelist_check: checks that physics namelist parameters are defined correctly.
! * physics_registry_init : initializes thicknesses of soil layers for NOAH land-surface scheme.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the namelist option config_eddy_scheme and associated sourcecode.
! * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
! * removed controls to the updated Kain-Fritsch convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * removed subroutine physics_idealized_init, also available in mpas_init_atm_cases.F in core_init_atmosphere.
!   Laura D. Fowler (laura@ucar.edu) / 2014-08-11.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * renamed "wsm6" to "mp_wsm6" and "kessler" to "mp_kessler".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-09.
! * renamed "kain_fritsch" to "cu_kain_fritsch" and "tiedtke" to "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * renamed "ysu" to "bl_ysu", "ysu_gwdo" to "bl_gwdo_ysu", and "monin_obukhov" to "sf_monin_obukhov".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25. 
! * added the option mp_thompson.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25.
! * added the option cu_grell_freitas.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-31.
! * added the options sf_mynn and bl_mynn and for the MYNN parameterization from WRF version 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.
! * added the option cu_ntiedtke for the "new" Tiedtke parameterization of convection from WRF version 3.8.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-19.
! * added the physics suite "convection_scale_aware" (see below for the physics options used in the suite).
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-28.
! * added the subroutine physics_tables_init which checks if the files containing the lokk-up tables for the
!   Thompson cloud microphysics are available or not.
!   Laura D. Fowler (laura@ucar.edu) / 2016-11-01.
! * modified checking the config_gwdo_scheme option to allow bl_ysu_gwdo to be run when the MYNN pbl and surface
!   layer scheme options are chosen.
!   Laura D. Fowler (laura@ucar.edu) / 2016-12-22.
! * modified logic in subroutine physics_tables_init so that the Thompson microphysics tables are read in each
!   MPI task.
!   Laura D. Fowler (laura@ucar.edu) / 2016-12-30.
! * added the option mp_thompson_aerosols.
!   Laura D. Fowler (laura@ucar.edu) / 2018-01-31.
! * added the option sf_monin_obukhov_rev to run the revised surface layer scheme with the YSU PBL scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2023-05-15.
! * replaced the option "noah" with "sf_noah" to run the NOAH land surface scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2022-02-18.
! * added the option "sf_noahmp" to run the NOAH-MP land surface scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2022-07-15.
! * in the mesoscale_reference suite, replaced the MM5 surface layer scheme with the MM5 revised surface layer
!   scheme as the default option for config_sfclayer_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2024-06-18.


 contains


!=================================================================================================================
 subroutine physics_namelist_check(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: config_physics_suite,     &
                                  config_microp_scheme,     &
                                  config_convection_scheme, &
                                  config_lsm_scheme,        &
                                  config_pbl_scheme,        &
                                  config_gwdo_scheme,       &
                                  config_radt_cld_scheme,   &
                                  config_radt_lw_scheme,    &
                                  config_radt_sw_scheme,    &
                                  config_sfclayer_scheme

 character(len=StrKIND),pointer:: config_o3_source,            &
                                  config_photochem_o3_scheme,  &
                                  config_photochem_h2o_scheme

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine physics_namelist_check:')

 call mpas_pool_get_config(configs,'config_physics_suite'    ,config_physics_suite    )
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )
 call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call mpas_pool_get_config(configs,'config_pbl_scheme'       ,config_pbl_scheme       )
 call mpas_pool_get_config(configs,'config_gwdo_scheme'      ,config_gwdo_scheme      )
 call mpas_pool_get_config(configs,'config_radt_cld_scheme'  ,config_radt_cld_scheme  )
 call mpas_pool_get_config(configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )
 call mpas_pool_get_config(configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )

 call mpas_pool_get_config(configs,'config_o3_source'            ,config_o3_source            )
 call mpas_pool_get_config(configs,'config_photochem_o3_scheme'  ,config_photochem_o3_scheme  )
 call mpas_pool_get_config(configs,'config_photochem_h2o_scheme' ,config_photochem_h2o_scheme )

 call mpas_log_write('')
 call mpas_log_write('----- Setting up physics suite '''//trim(config_physics_suite)//''' -----')

 !
 !setup schemes according to the selected physics suite:
 !
 if (trim(config_physics_suite) == 'mesoscale_reference') then

    if (trim(config_microp_scheme)     == 'suite') config_microp_scheme     = 'mp_wsm6'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'cu_ntiedtke'
    if (trim(config_pbl_scheme)        == 'suite') config_pbl_scheme        = 'bl_ysu'
    if (trim(config_gwdo_scheme)       == 'suite') config_gwdo_scheme       = 'bl_ysu_gwdo'
    if (trim(config_radt_lw_scheme)    == 'suite') config_radt_lw_scheme    = 'rrtmg_lw'
    if (trim(config_radt_sw_scheme)    == 'suite') config_radt_sw_scheme    = 'rrtmg_sw'
    if (trim(config_radt_cld_scheme)   == 'suite') config_radt_cld_scheme   = 'cld_fraction'
    if (trim(config_sfclayer_scheme)   == 'suite') config_sfclayer_scheme   = 'sf_monin_obukhov_rev'
    if (trim(config_lsm_scheme)        == 'suite') config_lsm_scheme        = 'sf_noah'

 else if (trim(config_physics_suite) == 'convection_permitting') then

    if (trim(config_microp_scheme)     == 'suite') config_microp_scheme     = 'mp_thompson'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'cu_grell_freitas'
    if (trim(config_pbl_scheme)        == 'suite') config_pbl_scheme        = 'bl_mynn'
    if (trim(config_gwdo_scheme)       == 'suite') config_gwdo_scheme       = 'bl_ysu_gwdo'
    if (trim(config_radt_lw_scheme)    == 'suite') config_radt_lw_scheme    = 'rrtmg_lw'
    if (trim(config_radt_sw_scheme)    == 'suite') config_radt_sw_scheme    = 'rrtmg_sw'
    if (trim(config_radt_cld_scheme)   == 'suite') config_radt_cld_scheme   = 'cld_fraction'
    if (trim(config_sfclayer_scheme)   == 'suite') config_sfclayer_scheme   = 'sf_mynn'
    if (trim(config_lsm_scheme)        == 'suite') config_lsm_scheme        = 'sf_noah'

 else if (trim(config_physics_suite) == 'none') then

    if (trim(config_microp_scheme)     == 'suite') config_microp_scheme     = 'off'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'off'
    if (trim(config_pbl_scheme)        == 'suite') config_pbl_scheme        = 'off'
    if (trim(config_gwdo_scheme)       == 'suite') config_gwdo_scheme       = 'off'
    if (trim(config_radt_lw_scheme)    == 'suite') config_radt_lw_scheme    = 'off'
    if (trim(config_radt_sw_scheme)    == 'suite') config_radt_sw_scheme    = 'off'
    if (trim(config_radt_cld_scheme)   == 'suite') config_radt_cld_scheme   = 'off'
    if (trim(config_sfclayer_scheme)   == 'suite') config_sfclayer_scheme   = 'off'
    if (trim(config_lsm_scheme)        == 'suite') config_lsm_scheme        = 'off'

 else

    write(mpas_err_message,'(A)') 'Unrecognized choice of physics suite: config_physics_suite = '''// &
                                  trim(config_physics_suite)//''''
    call physics_error_fatal(mpas_err_message)

 end if

!cloud microphysics scheme:
 if(.not. (config_microp_scheme .eq. 'off'                  .or. &
           config_microp_scheme .eq. 'mp_kessler'           .or. &
           config_microp_scheme .eq. 'mp_thompson'          .or. &
           config_microp_scheme .eq. 'mp_thompson_aerosols' .or. &
           config_microp_scheme .eq. 'mp_wsm6')) then

    write(mpas_err_message,'(A,A20)') 'illegal value for config_microp_scheme:', &
          trim(config_microp_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!convection scheme:
 if(.not. (config_convection_scheme .eq. 'off'              .or. &
           config_convection_scheme .eq. 'cu_grell_freitas' .or. &
           config_convection_scheme .eq. 'cu_kain_fritsch'  .or. &
           config_convection_scheme .eq. 'cu_tiedtke'       .or. &
           config_convection_scheme .eq. 'cu_ntiedtke')) then

    write(mpas_err_message,'(A,A20)') 'illegal value for config_convection_scheme: ', &
          trim(config_convection_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!pbl scheme:
 if(.not. (config_pbl_scheme .eq. 'off'     .or. &
           config_pbl_scheme .eq. 'bl_mynn' .or. &
           config_pbl_scheme .eq. 'bl_ysu')) then

    write(mpas_err_message,'(A,A20)') 'illegal value for pbl_scheme: ', &
          trim(config_pbl_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!gravity wave drag over orography scheme:
 if(.not. (config_gwdo_scheme .eq. 'off' .or. &
           config_gwdo_scheme .eq. 'bl_ysu_gwdo' .or. &
           config_gwdo_scheme .eq. 'bl_ugwp_gwdo')) then

    write(mpas_err_message,'(A,A20)') 'illegal value for gwdo_scheme: ', &
          trim(config_gwdo_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!lw radiation scheme:
 if(.not. (config_radt_lw_scheme .eq. 'off'    .or. &
           config_radt_lw_scheme .eq. 'cam_lw' .or. &
           config_radt_lw_scheme .eq. 'rrtmg_lw')) then
 
    write(mpas_err_message,'(A,A20)') 'illegal value for longwave radiation scheme: ', &
          trim(config_radt_lw_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!sw radiation scheme:
 if(.not. (config_radt_sw_scheme .eq. 'off'    .or. &
           config_radt_sw_scheme .eq. 'cam_sw' .or. &
           config_radt_sw_scheme .eq. 'rrtmg_sw')) then
 
    write(mpas_err_message,'(A,A20)') 'illegal value for shortwave radiation _scheme: ', &
          trim(config_radt_sw_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!cloud fraction for radiation schemes:
 if(.not. (config_radt_cld_scheme .eq. 'off'           .or. &
           config_radt_cld_scheme .eq. 'cld_incidence' .or. &
           config_radt_cld_scheme .eq. 'cld_fraction'  .or. &
           config_radt_cld_scheme .eq. 'cld_fraction_thompson')) then

    write(mpas_err_message,'(A,A20)') 'illegal value for calculation of cloud fraction: ', &
          trim(config_radt_cld_scheme)
    call physics_error_fatal(mpas_err_message)

 endif
 if((config_radt_lw_scheme.ne.'off' .and. config_radt_cld_scheme.eq.'off') .or. &
    (config_radt_sw_scheme.ne.'off' .and. config_radt_cld_scheme.eq.'off')) then

    call mpas_log_write('')
    write(mpas_err_message,'(A,A20)') &
       '    config_radt_cld_scheme is not set for radiation calculation'
    call physics_message(mpas_err_message)
    write(mpas_err_message,'(A,A20)') &
       '    switch calculation of cloud fraction to config_radt_cld_scheme = cld_incidence'
    call physics_message(mpas_err_message)
    config_radt_cld_scheme = "cld_incidence"

 endif

!surface-layer scheme:
 if(.not. (config_sfclayer_scheme .eq. 'off'              .or. &
           config_sfclayer_scheme .eq. 'sf_mynn'          .or. &
           config_sfclayer_scheme .eq. 'sf_monin_obukhov' .or. &
           config_sfclayer_scheme .eq. 'sf_monin_obukhov_rev')) then
 
    write(mpas_err_message,'(A,A20)') 'illegal value for surface layer scheme: ', &
          trim(config_sfclayer_scheme)
    call physics_error_fatal(mpas_err_message)
 else
    if(config_pbl_scheme == 'bl_mynn') then
       config_sfclayer_scheme = 'sf_mynn'
    elseif(config_pbl_scheme == 'bl_ysu') then
       if(config_sfclayer_scheme /= 'sf_monin_obukhov' .and. &
          config_sfclayer_scheme /= 'sf_monin_obukhov_rev') then
          write(mpas_err_message,'(A,A20)') 'wrong choice for surface layer scheme with YSU PBL: ', &
                trim(config_sfclayer_scheme)
          call physics_error_fatal(mpas_err_message)
       endif
    endif
 endif

!land-surface scheme: note that config_sfclayer_scheme must be defined for the land-surface
!scheme to be called:
 if(config_lsm_scheme .ne. 'off' .and. config_sfclayer_scheme .eq. 'off') then
 
    call physics_error_fatal('land surface scheme: ' // &
                             'set config_sfclayer_scheme different than off')
    
 elseif(.not. (config_lsm_scheme .eq. 'off '    .or. &
               config_lsm_scheme .eq. 'sf_noah' .or. &
               config_lsm_scheme .eq. 'sf_noahmp')) then
 
    write(mpas_err_message,'(A,A20)') 'illegal value for land surface scheme: ', &
          trim(config_lsm_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!checks if any physics process is called. if not, return:
 moist_physics = .true.
 
 if(config_microp_scheme     .eq. 'off' .and. &
    config_convection_scheme .eq. 'off' .and. &
    config_lsm_scheme        .eq. 'off' .and. &
    config_pbl_scheme        .eq. 'off' .and. &
    config_radt_lw_scheme    .eq. 'off' .and. &
    config_radt_sw_scheme    .eq. 'off' .and. &
    config_sfclayer_scheme   .eq. 'off' .and. &
    config_photochem_o3_scheme  .eq. 'off' .and. &
    config_photochem_h2o_scheme .eq. 'off') moist_physics = .false.

 call mpas_log_write('')
 call mpas_log_write('    config_microp_scheme       = '//trim(config_microp_scheme))
 call mpas_log_write('    config_convection_scheme   = '//trim(config_convection_scheme))
 call mpas_log_write('    config_pbl_scheme          = '//trim(config_pbl_scheme))
 call mpas_log_write('    config_gwdo_scheme         = '//trim(config_gwdo_scheme))
 call mpas_log_write('    config_radt_cld_scheme     = '//trim(config_radt_cld_scheme))
 call mpas_log_write('    config_radt_lw_scheme      = '//trim(config_radt_lw_scheme))
 call mpas_log_write('    config_radt_sw_scheme      = '//trim(config_radt_sw_scheme))
 call mpas_log_write('    config_sfclayer_scheme     = '//trim(config_sfclayer_scheme))
 call mpas_log_write('    config_lsm_scheme          = '//trim(config_lsm_scheme))
 call mpas_log_write('')

!
!the ozone and stratospheric water vapor photochemistry:
!
!config_o3_source selects which ozone the RRTMG long- and short-wave codes see (o3input = 2 or 3 in
!mpas_atmphys_driver_radiation_lw/sw.F). An unrecognized value would silently leave o3input = 0, that is
!radiation with no ozone at all, so it is rejected here.
 if(trim(config_o3_source) /= 'climatology' .and. trim(config_o3_source) /= 'prognostic') then
    call mpas_log_write('    config_o3_source           = '//trim(config_o3_source),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('config_o3_source must be either ''climatology'' or ''prognostic''.')
 endif

!the CCPP also provides oz_phys_2006, and it is accepted here, but it is deliberately left out of the
!possible_values of config_photochem_o3_scheme in the Registry: run_o3prog_2006 only updates ozone for a
!forcing file with 2 or 4 coefficients, and the file distributed with the GFS has 6. init_ozphys checks
!that, so an unusable combination stops the model rather than running silently.
 select case(trim(config_photochem_o3_scheme))
    case('off','oz_phys_2015')
    case('oz_phys_2006')
       call mpas_log_write('')
       call mpas_log_write('*** config_photochem_o3_scheme = ''oz_phys_2006'' is not among the ' // &
                           'supported values of that option.')
       call mpas_log_write('*** It needs an ozone forcing file with 2 or 4 coefficients, which is not ' // &
                           'the one distributed with the GFS.')
       call mpas_log_write('')
    case default
       call mpas_log_write('    config_photochem_o3_scheme = '//trim(config_photochem_o3_scheme), &
                           messageType=MPAS_LOG_ERR)
       call physics_error_fatal('config_photochem_o3_scheme must be ''off'' or ''oz_phys_2015''.')
 end select

 select case(trim(config_photochem_h2o_scheme))
    case('off','h2o_phys')
    case default
       call mpas_log_write('    config_photochem_h2o_scheme = '//trim(config_photochem_h2o_scheme), &
                           messageType=MPAS_LOG_ERR)
       call physics_error_fatal('config_photochem_h2o_scheme must be ''off'' or ''h2o_phys''.')
 end select

 if(trim(config_o3_source) == 'climatology' .and. trim(config_photochem_o3_scheme) /= 'off') then
    call mpas_log_write('    config_o3_source           = '//trim(config_o3_source),messageType=MPAS_LOG_ERR)
    call mpas_log_write('    config_photochem_o3_scheme = '//trim(config_photochem_o3_scheme), &
                        messageType=MPAS_LOG_ERR)
    call physics_error_fatal('(1) Climatological ozone requires config_o3_source = ''climatology'' and config_photochem_o3_scheme = ''off'' ; ' // &
                             '(2) Ozone photochemistry requires config_o3_source = ''prognostic''.')
 endif

!--- Ozone photochemistry is intended to advance the prognostic qo3 field that is also seen by radiation.
!    CAM radiation does not consume qo3; allowing CAM in either waveband would therefore create a partially
!    or completely uncoupled ozone-radiation system. Fail fast instead of silently running that scientifically
!    inconsistent configuration. This restriction applies to O3 photochemistry only; H2O photochemistry may
!    still be used with CAM because CAM consumes the model water-vapor field.
 if(trim(config_photochem_o3_scheme) /= 'off' .and. &
    (trim(config_radt_lw_scheme) == 'cam_lw' .or. trim(config_radt_sw_scheme) == 'cam_sw')) then
    call mpas_log_write('    config_radt_lw_scheme       = '//trim(config_radt_lw_scheme),messageType=MPAS_LOG_ERR)
    call mpas_log_write('    config_radt_sw_scheme       = '//trim(config_radt_sw_scheme),messageType=MPAS_LOG_ERR)
    call mpas_log_write('    config_photochem_o3_scheme  = '//trim(config_photochem_o3_scheme),messageType=MPAS_LOG_ERR)
    call physics_error_fatal('ozone photochemistry cannot be combined with CAM radiation because CAM does not ' // &
                             'use the prognostic qo3 field. Use RRTMG for both wavebands, or turn ozone ' // &
                             'photochemistry off.')
 endif

 call mpas_log_write('    config_o3_source            = '//trim(config_o3_source))
 call mpas_log_write('    config_photochem_o3_scheme  = '//trim(config_photochem_o3_scheme))
 call mpas_log_write('    config_photochem_h2o_scheme = '//trim(config_photochem_h2o_scheme))
 call mpas_log_write('')

 if(trim(config_o3_source) == 'prognostic' .and. &
    trim(config_photochem_o3_scheme) == 'off') then
    call mpas_log_write('*** config_o3_source = ''prognostic'' with config_photochem_o3_scheme = ''off'':')
    call mpas_log_write('*** qo3 is read from the initial conditions and transported by the dynamics, but')
    call mpas_log_write('*** nothing changes it by the photochem scheme. The radiation uses that transported field.')
    call mpas_log_write('')
 endif

 end subroutine physics_namelist_check

!=================================================================================================================
 subroutine physics_registry_init(mesh,configs,sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: config_do_restart
 character(len=StrKIND),pointer:: config_lsm_scheme
 integer,pointer:: nCells
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:,:),pointer:: dzs

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)
 call mpas_pool_get_config(configs,'config_lsm_scheme',config_lsm_scheme)

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(sfc_input,'landmask',landmask)
 call mpas_pool_get_array(sfc_input,'dzs'     , dzs    )

!initialization of input variables, if needed:

 if(.not. config_do_restart) then

    lsm_select: select case(trim(config_lsm_scheme))

       case("sf_noah","sf_noahmp")
       !initialize the thickness of the soil layers for the Noah scheme:
          do iCell = 1, nCells
             dzs(1,iCell) = 0.10_RKIND
             dzs(2,iCell) = 0.30_RKIND
             dzs(3,iCell) = 0.60_RKIND
             dzs(4,iCell) = 1.00_RKIND
          enddo

       case default
    
    end select lsm_select
    
 endif

!call mpas_log_write('--- enter subroutine physics_namelist_check.')
!call mpas_log_write('')

 end subroutine physics_registry_init

!=================================================================================================================
 subroutine physics_tables_init(dminfo,configs)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: configs

!local variables:
 character(len=StrKIND),pointer:: config_microp_scheme
 logical:: l_qr_acr_qg,l_qr_acr_qs,l_qi_aut_qs,l_freezeH2O

!-----------------------------------------------------------------------------------------------------------------

 l_mp_tables = .true.

 if(dminfo % my_proc_id == IO_NODE) then

    call mpas_pool_get_config(configs,'config_microp_scheme',config_microp_scheme)
    if(config_microp_scheme /= "mp_thompson" .or. &
       config_microp_scheme /= "mp_thompson_aerosols") return

    l_qr_acr_qg = .false.
    l_qr_acr_qs = .false.
    l_qi_aut_qs = .false.
    l_freezeH2O = .false.

    inquire(file='MP_THOMPSON_QRacrQG_DATA.DBL'  ,exist=l_qr_acr_qg)
    inquire(file='MP_THOMPSON_QRacrQS_DATA.DBL'  ,exist=l_qr_acr_qs)
    inquire(file='MP_THOMPSON_QIautQS_DATA.DBL'  ,exist=l_qi_aut_qs)
    inquire(file='MP_THOMPSON_freezeH2O_DATA.DBL',exist=l_freezeH2O)

!   call mpas_log_write('')
!   call mpas_log_write('--- enter subroutine physics_tables_init:')
!   call mpas_log_write('l_qr_acr_qg = $l',logicArgs=(/l_qr_acr_qg/))
!   call mpas_log_write('l_qr_acr_qs = $l',logicArgs=(/l_qr_acr_qs/))
!   call mpas_log_write('l_qi_aut_qs = $l',logicArgs=(/l_qi_aut_qs/))
!   call mpas_log_write('l_freezeH2O = $l',logicArgs=(/l_freezeH2O/))

    if(.not. (l_qr_acr_qg .and. l_qr_acr_qs .and. l_qi_aut_qs .and. l_freezeH2O)) then
       write(mpas_err_message,'(A)') &
          '--- tables to run the Thompson cloud microphysics scheme do not exist: run build_tables first.'
       call physics_error_fatal(mpas_err_message)
    endif
!   call mpas_log_write('l_mp_tables = $l',logicArgs=(/l_mp_tables/))

 endif

 end subroutine physics_tables_init

!=================================================================================================================
!  routine physics_compatibility_check()
!
!> \brief Checks physics input fields and options for compatibility
!> \author Miles Curry and Michael Duda
!> \date   25 October 2018
!> \details
!>  This routine checks the input fields and run-time options provided
!>  by the user for compatibility. For example, two run-time options may
!>  be mutually exclusive, or an option may require that a certain input
!>  field is provided. The checks performed by this routine are only for
!>  physics related fields and options.
!>
!>  A value of 0 is returned if there are no incompatibilities among
!>  the provided input fields and run-time options, and a non-zero value
!>  otherwise.
!>
 subroutine physics_compatibility_check(dminfo, blockList, streamManager, ierr)
!=================================================================================================================

 implicit none

 type (dm_info), pointer :: dminfo
 type (block_type), pointer :: blockList
 type (MPAS_streamManager_type), pointer :: streamManager
 integer, intent(out) :: ierr

 real (kind=RKIND) :: maxvar2d_local, maxvar2d_global
 real (kind=RKIND), dimension(:), pointer :: var2d
 integer, pointer :: nCellsSolve
 integer, pointer :: iswater_lu
 integer, pointer, dimension(:) :: ivgtyp
 integer :: all_water, iall_water
 character (len=StrKIND), pointer :: gwdo_scheme

 ! Prognostic ozone checks. The state scalar field is used for both
 ! cold-start initial conditions and restart states.
 character (len=StrKIND), pointer :: o3_source
 real (kind=RKIND) :: maxqo3_local, maxqo3_global
 real (kind=RKIND), dimension(:,:), pointer :: qo3
 real (kind=RKIND), dimension(:,:,:), pointer :: scalars
 integer, pointer :: nVertLevels
 integer, pointer :: index_qo3

 type (block_type), pointer :: block
 type (mpas_pool_type), pointer :: meshPool
 type (mpas_pool_type), pointer :: sfc_inputPool
 type (mpas_pool_type), pointer :: statePool

 ierr = 0

 call mpas_pool_get_config(blocklist % configs, 'config_gwdo_scheme', gwdo_scheme)

 if (trim(gwdo_scheme) == 'bl_ysu_gwdo') then
     maxvar2d_local = -huge(maxvar2d_local)
     block => blockList
     do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'sfc_input', sfc_inputPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(sfc_inputPool, 'var2d', var2d)

         maxvar2d_local = max(maxvar2d_local, maxval(var2d(1:nCellsSolve)))

         block => block % next
     end do

     call mpas_dmpar_max_real(dminfo, maxvar2d_local, maxvar2d_global)

     !
     ! The GWDO check below can fail on regional simulations that are completely above
     ! water. So, check to see if the simulation is completely above water and do not
     ! throw the error if it is.
     !
     call mpas_pool_get_array(sfc_inputPool, 'iswater', iswater_lu)
     call mpas_pool_get_array(sfc_inputPool, 'ivgtyp', ivgtyp)
     if (all(ivgtyp(1:nCellsSolve) == iswater_lu)) then
         all_water = 1 ! All water
     else
         all_water = 0 ! Land present
     end if

     call mpas_dmpar_min_int(dminfo, all_water, iall_water)

     if (maxvar2d_global <= 0.0_RKIND .and. iall_water /= 1) then
         call mpas_log_write('*******************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('The YSU GWDO scheme requires valid var2d, con, oa{1,2,3,4}, and ol{1,2,3,4} fields,', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('but these fields appear to be zero everywhere in the model input.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('Either set config_gwdo_scheme = ''off'' in the &physics namelist, or generate', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('the GWDO static fields with the init_atmosphere core.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('*******************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         ierr = ierr + 1
     end if

 end if

 !
 ! When config_o3_source = 'prognostic', qo3 is used directly by RRTMG.
 ! After initial-condition/restart loading, reject a state in which qo3 is effectively
 ! zero everywhere. This catches missing/uninitialized O3MR in an initial condition
 ! or missing qo3 in a restart file before radiation runs without ozone.
 ! Users must provide valid prognostic ozone or use config_o3_source = 'climatology'.
 !
 call mpas_pool_get_config(blockList % configs, 'config_o3_source', o3_source)

 if (trim(o3_source) == 'prognostic') then
     maxqo3_local = 0.0_RKIND

     block => blockList
     do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(statePool, 'index_qo3', index_qo3)
         call mpas_pool_get_array(statePool, 'scalars', scalars, 1)

         qo3 => scalars(index_qo3,:,:)

         maxqo3_local = max(maxqo3_local, &
                            maxval(abs(qo3(1:nVertLevels,1:nCellsSolve))))

         block => block % next
     end do

     call mpas_dmpar_max_real(dminfo, maxqo3_local, maxqo3_global)

     if (maxqo3_global <= 1.0e-20_RKIND) then
         call mpas_log_write('')
         call mpas_log_write('*******************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('The option config_o3_source = ''prognostic'' requires a valid non-zero ozone', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('mixing-ratio field (qo3), but qo3 is zero everywhere in the model input.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('This may indicate that O3MR was missing or was not correctly processed when', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('the MPAS initial-condition or restart file was generated.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('Either provide valid prognostic ozone initial conditions or a restart file,', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('or set config_o3_source = ''climatology'' to use monthly climatology for ozone', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('To include qo3 in the initial-condition file, set config_use_ozone=true when running init_atmosphere.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('Also, verify that the intermediate file used by init_atmosphere contains the O3MR field.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('*******************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('')
         ierr = ierr + 1
     end if

 end if

 end subroutine physics_compatibility_check

!=================================================================================================================
 end module mpas_atmphys_control
!=================================================================================================================

