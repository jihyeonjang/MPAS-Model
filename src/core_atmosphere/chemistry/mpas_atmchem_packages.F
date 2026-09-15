! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmchem_packages
 use mpas_kind_types
 use mpas_derived_types,only : mpas_pool_type,mpas_io_context_type,MPAS_LOG_ERR
 use mpas_pool_routines,only : mpas_pool_get_config,mpas_pool_get_package
 use mpas_log,only : mpas_log_write

 implicit none
 private
 public:: chemistry_setup_packages

!mpas_atmchem_packages contains the definitions of all chemistry packages.
!Laura D. Fowler (laura@ucar.edu) / 2022-07-27.


 contains


!=================================================================================================================
 function chemistry_setup_packages(configs,packages,iocontext) result(ierr)
!=================================================================================================================

!inout arguments:
 type (mpas_pool_type), intent(inout) :: configs
 type (mpas_pool_type), intent(inout) :: packages
 type (mpas_io_context_type), intent(inout) :: iocontext

!local variables:
 character(len=StrKIND),pointer:: o3_source,o3_scheme,h2o_scheme
 logical,pointer:: o3_prognostic_in,opp_in,photochem_tendencies_in

 integer :: ierr

!-----------------------------------------------------------------------------------------------------------------
 ierr = 0
 call mpas_pool_get_config(configs,'config_o3_source',o3_source)
 call mpas_pool_get_config(configs,'config_photochem_o3_scheme',o3_scheme)
 call mpas_pool_get_config(configs,'config_photochem_h2o_scheme',h2o_scheme)
 call mpas_pool_get_package(packages,'o3_prognostic_inActive',o3_prognostic_in)
 call mpas_pool_get_package(packages,'opp_inActive',opp_in)
 call mpas_pool_get_package(packages,'photochem_tendencies_inActive',photochem_tendencies_in)
 if(.not.associated(o3_prognostic_in).or..not.associated(opp_in).or. &
    .not.associated(photochem_tendencies_in)) then
    call mpas_log_write('Unable to obtain atmosphere chemistry package controls',messageType=MPAS_LOG_ERR)
    ierr=1
    return
 endif
 o3_prognostic_in=trim(o3_source)=='prognostic'
 opp_in=trim(o3_scheme)=='opp'
 photochem_tendencies_in=trim(o3_scheme)/='off'.or.trim(h2o_scheme)/='off'
 call mpas_log_write('    o3_prognostic_in        = $l',logicArgs=(/o3_prognostic_in/))
 call mpas_log_write('    opp_in                  = $l',logicArgs=(/opp_in/))
 call mpas_log_write('    photochem_tendencies_in = $l',logicArgs=(/photochem_tendencies_in/))

! call mpas_log_write('--- enter subroutine chemistry_setup_packages:')
! call mpas_log_write('')

 end function chemistry_setup_packages

!=================================================================================================================
 end module mpas_atmchem_packages
!=================================================================================================================



