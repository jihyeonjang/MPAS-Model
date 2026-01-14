! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_chemistry_packages
 use mpas_kind_types
 use mpas_derived_types,only : mpas_pool_type,mpas_io_context_type,MPAS_LOG_ERR
 use mpas_pool_routines,only : mpas_pool_get_config,mpas_pool_get_package
 use mpas_log,only : mpas_log_write

 implicit none
 private
 public:: chemistry_setup_packages

!mpas_chemistry_packages contains the definitions of all chemistry packages.
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
 logical,pointer:: o3prognostic
 logical,pointer:: chem2dopp_in

 integer :: ierr

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine chemistry_setup_packages:')

 ierr = 0

 call mpas_log_write('----- Setting up package variables -----')
 call mpas_log_write('')

!--- initialization of all packages for parameterizations of cloud microphysics:

 call mpas_pool_get_config(configs,'config_o3prognostic',o3prognostic)

 nullify(chem2dopp_in)
 call mpas_pool_get_package(packages,'chem2dopp_inActive',chem2dopp_in)

 if(.not.associated(chem2dopp_in)) then
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    call mpas_log_write('* Error while setting up packages chemistry options in atmosphere core.'             ,messageType=MPAS_LOG_ERR)
    call mpas_log_write('====================================================================================',messageType=MPAS_LOG_ERR)
    ierr = 1
    return
 endif

 chem2dopp_in   = .false.

 if(o3prognostic) then
    chem2dopp_in = .true.
 endif

 call mpas_log_write('    chem2dopp_in            = $l', logicArgs=(/chem2dopp_in/))

!call mpas_log_write('--- enter subroutine chemistry_setup_packages:')
 call mpas_log_write('')

 end function chemistry_setup_packages

!=================================================================================================================
 end module mpas_chemistry_packages
!=================================================================================================================



