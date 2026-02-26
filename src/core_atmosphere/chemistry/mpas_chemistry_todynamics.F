! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_chemistry_todynamics
 use mpas_log,only: mpas_log_write
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_dimension,mpas_pool_get_array
 use mpas_dmpar
 use mpas_atm_dimensions

 implicit none
 private
 public:: chemistry_get_tend


!computation of dynamics tendency for use in the dynamical core.

!add-ons and modifications to sourcecode:
!----------------------------------------
! * in subroutine chemistry_get_tend, hit return if config_o3prognostic is set to false.
!   Laura D. Fowler (laura@ucar.edu) / 2021-04-08.


 contains

 
!=================================================================================================================
 subroutine chemistry_get_tend(configs,mesh,state,tend,tend_chemistry)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: tend_chemistry

!inout variables:
 type(mpas_pool_type),intent(inout):: tend

!local variables:
 logical,pointer:: o3prognostic

 integer:: iCell,k
 integer,pointer:: nCells,nCellsSolve
 integer,pointer:: index_qo3

 real(kind=RKIND),dimension(:,:),pointer:: mass      ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:),pointer:: qo3oppten
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars

!-----------------------------------------------------------------------------------------------------------------

!--- hit return if config_o3prognostic is set to false:
 call mpas_pool_get_config(configs,'config_o3prognostic',o3prognostic)
 if(.not. o3prognostic) return

 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_get_tend:')

 call mpas_pool_get_dimension(state,'index_qo3'  ,index_qo3  )

 call mpas_pool_get_array(tend,'scalars_tend',tend_scalars)
!JHJ already applied in physics_get_tend tend_scalars(index_qo3,:,:) = 0._RKIND

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(state,'rho_zz',mass,2)

 call mpas_pool_get_array(tend_chemistry,'qo3oppten',qo3oppten)

!
!in case some variables are not allocated due to their associated packages, we need to make their pointers
!associated here to avoid triggering run-time checks when calling chemistry_get_tend_work:
!
 if (.not.associated(qo3oppten)) allocate(qo3oppten(0,0))

 call chemistry_get_tend_work(nCells,nCellsSolve,o3prognostic,index_qo3,mass, &
                              qo3oppten,tend_scalars)

!
!clean up any pointers that were allocated with zero size before the call to chemistry_get_tend_work:
!
 if (size(qo3oppten) == 0) deallocate(qo3oppten)

 call mpas_log_write('--- end subroutine chemistry_get_tend:')

 end subroutine chemistry_get_tend

!=================================================================================================================
 subroutine chemistry_get_tend_work(nCells,nCellsSolve,o3prognostic,index_qo3,mass, &
                                    qo3oppten,tend_scalars)
!=================================================================================================================

!input arguments:
 logical,intent(in) :: o3prognostic

 integer,intent(in):: nCells,nCellsSolve
 integer,intent(in):: index_qo3
 
 real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: mass
 real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: qo3oppten

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(num_scalars,nVertLevels,nCells+1):: tend_scalars

!local variables and arrays:
 integer:: i,k

!-----------------------------------------------------------------------------------------------------------------

!--- add coupled tendency due to CHEM2D Ozone Photochemistry Parameterization (CHEM2D-OPP):
 if(o3prognostic) then

    do i = 1,nCellsSolve
    do k = 1,nVertLevels
       tend_scalars(index_qo3,k,i) = tend_scalars(index_qo3,k,i) + qo3oppten(k,i)*mass(k,i)
    enddo
    enddo

 endif
 
 end subroutine chemistry_get_tend_work

!=================================================================================================================
 end module mpas_chemistry_todynamics
!=================================================================================================================
