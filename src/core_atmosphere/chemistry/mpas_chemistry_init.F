!=================================================================================================================
 module mpas_chemistry_init
 use mpas_kind_types
 use mpas_log,only          : mpas_log_write
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config

 use mpas_chemistry_opp,only: init_opp


 implicit none
 private
 public:: chemistry_init


!initialization and interpolation to MPAS mesh of input data needed to run the CHEM2D Ozone Photochemistry
!Parameterization (CHEM2D-OPP) of McCormack et al. (2006).

!add-ons and modifications to sourcecode:
!----------------------------------------
! * in subroutine chemistry_init, hit return if config_o3prognostic is set to false.
!   Laura D. Fowler (laura@ucar.edu) / 2021-04-08.


 contains


!=================================================================================================================
 subroutine chemistry_init(configs,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh


!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input


!local variables and pointers:
 logical,pointer:: o3prognostic

!-----------------------------------------------------------------------------------------------------------------

!--- hit return if config_o3prognostic is set to false:
 call mpas_pool_get_config(configs,'config_o3prognostic',o3prognostic)
 if(.not. o3prognostic) return

 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_init:')

!--- initializes the input data need to run the CHEM2D-OPP of McCormack et al. (2006):
 call init_opp(configs,mesh,chem_input)


 call mpas_log_write('--- end subroutine chemistry_init:')
 call mpas_log_write(' ')

 end subroutine chemistry_init

!=================================================================================================================
 end module mpas_chemistry_init
!=================================================================================================================
