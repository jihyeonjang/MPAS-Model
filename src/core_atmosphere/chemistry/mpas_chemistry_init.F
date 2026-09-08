!=================================================================================================================
 module mpas_chemistry_init
 use mpas_kind_types
 use mpas_log, only : mpas_log_write
 use mpas_derived_types, only : mpas_pool_type, MPAS_LOG_CRIT
 use mpas_pool_routines, only : mpas_pool_get_config
 use mpas_timekeeping

 use mpas_chemistry_opp, only : init_opp

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
 subroutine chemistry_init(configs,clock,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(MPAS_Clock_type):: clock
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input

!local variables and pointers:
 character(len=*), parameter:: oppAlarmID = 'opp'

 character(len=StrKIND),pointer:: OPP_update
 logical,pointer:: o3prognostic

 type(MPAS_Time_Type):: startTime,alarmStartTime
 type(MPAS_TimeInterval_Type):: alarmTimeStep

 integer:: ierr

!-----------------------------------------------------------------------------------------------------------------

!--- hit return if config_o3prognostic is set to false:
 call mpas_pool_get_config(configs,'config_o3prognostic',o3prognostic)
! if(.not. o3prognostic) return

 call mpas_pool_get_config(configs,'config_OPP_update',OPP_update)

! call mpas_log_write(' ')
! call mpas_log_write('--- enter subroutine chemistry_init:')

 if (o3prognostic) then
    !--- get startTime prior to initializating alarms:
    startTime = mpas_get_clock_time(clock,MPAS_START_TIME,ierr)

    !--- set alarm for updating the coefficients in the GFS linearized CHEM2D-OPP:
    call mpas_set_timeInterval(alarmTimeStep,timeString=OPP_update,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,oppAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) then
!       call physics_error_fatal('subroutine chemistry_init: error creating alarm GFS CHEM2D-OPP')
       call mpas_log_write('subroutine chemistry_init: error creating alarm GFS CHEM2D-OPP', &
                            messageType=MPAS_LOG_CRIT)
    end if

    !--- initializes the input data need to run the CHEM2D-OPP of McCormack et al. (2006):
    call init_opp(configs,mesh,chem_input)
 end if

! call mpas_log_write('--- end subroutine chemistry_init:')

 end subroutine chemistry_init

!=================================================================================================================
 end module mpas_chemistry_init
!=================================================================================================================
