!=================================================================================================================
 module mpas_chemistry_interface
 use mpas_log
 use mpas_kind_types
 use mpas_atmphys_constants,only: R_v,R_d
 use mpas_constants,only: gravity
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension

 implicit none

 type,public:: chems_vars
    integer:: its,ite,kts,kte,ktep1
    integer:: im,levs,no3,no3levs

    !--- ozone related fields:
    real(kind=RKIND),dimension(:),pointer    :: qozp      => null()
    real(kind=RKIND),dimension(:,:),pointer  :: qoz       => null()
    real(kind=RKIND),dimension(:,:),pointer  :: tend_qoz  => null()
    real(kind=RKIND),dimension(:,:),pointer  :: qozd1     => null()
    real(kind=RKIND),dimension(:,:),pointer  :: qozd2     => null()
    real(kind=RKIND),dimension(:,:),pointer  :: qozd3     => null()
    real(kind=RKIND),dimension(:,:),pointer  :: qozd4     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qozcoeffs => null()

    !--- local variables from the MPAS dynamical core:
    real(kind=RKIND):: dt_c
    real(kind=RKIND):: grav_c
    real(kind=RKIND),dimension(:,:),pointer:: temp_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: pres_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: dpres_c     => null()


    contains
       procedure:: chemistry_dims       => mpas_chemistry_dims
       procedure:: chemistry_allocate   => mpas_chemistry_allocate
       procedure:: chemistry_deallocate => mpas_chemistry_deallocate
       procedure:: from_MPAS            => mpas_chemistry_from_MPAS
       procedure:: to_MPAS              => mpas_chemistry_to_MPAS
 end type


 contains


!=================================================================================================================
 subroutine mpas_chemistry_dims(this,mesh)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 class(chems_vars),intent(inout):: this

!local variables:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: nOPPLevels,nOPPCoeffs

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_dims:')

 call mpas_pool_get_dimension(mesh,'nOPPLevels' ,nOPPLevels  )
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs' ,nOPPCoeffs  )
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels )

 this%its = 1 ; this%ite = nCellsSolve
 this%kts = 1 ; this%kte = nVertLevels ; this%ktep1 = nVertLevels+1
 call mpas_log_write('ITS = $i   ITE = $i',intArgs=(/this%its,this%ite/))
 call mpas_log_write('KTS = $i   KTE = $i',intArgs=(/this%kts,this%kte/))

 this%im   = this%ite-this%its+1
 this%levs = this%kte-this%kts+1

 this%no3     = nOPPCoeffs
 this%no3levs = nOPPLevels

 call mpas_log_write('--- end subroutine mpas_chemistry_dims:')

 end subroutine mpas_chemistry_dims

!=================================================================================================================
 subroutine mpas_chemistry_allocate(this)
!=================================================================================================================

!inout arguments:
 class(chems_vars),intent(inout):: this

!local variables and arrays:
 integer:: its,ite,kts,kte,ktep1
 integer:: no3,no3levs

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_allocate:')

 its     = this%its
 ite     = this%ite
 kts     = this%kts
 kte     = this%kte
 ktep1   = this%ktep1
 no3     = this%no3
 no3levs = this%no3levs


!--- allocate ozone-related fields:
 if(.not.associated(this%qoz)      ) allocate(this%qoz(its:ite,kts:kte)          )
 if(.not.associated(this%qozd1)    ) allocate(this%qozd1(its:ite,kts:kte)        )
 if(.not.associated(this%qozd2)    ) allocate(this%qozd2(its:ite,kts:kte)        )
 if(.not.associated(this%qozd3)    ) allocate(this%qozd3(its:ite,kts:kte)        )
 if(.not.associated(this%qozd4)    ) allocate(this%qozd4(its:ite,kts:kte)        )
 if(.not.associated(this%tend_qoz) ) allocate(this%tend_qoz(its:ite,kts:kte)     )

 if(.not.associated(this%qozp)     ) allocate(this%qozp(no3levs))
 if(.not.associated(this%qozcoeffs)) allocate(this%qozcoeffs(its:ite,no3levs,no3))


!--- allocate meterological fields:
 if(.not.associated(this%temp_c)   ) allocate(this%temp_c(its:ite,kts:kte)       )
 if(.not.associated(this%pres_c)   ) allocate(this%pres_c(its:ite,kts:kte)       )
 if(.not.associated(this%dpres_c)  ) allocate(this%dpres_c(its:ite,kts:kte)      )

 call mpas_log_write('--- end subroutine mpas_chemistry_allocate:')

 end subroutine mpas_chemistry_allocate

!=================================================================================================================
 subroutine mpas_chemistry_deallocate(this)
!=================================================================================================================

!inout arguments:
 class(chems_vars),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_deallocate:')

!--- deallocate ozone-related fields:
 if(associated(this%qoz)      ) deallocate(this%qoz      )
 if(associated(this%qozd1)    ) deallocate(this%qozd1    )
 if(associated(this%qozd2)    ) deallocate(this%qozd2    )
 if(associated(this%qozd3)    ) deallocate(this%qozd3    )
 if(associated(this%qozd4)    ) deallocate(this%qozd4    )
 if(associated(this%tend_qoz) ) deallocate(this%tend_qoz )

 if(associated(this%qozp)     ) deallocate(this%qozp     )
 if(associated(this%qozcoeffs)) deallocate(this%qozcoeffs)


!--- deallocate meterological fields:
 if(associated(this%temp_c)  ) deallocate(this%temp_c    )
 if(associated(this%pres_c)  ) deallocate(this%pres_c    )
 if(associated(this%dpres_c) ) deallocate(this%dpres_c   )

 call mpas_log_write('--- end subroutine mpas_chemistry_deallocate:')

 end subroutine mpas_chemistry_deallocate

!=================================================================================================================
 subroutine mpas_chemistry_from_MPAS(this,configs,mesh,chem_input,diag,state,time_lev)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: chem_input
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: state

 integer,intent(in):: time_lev


!inout arguments:
 class(chems_vars),intent(inout):: this


!local variables and arrays:
 character(len=StrKIND):: errmsg

 integer:: i,its,ite,k,kts,kte,n
 integer:: no3,no3levs
 integer:: errflg
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: index_qo3,index_qv

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: qo3,qv,dpres,pres
 real(kind=RKIND),dimension(:,:),pointer  :: exner,pressure_b,pressure_p,theta_m
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND):: fzm,fzp,tem,z0,z1,z2
 real(kind=RKIND),dimension(:,:),allocatable:: pres2

!local dimensions and arrays specific to CHEM2D-OPP:
 integer,pointer:: nOPPLevels,nOPPCoeffs
 real(kind=RKIND),dimension(:),pointer    :: OPPpres
 real(kind=RKIND),dimension(:,:,:),pointer:: OPPcoeffs

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_from_MPAS:')


!--- dimensions and arrays specific to CHEM2D-OPP:
 call mpas_pool_get_dimension(mesh,'nOPPLevels',nOPPLevels)
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs',nOPPCoeffs)

 call mpas_pool_get_array(chem_input,'OPPpres'  ,OPPpres  )
 call mpas_pool_get_array(chem_input,'OPPcoeffs',OPPcoeffs)


!--- dimensions and arrays for mesh and meteorological fields:
 call mpas_pool_get_config(configs,'config_dt',dt)

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'zgrid',zgrid)

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)

 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)

 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)
 call mpas_pool_get_dimension(state,'index_qv' ,index_qv )
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qo3 => scalars(index_qo3,:,:)
 qv  => scalars(index_qv,:,: )


 its     = this%its
 ite     = this%ite
 kts     = this%kts
 kte     = this%kte
 no3     = this%no3
 no3levs = this%no3levs


!--- initializes input CHEM2D-OPP data:
 do k = 1,no3levs
    this%qozp(k) = OPPpres(k)
 enddo
 do n = 1,no3
    do k = 1,no3levs
       do i = its,ite
          this%qozcoeffs(i,k,n) = OPPcoeffs(n,k,i)
       enddo
    enddo
 enddo


!--- initializes local variables and arrays before call to o3_gfs_run:
 this%dt_c   = dt
 this%grav_c = gravity

 do k = kts,kte
    do i = its,ite
       this%temp_c(i,k)   = theta_m(k,i)/(1._RKIND + R_v/R_d*qv(k,i))
       this%temp_c(i,k)   = this%temp_c(i,k)*exner(k,i)
       this%pres_c(i,k)   = pressure_p(k,i) + pressure_b(k,i)
       this%qoz(i,k)      = qo3(k,i)

       !initialization of tendency:
       this%tend_qoz(i,k) = qo3(k,i)

       !output diagnostics:
       this%qozd1(i,k)    = 0._RKIND
       this%qozd2(i,k)    = 0._RKIND
       this%qozd3(i,k)    = 0._RKIND
       this%qozd4(i,k)    = 0._RKIND
    enddo
 enddo


!--- interpolation of pressure from theta points to w points:
 pres  => this%pres_c
 dpres => this%dpres_c

 if(.not.allocated(pres2)) allocate(pres2(its:ite,kts:kte+1))
 do i = its,ite

    !pressure at the surface:
    k = kts
    z0  = zgrid(k,i)
    z1  = 0.5*(zgrid(k,i)+zgrid(k+1,i))
    z2  = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
    fzp = (z0-z2)/(z1-z2)
    fzm = 1.-fzp
    pres2(i,k) = fzp*pres(i,k) + fzm*pres(i,k+1)
    !pressure at interface between layers above the surface (k=kts) and below the model-top (k=kte+1):
    do k = kts+1,kte
       tem = 1./(zgrid(k+1,i)-zgrid(k-1,i))
       fzm = (zgrid(k,i)-zgrid(k-1,i))*tem
       fzp = (zgrid(k+1,i)-zgrid(k,i))*tem
       pres2(i,k) = fzm*pres(i,k) + fzp*pres(i,k-1)
    enddo
    !pressure at the model-top (k=kte+1):
    k = kte+1
    z0  = zgrid(k,i)
    z1  = 0.5*(zgrid(k,i)+zgrid(k-1,i))
    z2  = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
    fzm = (z0-z2)/(z1-z2)
    fzp = 1.-fzm
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pres2(i,k) = exp(fzm*log(pres(i,k-1)) + fzp*log(pres(i,k-2)))

    !pressure thickness:
    do k = kts,kte
       dpres(i,k) = pres2(i,k)-pres2(i,k+1)
    enddo
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r',intArgs=(/i,k/),realArgs=(/pres2(i,k),pres(i,k),dpres(i,k)/))
!   enddo
!   k = kte+1
!   call mpas_log_write('$i $i $r',intArgs=(/i,k/),realArgs=(/pres2(i,k)/))
!   call mpas_log_write(' ')
 enddo
 if(allocated(pres2)) deallocate(pres2)

 call mpas_log_write('--- end subroutine mpas_chemistry_from_MPAS:')

 end subroutine mpas_chemistry_from_MPAS

!=================================================================================================================
 subroutine mpas_chemistry_to_MPAS(this,mesh,diag_chemistry,tend_chemistry)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh


!inout arguments:
 class(chems_vars),intent(inout):: this
 type(mpas_pool_type),intent(inout):: diag_chemistry
 type(mpas_pool_type),intent(inout):: tend_chemistry
 

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,n
 integer:: no3,no3levs
 integer,pointer:: nCells,nVertLevels

 real(kind=RKIND):: dt_c
 real(kind=RKIND),dimension(:,:),pointer:: qoz,tend_qoz,qozd1,qozd2,qozd3,qozd4
 real(kind=RKIND),dimension(:,:),pointer:: oppdiag1,oppdiag2,oppdiag3,oppdiag4
 real(kind=RKIND),dimension(:,:),pointer:: qo3oppten

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_to_MPAS:')

 call mpas_pool_get_array(diag_chemistry,'oppdiag1',oppdiag1)
 call mpas_pool_get_array(diag_chemistry,'oppdiag2',oppdiag2)
 call mpas_pool_get_array(diag_chemistry,'oppdiag3',oppdiag3)
 call mpas_pool_get_array(diag_chemistry,'oppdiag4',oppdiag4)

 call mpas_pool_get_array(tend_chemistry,'qo3oppten',qo3oppten)

 its     = this%its
 ite     = this%ite
 kts     = this%kts
 kte     = this%kte
 no3     = this%no3
 no3levs = this%no3levs

!--- updates the local tendency and diagnostics:
 dt_c = this%dt_c

 qoz      => this%qoz
 tend_qoz => this%tend_qoz
 qozd1    => this%qozd1
 qozd2    => this%qozd2
 qozd3    => this%qozd3
 qozd4    => this%qozd4

 do i = its,ite
    do k = kts,kte
       tend_qoz(i,k) = (qoz(i,k)-tend_qoz(i,k))/dt_c
    enddo
 enddo


!---
 do i = its,ite
    do k = kts,kte
       oppdiag1(k,i) = qozd1(i,k)
       oppdiag2(k,i) = qozd2(i,k)
       oppdiag3(k,i) = qozd3(i,k)
       oppdiag4(k,i) = qozd4(i,k)
       qo3oppten(k,i) = tend_qoz(i,k)
    enddo
 enddo

 call mpas_log_write('--- end subroutine mpas_chemistry_to_MPAS:')

 end subroutine mpas_chemistry_to_MPAS

!=================================================================================================================
 end module mpas_chemistry_interface
!=================================================================================================================

