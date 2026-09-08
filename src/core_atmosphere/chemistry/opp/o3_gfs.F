!=================================================================================================================
 module o3_gfs
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log,only: mpas_log_write
 use mpas_io_units,only: mpas_new_unit

 use mpas_atmphys_utilities,only: physics_error_fatal


 implicit none
 private
 public:: read_o3data,o3_gfs_run

 integer,public:: latsozp,levozp,oz_coeff,timeoz
 real(kind=kind_phys),dimension(:),allocatable,public:: oz_lat,oz_pres,oz_time
 real(kind=kind_phys),dimension(:,:,:,:),allocatable,public:: ozplin


 contains


!=================================================================================================================
 subroutine read_o3data()
!=================================================================================================================

!local variables and arrays:
 integer:: i,k,ll,n,o3unit
 integer:: latsozc,levozc,timeozc

 real(kind=4),dimension(:),allocatable:: oz_lat4,oz_pres4,oz_time4,tempin

 integer,dimension(:),allocatable:: coeff_min_loc,coeff_max_loc
 real(kind=kind_phys):: tau_min,tau_max
 real(kind=kind_phys),dimension(:),allocatable:: coeff_min,coeff_max

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine read_o3data:')

!--- define unit:
 call mpas_new_unit(o3unit)
 if(o3unit < 0) &
    call physics_error_fatal('GFSprogO3 for MPAS: All file units are taken. Change maxUnits in mpas_io_units.F')


!--- prognostic ozone input data file:
 open(unit=o3unit,file='global_o3prdlos.f77',form='unformatted')


!--- read in indices:
 read (o3unit) oz_coeff,latsozp,levozp,timeoz
 call mpas_log_write('reading in o3data from global_o3prdlos.f77: ')
 call mpas_log_write('oz_coeff = $i',intArgs=(/oz_coeff/))
 call mpas_log_write('latsozp  = $i',intArgs=(/latsozp/))
 call mpas_log_write('levozp   = $i',intArgs=(/levozp/))
 call mpas_log_write('timeoz   = $i',intArgs=(/timeoz/))


!--- read in data
!---   oz_lat   -  latitude of data        (-90 to 90)
!---   oz_pres  -  vertical pressure level (mb)
!---   oz_time  -  time coordinate         (days)
!---

 if(.not.allocated(oz_lat) ) allocate(oz_lat(latsozp)    )
 if(.not.allocated(oz_pres)) allocate(oz_pres(levozp)    )
 if(.not.allocated(oz_time)) allocate(oz_time(timeoz+1)  )

 if(.not.allocated(oz_lat4) ) allocate(oz_lat4(latsozp)  )
 if(.not.allocated(oz_pres4)) allocate(oz_pres4(levozp)  )
 if(.not.allocated(oz_time4)) allocate(oz_time4(timeoz+1))
 if(.not.allocated(tempin)  ) allocate(tempin(latsozp)   )

 rewind(o3unit)
 read(o3unit) oz_coeff,latsozp,levozp,timeoz,oz_lat4,oz_pres4,oz_time4

!--- convert pressure levels from mb to ln(Pa):
 oz_pres(:) = oz_pres4(:)
 oz_pres(:) = log(100.0*oz_pres(:))
 oz_lat(:)  = oz_lat4(:)
 oz_time(:) = oz_time4(:)

call mpas_log_write(' ')
call mpas_log_write('--- pressure:')
do k = 1,levozp
   call mpas_log_write('$i $r $r',intArgs=(/k/),realArgs=(/oz_pres4(k),oz_pres(k)/))
enddo
call mpas_log_write(' ')
call mpas_log_write('--- latitude:')
do k = 1,latsozp
   call mpas_log_write('$i $r',intArgs=(/k/),realArgs=(/oz_lat(k)/))
enddo
call mpas_log_write(' ')
call mpas_log_write('--- oz_time:')
do k = 1,timeoz+1
   call mpas_log_write('$i $r',intArgs=(/k/),realArgs=(/oz_time(k)/))
enddo

!--- read in ozplin which is in order of (latitudes,ozone levels,coeff number,time) where latitudes are on
!    a uniform gaussian grid:
 if(.not.allocated(ozplin)) allocate(ozplin(latsozp,levozp,oz_coeff,timeoz))

 do i = 1,timeoz
    do n = 1,oz_coeff
       do k = 1,levozp
          read(o3unit) tempin
          ozplin(:,k,n,i) = tempin(:)
       enddo
    enddo
 enddo

!if(.not.allocated(coeff_min_loc)) allocate(coeff_min_loc(oz_coeff))
!if(.not.allocated(coeff_max_loc)) allocate(coeff_max_loc(oz_coeff))
!if(.not.allocated(coeff_min)    ) allocate(coeff_min(oz_coeff)    )
!if(.not.allocated(coeff_max)    ) allocate(coeff_max(oz_coeff)    )
!i = 7
!do ll = 1,latsozp
!   do n =1,oz_coeff
!      coeff_min(n)     = minval(ozplin(ll,:,n,i))
!      coeff_max(n)     = maxval(ozplin(ll,:,n,i))
!   enddo
!   tau_min = 0.
!   tau_max = 0.
!   if(abs(coeff_min(2)) .gt. 0_kind_phys) tau_min = -1./coeff_min(2)
!   if(abs(coeff_max(2)) .gt. 0_kind_phys) tau_max = -1./coeff_max(2)
!   tau_min = tau_min / 86400.
!   tau_max = tau_max / 86400.
!   call mpas_log_write('$i $r $r $r $r $r',intArgs=(/ll/),realArgs=(/oz_lat(ll), &
!                   (coeff_min(n),n=2,2),(coeff_max(n),n=2,2),tau_min,tau_max/))
!enddo
!if(allocated(coeff_min_loc)) deallocate(coeff_min_loc)
!if(allocated(coeff_max_loc)) deallocate(coeff_max_loc)
!if(allocated(coeff_min)    ) deallocate(coeff_min    )
!if(allocated(coeff_max)    ) deallocate(coeff_max    )

 if(allocated(oz_lat4) ) deallocate(oz_lat4 )
 if(allocated(oz_pres4)) deallocate(oz_pres4)
 if(allocated(oz_time4)) deallocate(oz_time4)
 if(allocated(tempin)  ) deallocate(tempin  )


 close(o3unit)
 call mpas_log_write('--- end subroutine read_o3data:')

 end subroutine read_o3data

!=================================================================================================================
 subroutine o3_gfs_run(im,levs,ko3,dt,oz,tin,po3,prsl,prdout,pl_coeff,delp,ldiag3d,qdiag3d, &
                       ozp1,ozp2,ozp3,ozp4,grav,errmsg,errflg)
!=================================================================================================================

!input arguments:
 logical,intent(in):: &
    ldiag3d, &! flag for calculating 3-D diagnostic fields.
    qdiag3d   ! flag for 3d tracer diagnostic fields.

 integer,intent(in):: &
    im,      &!
    levs,    &!
    ko3,     &!
    pl_coeff  !

 real(kind=kind_phys),intent(in):: &
    dt,      &! time-step (s).
    grav      ! gravitational acceleration (m s-2).

 real(kind=kind_phys),intent(in),dimension(im,levs):: &
    delp,    &! layer pressure thickness (Pa).
    prsl,    &! air pressure (Pa).
    tin       ! air temperature (K).

 real(kind=kind_phys),intent(in),dimension(ko3):: &
    po3       ! natural log of ozone forcing data pressure levels (ln(Pa))

 real(kind=kind_phys),intent(in),dimension(im,ko3,pl_coeff):: &
    prdout    ! ozone forcing data (various).


!inout arguments:
 real(kind=kind_phys),intent(inout),dimension(im,levs):: &
    ozp1,    &! cumulative change in ozone concentration due to production and loss rate (kg kg-1).
    ozp2,    &! cumulative change in ozone concentration due to ozone mixing ratio (kg kg-1).
    ozp3,    &! cumulative change in ozone concentration due to temperature (kg kg-1).
    ozp4      ! cumulative change in ozone concentration due to overhead ozone column (kg kg-1).

 real(kind=kind_phys),intent(inout),dimension(im,levs):: &
    oz        ! ozone mixing ratio (kg/kg).


!output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg


!local variables and arrays:
 logical,dimension(im):: flg

 integer:: i,j,k,l
 integer:: kmin,kmax

 real(kind=kind_phys):: gravi
 real(kind=kind_phys):: pmax,pmin,tem,temp
 real(kind=kind_phys),dimension(im)         :: ozib,wk1,wk2,wk3
 real(kind=kind_phys),dimension(im,levs)    :: ozi
 real(kind=kind_phys),dimension(im,levs+1)  :: colo3,coloz
 real(kind=kind_phys),dimension(im,pl_coeff):: prod

!-----------------------------------------------------------------------------------------------------------------
! call mpas_log_write(' ')
! call mpas_log_write('--- enter subroutine o3_gfs_run:')

 gravi = 1./grav

!--- initialize CCPP error handling variables:
 errmsg = ' '
 errflg = 0

!--- save input ozone profiles:
 do k = 1,levs
    do i = 1,im
       ozi(i,k)  = oz(i,k)
       ozp1(i,k) = 0.
       ozp2(i,k) = 0.
       ozp3(i,k) = 0.
       ozp4(i,k) = 0.
    enddo
 enddo

!--- initialization of the overlying vertically-integrated mixing ratios:
 do k = 1,levs+1
    do i = 1,im
       colo3(i,k) = 0.
       coloz(i,k) = 0.
    enddo
 enddo

 do l = levs,1,-1
    pmin =  1.0e10
    pmax = -1.0e10
    do i = 1,im
       wk1(i) = log(prsl(i,l))
       pmin   = min(wk1(i),pmin)
       pmax   = max(wk1(i),pmax)
       prod(i,:) = 0.0
    enddo

    kmax = 1
    kmin = 1
    do k = 1,ko3-1
       if(pmin < po3(k)) kmax = k
       if(pmax < po3(k)) kmin = k
    enddo

    do k = kmin,kmax
       temp = 1.0 / (po3(k) - po3(k+1))
       do i = 1,im
          flg(i) = .false.
          if(wk1(i) < po3(k) .and. wk1(i) >= po3(k+1)) then
             flg(i) = .true.
             wk2(i) = (wk1(i) - po3(k+1)) * temp
             wk3(i) = 1.0 - wk2(i)
          endif
       enddo

       do j = 1,pl_coeff
          do i = 1,im
             if(flg(i)) then
                prod(i,j) = wk2(i)*prdout(i,k,j) &
                          + wk3(i)*prdout(i,k+1,j)
             endif
          enddo
       enddo
    enddo

    do j = 1,pl_coeff
       do i = 1,im
          if(wk1(i) < po3(ko3)) then
             prod(i,j) = prdout(i,ko3,j)
          endif
          if(wk1(i) >= po3(1)) then
             prod(i,j) = prdout(i,1,j)
          endif
       enddo
    enddo
    do i = 1,im
       colo3(i,l) = colo3(i,l+1) + ozi(i,l)*delp(i,l)*gravi
       coloz(i,l) = coloz(i,l+1) + prod(i,6)*delp(i,l)*gravi
       prod(i,2)  = min(prod(i,2),0.0)
    enddo

    do i = 1,im
       ozib(i) = ozi(i,l)
       tem = prod(i,1) - prod(i,2) * prod(i,6)   &
           + prod(i,3) * (tin(i,l) - prod(i,5))  &
           + prod(i,4) * (colo3(i,l)-coloz(i,l))
       oz(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
    enddo
    if(ldiag3d .and. qdiag3d) then ! ozone change diagnostics
       do i = 1,im
          ozp1(i,l) = ozp1(i,l) + (prod(i,1)-prod(i,2)*prod(i,6))*dt
          ozp2(i,l) = ozp2(i,l) + (oz(i,l) - ozib(i))
          ozp3(i,l) = ozp3(i,l) + prod(i,3)*(tin(i,l)-prod(i,5))*dt
          ozp4(i,l) = ozp4(i,l) + prod(i,4)*(colo3(i,l)-coloz(i,l))*dt
       enddo
    endif

 enddo

! call mpas_log_write('--- end subroutine o3_gfs_run:')

 end subroutine o3_gfs_run

!=================================================================================================================
 end module o3_gfs
!=================================================================================================================
