!=================================================================================================================
!cld_fraction.F modification started from 9/26/2019 Jihyeon Jang
! add interface (or driver?) for the mpas run 10/04/2019
!=================================================================================================================
!
module cld_fraction

  use mpas_log
  !USE ccpp_kinds, ONLY: kind_phys

  IMPLICIT NONE

  PUBLIC :: cld_fraction_init
  PUBLIC :: cld_fraction_run
  PUBLIC :: cld_fraction_finalize

  PUBLIC :: cldfraction

contains
!
!-------------------------------------------------------------------------------
   subroutine cldfraction(cldfrac_3d,t_p_3d,pres_p_3d,                        &
                           qv_3d,qc_3d,qi_3d,qs_3d,                            &
                           ids,ide, jds,jde, kds,kde,                          &
                           ims,ime, jms,jme, kms,kme,                          &
                           its,ite, jts,jte, kts,kte,                          &
                           errmsg,errflg                                       &
                               )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!  input arguments:
!
   integer,  intent(in)      ::      ids,ide, jds,jde, kds,kde,                &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte

   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in)      ::                                         qv_3d, &
                                                                        qc_3d, &
                                                                        qi_3d, &
                                                                        qs_3d
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in)      ::                                        t_p_3d, &
                                                                    pres_p_3d
!
!  output arguments:
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(out)     ::                                    cldfrac_3d
!
   character(len=*), intent(out)   ::                                  errmsg
   integer,          intent(out)   ::                                  errflg
!
!  local variables:
!
   real,     dimension( its:ite, kts:kte )                 ::                  &
                                                                           qv, &
                                                                           qc, &
                                                                           qi, &
                                                                           qs, &
                                                                          t_p, &
                                                                       pres_p, &
                                                                      cldfrac
   integer :: i,j,k
!-------------------------------------------------------------------------------
!
   do j = jts,jte
!
      do k = kts, kte
         do i = its, ite
            t_p(i,k)     = t_p_3d(i,k,j)
            pres_p(i,k)  = pres_p_3d(i,k,j)
            qv(i,k)      = qv_3d(i,k,j)
            qc(i,k)      = qc_3d(i,k,j)
            qi(i,k)      = qi_3d(i,k,j)
            qs(i,k)      = qs_3d(i,k,j)
         end do
      end do
!
   call cld_fraction_run(cldfrac,t_p,pres_p,qv,qc,qi,qs,                       &
                               its,ite,kts,kte,                                &
                               errmsg,errflg)
!
      do k = kts, kte
         do i = its, ite
            cldfrac_3d(i,k,j) = cldfrac(i,k)
         end do
      end do
!
   end do

   end subroutine cldfraction
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!> \section arg_table_cld_fraction_run  Argument Table
!! \htmlinclude arg_table_cld_fraction_run.html
!!
   subroutine cld_fraction_run(cldfrac,t_p,pres_p,qv,qc,qi,qs,                 &
                               its,ite,kts,kte,                                &
                               errmsg,errflg                                   &
                               )
!-------------------------------------------------------------------------------
   use mpas_atmphys_constants, only: ep_2
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!  references: 
!-------------------------------------------------------------------------------
!
!  parameters:
!
   real,parameter    ::  alpha0  = 100.
   real,parameter    ::  gamma   = 0.49
   real,parameter    ::  qcldmin = 1.e-12
   real,parameter    ::  pexp    = 0.25
   real,parameter    ::  rhgrid  = 1.0

   real,parameter    ::  svp1  = 0.61078
   real,parameter    ::  svp2  = 17.2693882
   real,parameter    ::  svpi2 = 21.8745584
   real,parameter    ::  svp3  = 35.86
   real,parameter    ::  svpi3 = 7.66
   real,parameter    ::  svpt0 = 273.15
!
!  input arguments:
!
   integer,  intent(in)      :: its,ite,kts,kte

   real,     dimension( its:ite, kts:kte ),                                    &
             intent(in)      ::                                            qv, &
                                                                           qc, &
                                                                           qi, &
                                                                           qs
!
   real,     dimension( its:ite, kts:kte ),                                    &
             intent(in)      ::                                           t_p, &
                                                                       pres_p
!
!  output arguments:
!
   real,     dimension( its:ite, kts:kte ),                                    &
             intent(out)     ::                                       cldfrac
!
   character(len=*), intent(out)   ::                                  errmsg
   integer,          intent(out)   ::                                  errflg
!
!  local variables:
!
   integer           :: i,k
   real              :: esi,esw,qvsi,qvsw
   real              :: arg,denom,qcld,qvs,rhum,subsat,weight
!
!-------------------------------------------------------------------------------
!
   do k = kts,kte
     do i = its,ite
       cldfrac(i,k) = 0.
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite

!... calculation of the saturation mixing ratios over water and over ice (Murray, 1966):
       esw = 1000. * svp1 * exp(svp2 * (t_p(i,k) - svpt0) / (t_p(i,k) - svp3))
       esi = 1000. * svp1 * exp(svpi2 * (t_p(i,k) - svpt0) / (t_p(i,k) - svpi3))

       qvsw = ep_2 * esw / (pres_p(i,k) - esw)
       qvsi = ep_2 * esi / (pres_p(i,k) - esi)

       qcld = qc(i,k) + qi(i,k) + qs(i,k)
       if(qcld .lt. qcldmin) then
         weight = 0.
       else
         weight = (qi(i,k) + qs(i,k)) / qcld
       endif

       qvs = (1-weight) * qvsw + weight * qvsi
       rhum = qv(i,k) / qvs

       if(qcld .lt. qcldmin) then
         !assume that the cloud fraction is equal to 0. when the cloudy mixing ratio equals 0.
         cldfrac(i,k) = 0.

       elseif(rhum .ge. rhgrid) then
         !assume that the cloud fraction is equal to 1. when the relative humidity equal 100%.
         cldfrac(i,k) = 1.

       else
         !computation of the cloud fraction:
         subsat = max(1.e-10,rhgrid*qvs-qv(i,k))
         denom  = subsat**gamma
         arg    = max(-6.9,-alpha0*qcld/denom) ! exp(-6.9) = 0.001

         rhum = max(1.e-10,rhum)
         cldfrac(i,k) = (rhum/rhgrid)**pexp*(1.-exp(arg))
         if(cldfrac(i,k) .lt. 0.01) cldfrac(i,k) = 0.

       endif

     enddo
   enddo
!
   errmsg = 'cld_fraction_run OK'
   errflg = 0
!
   end subroutine cld_fraction_run

!> \section arg_table_cld_fraction_init  Argument Table
!! \htmlinclude arg_table_cld_fraction_init.html
!!
   subroutine cld_fraction_init (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

   end subroutine cld_fraction_init

!> \section arg_table_cld_fraction_finalize  Argument Table
!! \htmlinclude arg_table_cld_fraction_finalize.html
!!
   subroutine cld_fraction_finalize (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

   end subroutine cld_fraction_finalize
!=================================================================================================================
 end module cld_fraction
!=================================================================================================================
