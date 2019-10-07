! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_cloudiness
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants, only: ep_2
 use mpas_atmphys_vars
 use module_mp_thompson_cldfra3

 use cld_fraction, only: calc_cldfraction
 
 implicit none
 private
 public:: allocate_cloudiness,   &
          deallocate_cloudiness, &
          driver_cloudiness

!MPAS driver for parameterization of the diagnostic cloud fraction.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_cloudiness:
! ----------------------------------------------
! allocate_cloudiness  : allocate local arrays for parameterization of diagnostic cloudiness.
! deallocate_cloudiness: deallocate local arrays for parameterization of diagnostic cloudiness.
! cloudiness_from_MPAS : initialize local arrays.
! cloudiness_to_MPAS   : copy local arrays to MPAS arrays.
! driver_cloudiness    : main driver (called from subroutine physics_driver).
! calc_cldincidence    : calculates the cloud fraction as 0 or 1, depending on cloud condensates.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * added subroutine cal_cldfra3 which contains the calculation of the radiative cloud fraction in conjunction
!   with the Thompson cloud microphysics.
!   Laura D. Fowler (laura@ucar.edu) / 2016-06-04.
! * initialize the local "radiation" water vapor, cloud water, cloud ice, and snow mixing ratios. When using the
!   option "cld_fraction_thompson", the "radiative" cloud water and ice paths are updated in conjunction with
!   cloud formation, but changes to the cloud water and cloud ice mixing ratios only affect the long wave and
!   short wave radiation codes.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-05.
! * since we removed the local variable radt_cld_scheme from mpas_atmphys_vars.F, now defines radt_cld_scheme
!   as a pointer to config_radt_cld_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2917-02-16.


 contains


!=================================================================================================================
 subroutine allocate_cloudiness
!=================================================================================================================

 if(.not.allocated(dx_p)     ) allocate(dx_p(ims:ime,jms:jme)             ) 
 if(.not.allocated(xland_p)  ) allocate(xland_p(ims:ime,jms:jme)          )
 if(.not.allocated(cldfrac_p)) allocate(cldfrac_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qvrad_p)  ) allocate(qvrad_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(qcrad_p)  ) allocate(qcrad_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(qirad_p)  ) allocate(qirad_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(qsrad_p)  ) allocate(qsrad_p(ims:ime,kms:kme,jms:jme)  )

 end subroutine allocate_cloudiness

!=================================================================================================================
 subroutine deallocate_cloudiness
!=================================================================================================================

 if(allocated(dx_p)     ) deallocate(dx_p     )
 if(allocated(xland_p)  ) deallocate(xland_p  )
 if(allocated(cldfrac_p)) deallocate(cldfrac_p)
 if(allocated(qvrad_p)  ) deallocate(qvrad_p  )
 if(allocated(qcrad_p)  ) deallocate(qcrad_p  )
 if(allocated(qirad_p)  ) deallocate(qirad_p  )
 if(allocated(qsrad_p)  ) deallocate(qsrad_p  )

 end subroutine deallocate_cloudiness

!=================================================================================================================
 subroutine cloudiness_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite

!input and inout arguments:
 type(mpas_pool_type),intent(in)   :: configs
 type(mpas_pool_type),intent(in)   :: mesh
 type(mpas_pool_type),intent(in)   :: sfc_input

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer:: i,j,k

 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer:: areaCell,meshDensity
 real(kind=RKIND),dimension(:),pointer:: xland

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_len_disp',len_disp)

 call mpas_pool_get_array(mesh,'areaCell'   ,areaCell   )
 call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

 call mpas_pool_get_array(sfc_input,'xland',xland)

 do j = jts,jte
    do i = its,ite
       dx_p(i,j)    = len_disp / meshDensity(i)**0.25
       xland_p(i,j) = xland(i)
    enddo

    !--- initialize the radiative cloud fraction and water vapor, cloud water, cloud ice,
    !    and snow mixing ratios:
    do k = kts,kte
    do i = its,ite
       qvrad_p(i,k,j)   = qv_p(i,k,j)
       qcrad_p(i,k,j)   = qc_p(i,k,j)
       qirad_p(i,k,j)   = qi_p(i,k,j)
       qsrad_p(i,k,j)   = qs_p(i,k,j)
       cldfrac_p(i,k,j) = 0._RKIND
    enddo
    enddo
 enddo

 end subroutine cloudiness_from_MPAS

!=================================================================================================================
 subroutine cloudiness_to_MPAS(diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer:: i,j,k

 real(kind=RKIND),dimension(:,:),pointer:: cldfrac

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'cldfrac',cldfrac)

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(k,i) = cldfrac_p(i,k,j)
 enddo
 enddo
 enddo
 
 end subroutine cloudiness_to_MPAS

!=================================================================================================================
 subroutine driver_cloudiness(configs,mesh,diag_physics,sfc_input,its,ite, &
                       errmsg,errflg)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)   :: configs
 type(mpas_pool_type),intent(in)   :: mesh
 type(mpas_pool_type),intent(in)   :: sfc_input

 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!out arguments:
 character(len=*),intent(out):: errmsg
 integer,         intent(out):: errflg

!local variables and pointers:
 character(len=StrKIND),pointer:: radt_cld_scheme

 integer:: i,j,k

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_cloudiness:')

 call mpas_pool_get_config(configs,'config_radt_cld_scheme',radt_cld_scheme)

!copy MPAS arrays to local arrays:
 call cloudiness_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)

 cld_fraction_select: select case (trim(radt_cld_scheme))
    case("cld_incidence")
      call mpas_timer_start('calc_cldincidence')
      call calc_cldincidence(cldfrac_p,qcrad_p,qirad_p,f_qc,f_qi,its,ite)
      call mpas_timer_stop('calc_cldincidence')

    case("cld_fraction")
      call mpas_timer_start('calc_cldfraction')
      call calc_cldfraction( &
           cldfrac_3d = cldfrac_p , t_p_3d  = t_p , pres_p_3d = pres_p ,            &
           qv_3d = qvrad_p , qc_3d = qcrad_p , qi_3d = qirad_p , qs_3d = qs_p ,     &
           ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,  &
           ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme ,  &
           its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte ,  &
           errmsg   = errmsg     , errflg   = errflg                                &
                      )
      call mpas_timer_stop('calc_cldfraction')

    case("cld_fraction_thompson")
      call mpas_timer_start('cal_cldfra3')
      call cal_cldfra3( &
           cldfra = cldfrac_p , qv     = qvrad_p    , qc = qcrad_p , qi  = qirad_p , &
           qs     = qsrad_p   , p      = pres_hyd_p , t  = t_p     , rho = rho_p   , &
           xland  = xland_p   , gridkm = dx_p       ,                                &
           ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,   &
           ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme ,   &
           its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                      )
      call mpas_timer_stop('cal_cldfra3')

    case default

 end select cld_fraction_select

!copy local arrays to MPAS grid:
 call cloudiness_to_MPAS(diag_physics,its,ite)

!call mpas_log_write('--- end subroutine driver_cloudiness.')

 end subroutine driver_cloudiness

!=================================================================================================================
 subroutine calc_cldincidence(cldfrac,qc,qi,f_qc,f_qi,its,ite)
!=================================================================================================================

!input arguments:
 logical,intent(in):: f_qc,f_qi
 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: qc,qi

!output arguments:
 real(kind=RKIND),intent(out),dimension(ims:ime,kms:kme,jms:jme):: cldfrac

!local variables:
 integer:: i,j,k

 real(kind=RKIND),parameter:: thresh = 1.e-06

!-----------------------------------------------------------------------------------------------------------------

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo

 if(f_qc .and. f_qi) then
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       if(qc(i,k,j)+qi(i,k,j) .gt. thresh) cldfrac(i,k,j) = 1.0_RKIND
    enddo
    enddo
    enddo
 elseif(f_qc) then
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       if(qc(i,k,j) .gt. thresh) cldfrac(i,k,j) = 1.0_RKIND
    enddo
    enddo
    enddo
 endif
    
 end subroutine calc_cldincidence

!=================================================================================================================
 end module mpas_atmphys_driver_cloudiness
!=================================================================================================================
