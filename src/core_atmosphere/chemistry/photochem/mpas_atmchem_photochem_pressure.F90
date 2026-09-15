module mpas_atmchem_photochem_pressure
   use mpas_kind_types, only: RKIND
   implicit none
   private
   public :: photochem_hydrostatic_pressure
contains
   ! Port of NOAA-EMC/ufsatm 0c672e0 mpas/atmos_coupling.F90:
   ! hydrostatic_pressure.  zz is d(zeta)/dz; rho_zz is dry density/zz.
   subroutine photochem_hydrostatic_pressure(zz,zgrid,rho_zz,theta_m,exner,q,index_qv, &
         gravity,cp,rd,rv_over_rd,p0,pint,pintdry,pmiddry,temp,ierr)
      real(RKIND), intent(in) :: zz(:,:),zgrid(:,:),rho_zz(:,:),theta_m(:,:),exner(:,:)
      real(RKIND), intent(in) :: q(:,:,:),gravity,cp,rd,rv_over_rd,p0
      integer, intent(in) :: index_qv
      real(RKIND), intent(out) :: pint(:,:),pintdry(:,:),pmiddry(:,:),temp(:,:)
      integer, intent(out) :: ierr
      integer :: i,k,idx,nlev,ncell,qsize
      real(RKIND) :: dz,rhodryk,rhok,sum_water,dp,dpdry,pi_top,pmid,kappa
      nlev=size(zz,1); ncell=size(zz,2); qsize=size(q,1); kappa=rd/cp; ierr=0
      if (size(zgrid,1)/=nlev+1 .or. size(q,2)/=nlev .or. index_qv<1 .or. index_qv>qsize) then
         ierr=1; return
      end if
      do i=1,ncell
         dz=zgrid(nlev+1,i)-zgrid(nlev,i)
         pi_top=exner(nlev,i)-0.5_RKIND*gravity*dz/(cp*theta_m(nlev,i))
         if (pi_top<=0._RKIND) then; ierr=2; return; end if
         pint(nlev+1,i)=p0*pi_top**(cp/rd)
         sum_water=1._RKIND
         do idx=2,qsize; sum_water=sum_water+q(idx,nlev,i); end do
         pintdry(nlev+1,i)=pint(nlev+1,i)/sum_water
         do k=nlev,1,-1
            dz=zgrid(k+1,i)-zgrid(k,i)
            rhodryk=zz(k,i)*rho_zz(k,i)
            sum_water=1._RKIND
            do idx=2,qsize; sum_water=sum_water+q(idx,k,i); end do
            rhok=sum_water*rhodryk
            dp=gravity*dz*rhok; dpdry=gravity*dz*rhodryk
            if (dpdry<=0._RKIND) then; ierr=3; return; end if
            pint(k,i)=pint(k+1,i)+dp
            pintdry(k,i)=pintdry(k+1,i)+dpdry
            pmid=((pint(k,i)**(kappa+1._RKIND)-pint(k+1,i)**(kappa+1._RKIND))/ &
                  ((kappa+1._RKIND)*dp))**(1._RKIND/kappa)
            pmid=max(pint(k+1,i)+0.05_RKIND*dp,min(pint(k,i)-0.05_RKIND*dp,pmid))
            pmiddry(k,i)=pintdry(k+1,i)+(pmid-pint(k+1,i))*dpdry/dp
            temp(k,i)=theta_m(k,i)*exner(k,i)/(1._RKIND+rv_over_rd*q(index_qv,k,i))
         end do
      end do
   end subroutine photochem_hydrostatic_pressure
end module mpas_atmchem_photochem_pressure
