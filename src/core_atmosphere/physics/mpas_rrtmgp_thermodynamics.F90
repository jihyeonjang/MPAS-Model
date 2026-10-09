! MPAS-owned state/flux conversions for the future RRTMGP adapter.
module mpas_rrtmgp_thermodynamics
 use, intrinsic :: iso_fortran_env, only: real64
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(real64), parameter :: mdry=28.9647_real64, mh2o=18.01528_real64, mo3=47.9982_real64
 real(real64), parameter :: gravity=9.80665_real64, cp=1004.64_real64
 public :: mpas_rrtmgp_water_vmr, mpas_rrtmgp_ozone_vmr, mpas_rrtmgp_heating_from_flux
contains
 ! qv: kg water / kg moist air. Output: mol H2O / mol dry air.
 pure subroutine mpas_rrtmgp_water_vmr(qv, vmr, ok)
  real(real64), intent(in) :: qv(:,:)
  real(real64), intent(out) :: vmr(:,:)
  logical, intent(out) :: ok
  ok=all(shape(vmr)==shape(qv))
  if(.not.ok) return
  ok=all(ieee_is_finite(qv)).and.all(qv>=0._real64).and.all(qv<1._real64)
  if(.not.ok) return
  vmr=qv/(1._real64-qv)*mdry/mh2o
 end subroutine
 ! qo3: kg ozone / kg moist air; qv: moist-air specific humidity.
 ! Call only after confirming the MPAS tracer's exact mass basis.
 pure subroutine mpas_rrtmgp_ozone_vmr(qo3,qv,vmr,ok)
  real(real64), intent(in) :: qo3(:,:),qv(:,:)
  real(real64), intent(out) :: vmr(:,:)
  logical, intent(out) :: ok
  ok=all(shape(qo3)==shape(qv)).and.all(shape(qo3)==shape(vmr))
  if(.not.ok) return
  ok=all(ieee_is_finite(qo3)).and.all(ieee_is_finite(qv)).and. &
    all(qo3>=0._real64).and.all(qv>=0._real64).and.all(qv<1._real64)
  if(.not.ok) return
  vmr=qo3/(1._real64-qv)*mdry/mo3
 end subroutine
 ! All interfaces bottom-to-top. Flux (down - up) in W/m2.
 ! Output heating is K/s, positive for atmospheric energy convergence.
 pure subroutine mpas_rrtmgp_heating_from_flux(plev,fup,fdn,dtdt,ok)
  real(real64), intent(in) :: plev(:,:),fup(:,:),fdn(:,:)
  real(real64), intent(out) :: dtdt(:,:)
  logical, intent(out) :: ok
  integer :: k,nlev
  nlev=size(plev,2)
  ok=nlev>=2.and.all(shape(plev)==shape(fup)).and.all(shape(plev)==shape(fdn)) &
    .and.size(dtdt,1)==size(plev,1).and.size(dtdt,2)==nlev-1
  if(.not.ok) return
  ok=all(ieee_is_finite(plev)).and.all(ieee_is_finite(fup)).and.all(ieee_is_finite(fdn))
  if(.not.ok) return
  do k=1,nlev-1
   if(any(plev(:,k)<=plev(:,k+1))) then
    ok=.false.
    return
   endif
  enddo
  do k=1,nlev-1
   dtdt(:,k)=gravity/cp*((fdn(:,k+1)-fup(:,k+1))-(fdn(:,k)-fup(:,k))) &
     /(plev(:,k)-plev(:,k+1))
  enddo
 end subroutine
end module
