program test_mpas_rrtmgp_thermodynamics
 use iso_fortran_env, only: real64
 use mpas_rrtmgp_thermodynamics
 implicit none
 real(real64) :: q(1,2),o(1,2),v(1,2),p(1,3),up(1,3),dn(1,3),h(1,2)
 logical :: ok
 q=0.01_real64
 o=1.e-6_real64
 call mpas_rrtmgp_water_vmr(q,v,ok)
 if(.not.ok.or.abs(v(1,1)-0.016243_real64)>1.e-5_real64) error stop "H2O"
 call mpas_rrtmgp_ozone_vmr(o,q,v,ok)
 if(.not.ok.or.abs(v(1,1)-6.095e-7_real64)>1.e-9_real64) error stop "O3"
 p(1,:)=[100000._real64,90000._real64,80000._real64]
 up=0._real64
 dn(1,:)=[100._real64,90._real64,80._real64]
 call mpas_rrtmgp_heating_from_flux(p,up,dn,h,ok)
 if(.not.ok.or.any(h>=0._real64)) error stop "flux sign"
 p(1,2)=p(1,1)
 call mpas_rrtmgp_heating_from_flux(p,up,dn,h,ok)
 if(ok) error stop "pressure check"
 print *, "PASS: gas conversions and heating checks"
end program
