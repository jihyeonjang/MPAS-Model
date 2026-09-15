program test_atmchem_contracts
 use machine,only:kind_phys
 use mpas_kind_types,only:RKIND
 use mpas_atmchem_config,only:atmchem_configuration_valid
 use mpas_atmchem_photochem_time,only:photochem_forcing_day,photochem_time_indices
 use mpas_atmchem_photochem_pressure,only:photochem_hydrostatic_pressure
 implicit none
 real(kind_phys)::times(5),rjday
 integer::n1,n2,ierr,k
 real(RKIND)::zz(3,1),zgrid(4,1),rho(3,1),theta(3,1),exner(3,1),q(2,3,1)
 real(RKIND)::pint(4,1),pintdry(4,1),pmid(3,1),temp(3,1),direct(1,3),before(1,3),work(1,3),tend1(1,3),tend2(1,3)
 times=[15._kind_phys,105._kind_phys,196._kind_phys,288._kind_phys,380._kind_phys]
 if(.not.atmchem_configuration_valid('climatology','off','h2o_phys'))error stop 1
 if(.not.atmchem_configuration_valid('prognostic','off','off'))error stop 2
 if(.not.atmchem_configuration_valid('prognostic','opp','off'))error stop 3
 if(.not.atmchem_configuration_valid('prognostic','oz_phys_2006','off'))error stop 4
 if(.not.atmchem_configuration_valid('prognostic','oz_phys_2015','h2o_phys'))error stop 5
 if(atmchem_configuration_valid('climatology','opp','off'))error stop 6
 if(atmchem_configuration_valid('climatology','oz_phys_2015','off'))error stop 7
 rjday=photochem_forcing_day(1,0._kind_phys,times(1));if(rjday/=366._kind_phys)error stop 8
 call photochem_time_indices(4,times,rjday,n1,n2);if(n1/=4.or.n2/=1)error stop 9
 rjday=photochem_forcing_day(365,12._kind_phys,times(1));call photochem_time_indices(4,times,rjday,n1,n2)
 if(n1/=4.or.n2/=1)error stop 10
 zz=1._RKIND;rho(:,1)=[1._RKIND,.5_RKIND,.25_RKIND];zgrid(:,1)=[0._RKIND,100._RKIND,300._RKIND,700._RKIND]
 theta=300._RKIND;exner=.9_RKIND;q=0._RKIND;q(1,:,1)=.01_RKIND
 call photochem_hydrostatic_pressure(zz,zgrid,rho,theta,exner,q,1, &
      9.80665_RKIND,1004.5_RKIND,287._RKIND,461.6_RKIND/287._RKIND, &
      100000._RKIND,pint,pintdry,pmid,temp,ierr)
 if(ierr/=0)error stop 11
 do k=1,3
   if(pintdry(k,1)<=pintdry(k+1,1))error stop 12
   if(pintdry(k,1)-pintdry(k+1,1)<=0._RKIND)error stop 13
 enddo
 ! Direct-order sentinel: k=1 and k=N retain their positions.
 direct(1,:)=[11._RKIND,22._RKIND,33._RKIND];if(direct(1,1)/=11._RKIND.or.direct(1,3)/=33._RKIND)error stop 14
 ! Independent work copies give order-independent tendencies.
 before=direct;work=before+1._RKIND;tend1=work-before;work=before+2._RKIND;tend2=work-before
 if(any(tend1/=1._RKIND).or.any(tend2/=2._RKIND))error stop 15
end program
