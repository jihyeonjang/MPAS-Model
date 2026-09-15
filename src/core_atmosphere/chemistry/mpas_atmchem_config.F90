module mpas_atmchem_config
   implicit none
   private
   public :: atmchem_configuration_valid
contains
   pure logical function atmchem_configuration_valid(o3_source, o3_scheme, h2o_scheme)
      character(len=*), intent(in) :: o3_source, o3_scheme, h2o_scheme
      logical :: source_ok, o3_ok, h2o_ok
      source_ok = o3_source == 'climatology' .or. o3_source == 'prognostic'
      o3_ok = o3_scheme == 'off' .or. o3_scheme == 'opp' .or. &
              o3_scheme == 'oz_phys_2006' .or. o3_scheme == 'oz_phys_2015'
      h2o_ok = h2o_scheme == 'off' .or. h2o_scheme == 'h2o_phys'
      atmchem_configuration_valid = source_ok .and. o3_ok .and. h2o_ok .and. &
           (o3_scheme == 'off' .or. o3_source == 'prognostic')
   end function atmchem_configuration_valid
end module mpas_atmchem_config
