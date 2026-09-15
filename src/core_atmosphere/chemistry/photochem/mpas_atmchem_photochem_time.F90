module mpas_atmchem_photochem_time
   use machine, only: kind_phys
   implicit none
   private
   public :: photochem_forcing_day, photochem_time_indices
contains
   pure function photochem_forcing_day(jdoy, hour, first_time) result(rjday)
      integer, intent(in) :: jdoy
      real(kind_phys), intent(in) :: hour, first_time
      real(kind_phys) :: rjday
      rjday = real(jdoy, kind_phys) + hour / 24.0_kind_phys
      if (rjday < first_time) rjday = rjday + 365.0_kind_phys
   end function photochem_forcing_day

   ! Copied from the UFS photochemistry time-index contract.  TIME has
   ! NTIME+1 entries; NTIME is therefore explicit and must not be inferred.
   pure subroutine photochem_time_indices(ntime, forcing_time, rjday, n1, n2)
      integer, intent(in) :: ntime
      real(kind_phys), intent(in) :: forcing_time(:), rjday
      integer, intent(out) :: n1, n2
      integer :: j
      n2 = ntime + 1
      do j = 2, ntime
         if (rjday < forcing_time(j)) then
            n2 = j
            exit
         end if
      end do
      n1 = n2 - 1
      if (n2 > ntime) n2 = n2 - ntime
   end subroutine photochem_time_indices
end module mpas_atmchem_photochem_time
