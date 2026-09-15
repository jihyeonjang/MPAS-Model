module mpas_atmchem_photochem
   use machine, only: kind_phys
   use module_ozphys, only: ty_ozphys
   use module_h2ophys, only: ty_h2ophys
   use mpas_atmchem_photochem_time, only: photochem_forcing_day, photochem_time_indices
   implicit none
   private
   public :: photochem_context_type

   type :: photochem_context_type
      type(ty_ozphys) :: ozphys
      type(ty_h2ophys) :: h2ophys
      logical :: use_o3=.false., use_h2o=.false., initialized=.false.
      integer, allocatable :: jindx1_o3(:),jindx2_o3(:),jindx1_h(:),jindx2_h(:)
      real(kind_phys), allocatable :: ddy_o3(:),ddy_h(:)
      real(kind_phys), allocatable :: ozpl(:,:,:),h2opl(:,:,:)
      real(kind_phys), allocatable :: prsl(:,:),delp(:,:),temperature(:,:)
      real(kind_phys), allocatable :: o3_before(:,:),o3_work(:,:),h2o_before(:,:),h2o_work(:,:)
   contains
      procedure :: initialize => photochem_initialize
      procedure :: update => photochem_update
      procedure :: compute_tendencies => photochem_compute_tendencies
      procedure :: finalize => photochem_finalize
   end type
contains
   subroutine photochem_initialize(this,o3_scheme,h2o_scheme,o3_file,h2o_file,lat_rad, &
         nlev,no3coef,nh2ocoef,file_id)
      class(photochem_context_type),intent(inout)::this
      character(len=*),intent(in)::o3_scheme,h2o_scheme,o3_file,h2o_file
      real(kind_phys),intent(in)::lat_rad(:)
      integer,intent(in)::nlev,no3coef,nh2ocoef,file_id
      real(kind_phys),allocatable::lat_deg(:)
      integer::ncell
      ncell=size(lat_rad); this%use_o3=o3_scheme=='oz_phys_2006'.or.o3_scheme=='oz_phys_2015'
      this%use_h2o=h2o_scheme=='h2o_phys'
      allocate(lat_deg(ncell)); lat_deg=lat_rad*180.0_kind_phys/acos(-1.0_kind_phys)
      allocate(this%prsl(ncell,nlev),this%delp(ncell,nlev),this%temperature(ncell,nlev))
      allocate(this%o3_before(ncell,nlev),this%o3_work(ncell,nlev), &
               this%h2o_before(ncell,nlev),this%h2o_work(ncell,nlev))
      if(this%use_o3)then
         call this%ozphys%load_o3prog(o3_file,file_id)
         allocate(this%jindx1_o3(ncell),this%jindx2_o3(ncell),this%ddy_o3(ncell), &
                  this%ozpl(ncell,nlev,no3coef))
         call this%ozphys%setup_o3prog(lat_deg,this%jindx1_o3,this%jindx2_o3,this%ddy_o3)
      endif
      if(this%use_h2o)then
         call this%h2ophys%load(h2o_file,file_id)
         allocate(this%jindx1_h(ncell),this%jindx2_h(ncell),this%ddy_h(ncell), &
                  this%h2opl(ncell,nlev,nh2ocoef))
         call this%h2ophys%setup(lat_deg,this%jindx1_h,this%jindx2_h,this%ddy_h)
      endif
      deallocate(lat_deg); this%initialized=.true.
   end subroutine

   subroutine photochem_update(this,jdoy,hour)
      class(photochem_context_type),intent(inout)::this
      integer,intent(in)::jdoy
      real(kind_phys),intent(in)::hour
      real(kind_phys)::rjday
      integer::n1,n2
      if(this%use_o3)then
         rjday=photochem_forcing_day(jdoy,hour,this%ozphys%time(1))
         call photochem_time_indices(this%ozphys%ntime,this%ozphys%time,rjday,n1,n2)
         call this%ozphys%update_o3prog(this%jindx1_o3,this%jindx2_o3,this%ddy_o3,rjday,n1,n2,this%ozpl)
      endif
      if(this%use_h2o)then
         rjday=photochem_forcing_day(jdoy,hour,this%h2ophys%time(1))
         call photochem_time_indices(this%h2ophys%ntime,this%h2ophys%time,rjday,n1,n2)
         call this%h2ophys%update(this%jindx1_h,this%jindx2_h,this%ddy_h,rjday,n1,n2,this%h2opl)
      endif
   end subroutine

   subroutine photochem_compute_tendencies(this,o3_scheme,dt,inv_g,prsl,delp,temperature,o3,h2o,tend_o3,tend_h2o)
      class(photochem_context_type),intent(inout)::this
      character(len=*),intent(in)::o3_scheme
      real(kind_phys),intent(in)::dt,inv_g,prsl(:,:),delp(:,:),temperature(:,:),o3(:,:),h2o(:,:)
      real(kind_phys),intent(out)::tend_o3(:,:),tend_h2o(:,:)
      this%prsl=prsl;this%delp=delp;this%temperature=temperature
      this%o3_before=o3;this%h2o_before=h2o
      this%o3_work=this%o3_before;this%h2o_work=this%h2o_before
      ! Direct MPAS component contract: both sides are surface-to-TOA; never reverse k.
      if(this%use_o3)then
         select case(o3_scheme)
         case('oz_phys_2006')
            call this%ozphys%run_o3prog_2006(inv_g,dt,this%prsl,this%temperature,this%delp,this%ozpl,this%o3_work)
         case('oz_phys_2015')
            call this%ozphys%run_o3prog_2015(inv_g,dt,this%prsl,this%temperature,this%delp,this%ozpl,this%o3_work)
         end select
      endif
      if(this%use_h2o) call this%h2ophys%run(dt,this%prsl,this%h2opl,this%h2o_work)
      tend_o3=(this%o3_work-this%o3_before)/dt
      tend_h2o=(this%h2o_work-this%h2o_before)/dt
   end subroutine

   subroutine photochem_finalize(this)
      class(photochem_context_type),intent(inout)::this
      if(allocated(this%jindx1_o3))deallocate(this%jindx1_o3,this%jindx2_o3,this%ddy_o3,this%ozpl)
      if(allocated(this%jindx1_h))deallocate(this%jindx1_h,this%jindx2_h,this%ddy_h,this%h2opl)
      if(allocated(this%prsl))deallocate(this%prsl,this%delp,this%temperature, &
           this%o3_before,this%o3_work,this%h2o_before,this%h2o_work)
      this%initialized=.false.
   end subroutine
end module mpas_atmchem_photochem
