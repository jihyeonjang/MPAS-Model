! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
 module mpas_init_atm_o3mr
 use mpas_kind_types
 use mpas_log
 use mpas_dmpar
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_dimension,mpas_pool_get_array

 use init_atm_read_met
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_vinterp

 implicit none
 private
 public:: init_atm_o3mr


!mpas_init_atm_o3mr contains the subroutines needed for the interpolation of the firt-guess ozone mixing ratio
!to the MPAS mesh and vertical levels.
!Laura D. Fowler (laura@ucar.edu) / (2020-07-13)

!add-ons and modifications to sourcecode:
!----------------------------------------
! * switched the ozone mixing ratio (o3mr) from a diagnostic defined in the "diag" pool to a prognostic
!   variable (qo3) defined in the "state" pool as a scalar.
!   Laura D. Fowler (2022-02-09).


 contains


!=================================================================================================================
 subroutine init_atm_o3mr(mesh,configs,fg,diag,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg
 type(mpas_pool_type),intent(inout):: state

!local variables and pointers:
 character (len=StrKIND),pointer:: config_start_time
 character(len=StrKIND):: filename_gocart
 character(len=StrKIND):: initial_date,mess

 logical:: lexist

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_atm_o3mr:')

!--- horizontal interpolation of the first-guess ozone mixing ratio to the MPAS mesh:
 call init_hinterp_o3mr(configs,mesh,fg)

!--- vertical interpolation of the first-guess ozone mixing ratio from the first-guess fixed pressure levels
!    to the MPAS pressure levels:
 call init_vinterp_o3mr(configs,mesh,diag,fg,state)

 call mpas_log_write('--- end subroutine init_atm_o3mr.')
 call mpas_log_write(' ')

 end subroutine init_atm_o3mr

!=================================================================================================================
 subroutine init_hinterp_o3mr(configs,mesh,fg)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local variables:
 type(met_data) :: field !real*4 meteorological data.
 type(proj_info):: proj

 character(len=StrKIND),pointer:: config_met_prefix
 character(len=StrKIND),pointer:: config_start_time

 logical:: have_landmask

 integer,pointer:: nCells,nFGLevels
 
 integer:: icount
 integer:: i,j,k
 integer:: iCell,istatus,masked,ndims,nInterpPoints
 integer,dimension(5):: interp_list
 integer,dimension(:),pointer:: landmask
 integer,dimension(:),pointer:: mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer  :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer  :: latPoints,lonPoints
 real(kind=RKIND),dimension(:,:),pointer:: o3mr_fg,p_fg
 real(kind=RKIND),dimension(:,:),pointer:: destField2d

 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_hinterp_o3mr:')

 call mpas_pool_get_config(configs,'config_met_prefix',config_met_prefix)
 call mpas_pool_get_config(configs,'config_start_time',config_start_time)
 call mpas_pool_get_config(configs,'config_nfglevels' ,nFGLevels)

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )

 call mpas_pool_get_array(fg,'p' ,p_fg   )
 call mpas_pool_get_array(fg,'o3',o3mr_fg)

 call mpas_log_write('--- nCells    = $i',intArgs=(/nCells/)   )
 call mpas_log_write('--- nFGLevels = $i',intArgs=(/nFGLevels/))

 do i = 1,nCells
    do k = 1,nFGLevels
       o3mr_fg(k,i) = 0._RKIND
    enddo
 enddo

 istatus = 0
 call read_met_init(trim(config_met_prefix),.false.,config_start_time(1:13),istatus)
 if(istatus /= 0) then
    call mpas_log_write('**********************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Error opening initial meteorological data file ' &
                              //trim(config_met_prefix)//':'//config_start_time(1:13),messageType=MPAS_LOG_ERR)
    call mpas_log_write('**********************************************************', messageType=MPAS_LOG_CRIT)
 endif
 call read_next_met_field(field,istatus)

 do while (istatus == 0)
    if(trim(field%field) == 'LANDSEA') then
       allocate(maskslab(-2:field%nx+3,field%ny))
       maskslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       maskslab(0 ,1:field%ny) = field%slab(field%nx  ,1:field%ny)
       maskslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       maskslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       maskslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       maskslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       maskslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)
       call mpas_log_write('minval,maxval of LANDSEA = $r $r', realArgs=(/minval(maskslab), maxval(maskslab)/))
    endif
    deallocate(field%slab)
    call read_next_met_field(field, istatus)
 enddo
 call read_met_close()

 if(.not. allocated(maskslab)) then
    call mpas_log_write('**********************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('LANDSEA field not found in meteorological data file ' &
                              //trim(config_met_prefix)//':'//config_start_time(1:13),messageType=MPAS_LOG_ERR)
    call mpas_log_write('**********************************************************', messageType=MPAS_LOG_CRIT)
 endif

!horizontally interpolate O3 mixing ratio:
 istatus = 0
 call read_met_init(trim(config_met_prefix),.false.,config_start_time(1:13),istatus)
 if(istatus /= 0) then
    call mpas_log_write('**********************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Error opening initial meteorological data file ' &
                              //trim(config_met_prefix)//':'//config_start_time(1:13),messageType=MPAS_LOG_ERR)
    call mpas_log_write('**********************************************************', messageType=MPAS_LOG_CRIT)
 endif
 call read_next_met_field(field,istatus)

 icount = 0
 do while(istatus == 0)

    interp_list(1) = FOUR_POINT
    interp_list(2) = W_AVERAGE4
    interp_list(3) = W_AVERAGE16
    interp_list(4) = SEARCH
    interp_list(5) = 0

    maskval = -1.0
    masked  = -1
    fillval = 0.0
    msgval  = 0.0

    mask_array => landmask

    if(index(field%field,'O3MR') /= 0) then

       icount = icount+1
!      call mpas_log_write('--- icount = $i',intArgs=(/icount/))

       !
       !set up projection:
       !
       call map_init(proj)
   
       if(field%iproj == PROJ_LATLON) then
          call map_set(PROJ_LATLON,proj, &
                       latinc = real(field%deltalat,RKIND), &
                       loninc = real(field%deltalon,RKIND), &
                       knowni = 1.0_RKIND, &
                       knownj = 1.0_RKIND, &
                       lat1   = real(field%startlat,RKIND), &
                       lon1   = real(field%startlon,RKIND))
       elseif(field%iproj == PROJ_GAUSS) then
          call map_set(PROJ_GAUSS,proj, &
                       nlat = nint(field%deltalat), &
                       loninc = 360.0_RKIND / real(field%nx,RKIND), &
                       lat1 = real(field%startlat,RKIND), &
                       lon1 = real(field%startlon,RKIND))
       endif

       !
       !horizontally interpolate field at level k:
       !

       do k = 1,nFGLevels
          if(field%xlvl == p_fg(k,1)) exit
       enddo
       call mpas_log_write('Interpolating O3MR at $r $r $i',&
         realArgs=(/real(field%xlvl,RKIND),p_fg(k,1)/), intArgs=(/k/))
       nInterpPoints = nCells
       latPoints => latCell
       lonPoints => lonCell
       call mpas_pool_get_array(fg,'o3',destField2d)
       ndims = 2

       allocate(rslab(-2:field%nx+3,field%ny))
       rslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       rslab(0,1:field%ny)  = field%slab(field%nx  ,1:field%ny)
       rslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       rslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       rslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       rslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       rslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)

       do iCell = 1,nInterpPoints
          if(mask_array(iCell) /= masked) then
             lat = latPoints(iCell)*DEG_PER_RAD
             lon = lonPoints(iCell)*DEG_PER_RAD
             call latlon_to_ij(proj,lat,lon,x,y)
             if(x < 0.5) then
                lon = lon + 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             elseif(x > real(field%nx,kind=RKIND)+ 0.5) then
                lon = lon - 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             endif

             if(maskval /= -1.0) then
                destField2d(k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
                                              interp_list,1,maskval=maskval,mask_array=maskslab)
             else
                destField2d(k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
                                              interp_list,1)
             endif
          else
             destField2d(k,iCell) = fillval
          endif
       enddo
       deallocate(rslab)

    endif
    deallocate(field%slab)
    call read_next_met_field(field,istatus)

 enddo
 call read_met_close()

 call mpas_log_write('--- end subroutine init_hinterp_o3mr:')

 end subroutine init_hinterp_o3mr
 
!=================================================================================================================
 subroutine init_vinterp_o3mr(configs,mesh,diag,fg,state)
!=================================================================================================================
!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: fg

!inout arguments:
 type(mpas_pool_type),intent(inout):: state

!local variables and arrays:
 integer,pointer:: nCells,nFGLevels,nVertLevels
 integer,pointer:: index_qo3
 integer:: iCell,k,kk

 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),pointer  :: o3mr_fg,p_fg
 real(kind=RKIND),dimension(:,:),pointer  :: pressure,qo3
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_vinterp_o3mr:')

 call mpas_pool_get_config(configs,'config_nfglevels',nFGLevels)

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(fg,'p' ,p_fg   )
 call mpas_pool_get_array(fg,'o3',o3mr_fg)

 call mpas_pool_get_array(diag,'pressure_base',pressure)

 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)
 call mpas_pool_get_array(state,'scalars',scalars)
 qo3 => scalars(index_qo3,:,:)

 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nFGLevels))

!--- vertical interpolation to the MPAS levels:
 do iCell = 1, nCells
    sorted_arr(1,1:nFGLevels) = 0._RKIND
    sorted_arr(1,1:nFGLevels) = 0._RKIND
    do k = 1, nFGLevels
       kk = nFGLevels + 1 -k
       sorted_arr(1,kk) = p_fg(k,iCell)
       sorted_arr(2,kk) = o3mr_fg(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       qo3(k,iCell) = vertical_interp(target_p,nFGLevels-1, &
                                sorted_arr(:,1:nFGLevels-1),order=1,extrap=0)
       if(target_p >= p_fg(1,iCell) .and. k < nVertLevels) qo3(k,iCell) = qo3(k+1,iCell)
    enddo
 enddo

 if(allocated(sorted_arr)) deallocate(sorted_arr)

 call mpas_log_write('--- end subroutine init_vinterp_o3mr:')

 end subroutine init_vinterp_o3mr

!=================================================================================================================
 end module mpas_init_atm_o3mr
!=================================================================================================================
