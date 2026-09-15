 module mpas_atmchem_todynamics
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_dimension,mpas_pool_get_array
 use mpas_atm_dimensions
 implicit none
 private
 public:: chemistry_get_tend
 contains
 subroutine chemistry_get_tend(configs,mesh,state,tend,tend_chemistry)
 type(mpas_pool_type),intent(in)::configs,mesh,state,tend_chemistry
 type(mpas_pool_type),intent(inout)::tend
 integer,pointer::nCells,nCellsSolve,index_qo3,index_qv
 integer::i,k
 real(kind=RKIND),dimension(:,:),pointer::mass,photochem_tend_o3,photochem_tend_h2o
 real(kind=RKIND),dimension(:,:,:),pointer::tend_scalars
 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)
 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_array(state,'rho_zz',mass,2)
 call mpas_pool_get_array(tend,'scalars_tend',tend_scalars)
 call mpas_pool_get_array(tend_chemistry,'photochem_tend_o3',photochem_tend_o3)
 call mpas_pool_get_array(tend_chemistry,'photochem_tend_h2o',photochem_tend_h2o)
 ! Component tendencies are mixing-ratio tendencies. Add their mass-coupled
 ! contributions; initialization belongs to the time-integration assembly.
 if(associated(photochem_tend_o3).and.index_qo3>0)then
    do i=1,nCellsSolve;do k=1,nVertLevels
       tend_scalars(index_qo3,k,i)=tend_scalars(index_qo3,k,i)+photochem_tend_o3(k,i)*mass(k,i)
    enddo;enddo
 endif
 if(associated(photochem_tend_h2o).and.index_qv>0)then
    do i=1,nCellsSolve;do k=1,nVertLevels
       tend_scalars(index_qv,k,i)=tend_scalars(index_qv,k,i)+photochem_tend_h2o(k,i)*mass(k,i)
    enddo;enddo
 endif
 end subroutine chemistry_get_tend
 end module mpas_atmchem_todynamics
