# Component-local source closure, consumed by core_atmosphere/CMakeLists.txt.
set(ATMOSPHERE_CORE_PHOTOCHEM_SOURCES
    chemistry/photochem/mpas_atmchem_photochem_time.F90
    chemistry/photochem/mpas_atmchem_photochem_pressure.F90
    chemistry/photochem/mpas_atmchem_photochem.F90
    chemistry/photochem/ccpp-physics/physics/hooks/machine.F
    chemistry/photochem/ccpp-physics/physics/hooks/physcons.F90
    chemistry/photochem/ccpp-physics/physics/tools/funcphys.f90
    chemistry/photochem/ccpp-physics/physics/photochem/module_ozphys.F90
    chemistry/photochem/ccpp-physics/physics/photochem/module_h2ophys.F90)
