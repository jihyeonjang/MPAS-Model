# Photochemistry component build. Invoked by chemistry/Makefile with MPAS flags.
CCPP := ccpp-physics
UPSTREAM_OBJS := machine.o physcons.o funcphys.o module_ozphys.o module_h2ophys.o
MPAS_OBJS := mpas_atmchem_photochem_time.o mpas_atmchem_photochem_pressure.o mpas_atmchem_photochem.o
OBJS := $(UPSTREAM_OBJS) $(MPAS_OBJS)
all: $(OBJS)
machine.o: $(CCPP)/physics/hooks/machine.F
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@
physcons.o: $(CCPP)/physics/hooks/physcons.F90 machine.o
	$(FC) $(FFLAGS) -I. -c $< -o $@
funcphys.o: $(CCPP)/physics/tools/funcphys.f90 machine.o physcons.o
	$(FC) $(FFLAGS) -I. -c $< -o $@
module_ozphys.o: $(CCPP)/physics/photochem/module_ozphys.F90 machine.o funcphys.o
	$(FC) $(FFLAGS) -I. -c $< -o $@
module_h2ophys.o: $(CCPP)/physics/photochem/module_h2ophys.F90 machine.o
	$(FC) $(FFLAGS) -I. -c $< -o $@
mpas_atmchem_photochem_time.o: mpas_atmchem_photochem_time.F90 machine.o
	$(FC) $(FFLAGS) -I. -c $< -o $@
mpas_atmchem_photochem_pressure.o: mpas_atmchem_photochem_pressure.F90
	$(FC) $(FFLAGS) $(FCINCLUDES) -I../.. -c $< -o $@
mpas_atmchem_photochem.o: mpas_atmchem_photochem.F90 module_ozphys.o module_h2ophys.o mpas_atmchem_photochem_time.o
	$(FC) $(FFLAGS) $(FCINCLUDES) -I. -I../.. -c $< -o $@
clean:
	$(RM) *.o *.mod
