.SUFFIXES: .F .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

all: dummy core_chemistry_opp core_chemistry

dummy:
	echo "****** compiling chemistry ******"

OBJS = \
	mpas_chemistry_driver.o \
	mpas_chemistry_init.o \
	mpas_chemistry_interface.o \
	mpas_chemistry_manager.o \
	mpas_chemistry_opp.o \
	mpas_chemistry_packages.o \
	mpas_chemistry_todynamics.o

core_chemistry_opp:
	(cd chemistry_opp; $(MAKE) all COREDEF="$(COREDEF)")

core_chemistry: core_chemistry_opp
	($(MAKE) chem_interface COREDEF="$(COREDEF)")
	ar -ru libchem.a $(OBJS)

chem_interface: $(OBJS)

# DEPENDENCIES:
mpas_chemistry_driver.o: \
	mpas_chemistry_interface.o

mpas_chemistry_init.o: \
	mpas_chemistry_opp.o


clean:
	$(RM) *.o *.mod *.f90 libchem.a
	( cd chemistry_opp; $(MAKE) clean )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I. -I./chemistry_opp -I.. -I../dynamics -I../physics -I../physics/physics_wrf -I../../framework -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I. -I./chemistry_opp -I.. -I../dynamics -I../physics -I../physics/physics_wrf -I../../framework -I../../external/esmf_time_f90
endif
