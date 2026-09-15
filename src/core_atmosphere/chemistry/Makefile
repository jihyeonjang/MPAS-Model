.SUFFIXES: .F .o


ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

all: dummy core_chemistry

dummy:
	echo "****** compiling chemistry ******"
ifeq ($(MUSICA),true)
	echo "****** building with MUSICA ******"
else
	echo "****** building without MUSICA ******"
endif

# ----------------------------------------------------------------------
# Chemistry objects
#
# OPP capability is always compiled. Whether OPP is actually used is
# controlled at runtime by config_photochem_o3_scheme in namelist.atmosphere.
# ----------------------------------------------------------------------

OBJS = \
	mpas_atmchem.o \
	mpas_atmchem_driver.o \
	mpas_atmchem_init.o \
	mpas_atmchem_interface.o \
	mpas_atmchem_manager.o \
	opp/mpas_atmchem_opp.o \
	mpas_atmchem_packages.o \
	mpas_atmchem_todynamics.o

# ----------------------------------------------------------------------
# OPP chemistry
# ----------------------------------------------------------------------

core_chemistry_opp:
	(cd opp; $(MAKE) all COREDEF="$(COREDEF)")


core_chemistry_photochem:
	$(MAKE) -C photochem all COREDEF="$(COREDEF)" FCINCLUDES="$(FCINCLUDES)"

# ----------------------------------------------------------------------
# MUSICA
# ----------------------------------------------------------------------

ifeq ($(MUSICA),true)
core_mpas_musica:
	echo "****** building with MUSICA ******"
	(cd musica; $(MAKE) all)
endif

# ----------------------------------------------------------------------
# Main chemistry library
# ----------------------------------------------------------------------

core_chemistry: core_chemistry_opp core_chemistry_photochem $(if $(filter true,$(MUSICA)),core_mpas_musica)
	$(MAKE) chem_interface COREDEF="$(COREDEF)" PHYSICS="$(PHYSICS)" CHEMISTRY="$(CHEMISTRY)"
	ar -ru libchem.a $(OBJS) photochem/*.o

ifeq ($(MUSICA),true)
	$(MAKE) -C ./musica mpas_musica_lib
endif

chem_interface: $(OBJS)

# ----------------------------------------------------------------------
# DEPENDENCIES
# ----------------------------------------------------------------------

mpas_atmchem_driver.o: mpas_atmchem_interface.o

mpas_atmchem_init.o: opp/mpas_atmchem_opp.o

# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------

clean:
	$(RM) *.o *.mod *.f90 libchem.a
	(cd opp; $(MAKE) clean)
	$(MAKE) -C photochem clean
	(cd musica; $(MAKE) clean)
	@# Certain systems with Intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

# ----------------------------------------------------------------------
# Compilation rule
# ----------------------------------------------------------------------

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(PHYSICS) $(CHEMISTRY) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I. -I./opp -I./musica -I.. -I../dynamics -I../physics -I../physics/physics_wrf -I../../framework $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(PHYSICS) $(CHEMISTRY) $(FFLAGS) -c $*.F \
		$(CPPINCLUDES) $(FCINCLUDES) -I. -I./opp -I./musica -I.. -I../dynamics -I../physics -I../physics/physics_wrf -I../../framework $(MPAS_ESMF_INC)
endif

opp/mpas_atmchem_opp.o: opp/mpas_atmchem_opp.F
	$(FC) $(CPPFLAGS) $(COREDEF) $(PHYSICS) $(CHEMISTRY) $(FFLAGS) -c $< -o $@ \
		$(CPPINCLUDES) $(FCINCLUDES) -I. -I./opp -I.. -I../physics -I../../framework $(MPAS_ESMF_INC)
