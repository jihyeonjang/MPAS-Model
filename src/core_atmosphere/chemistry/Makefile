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
# The ozone photochemistry capability is always compiled. Whether it is
# actually used is controlled at runtime by config_o3prognostic in
# namelist.atmosphere.
# ----------------------------------------------------------------------

OBJS = \
	mpas_atm_chemistry.o \
	mpas_chemistry_driver.o \
	mpas_chemistry_init.o \
	mpas_chemistry_interface.o \
	mpas_chemistry_manager.o \
	mpas_chemistry_ozphys.o \
	mpas_chemistry_packages.o \
	mpas_chemistry_todynamics.o

# ----------------------------------------------------------------------
# CCPP ozone photochemistry
#
# The sources are taken, unmodified, from the CCPP physics repository that
# manage_externals checks out into ../physics/physics_ccpp (see
# ../Externals.cfg). ccpp_ozone/Makefile only holds the build rules.
# ----------------------------------------------------------------------

core_ccpp_ozone:
	(cd ccpp_ozone; $(MAKE) all COREDEF="$(COREDEF)")

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

core_chemistry: core_ccpp_ozone $(if $(filter true,$(MUSICA)),core_mpas_musica)
	$(MAKE) chem_interface COREDEF="$(COREDEF)" PHYSICS="$(PHYSICS)" CHEMISTRY="$(CHEMISTRY)"
	ar -ru libchem.a $(OBJS)

ifeq ($(MUSICA),true)
	$(MAKE) -C ./musica mpas_musica_lib
endif

chem_interface: $(OBJS)

# ----------------------------------------------------------------------
# DEPENDENCIES
# ----------------------------------------------------------------------

mpas_chemistry_driver.o: mpas_chemistry_interface.o mpas_chemistry_ozphys.o

mpas_chemistry_init.o: mpas_chemistry_ozphys.o

mpas_chemistry_manager.o: mpas_chemistry_ozphys.o

# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------

clean:
	$(RM) *.o *.mod *.f90 libchem.a
	(cd ccpp_ozone; $(MAKE) clean)
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
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I. -I./ccpp_ozone -I./musica -I.. -I../dynamics -I../physics -I../physics/physics_wrf -I../../framework $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(PHYSICS) $(CHEMISTRY) $(FFLAGS) -c $*.F \
		$(CPPINCLUDES) $(FCINCLUDES) -I. -I./ccpp_ozone -I./musica -I.. -I../dynamics -I../physics -I../physics/physics_wrf -I../../framework $(MPAS_ESMF_INC)
endif
