.SUFFIXES: .F .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

all:
	./../tools/manage_externals/checkout_externals --externals ./../Externals.cfg
	$(MAKE) lookup_tables core_physics_init core_physics_mmm core_UGWP core_physics_wrf core_physics_noahmp core_photochem core_physics

dummy:
	echo "****** compiling physics ******"

OBJS_init = \
	ccpp_kind_types.o        \
	mpas_atmphys_constants.o \
	mpas_atmphys_date_time.o \
	mpas_atmphys_functions.o \
	mpas_atmphys_utilities.o

OBJS = \
	mpas_atmphys_camrad_init.o         \
	mpas_atmphys_control.o             \
	mpas_atmphys_driver.o              \
	mpas_atmphys_driver_cloudiness.o   \
	mpas_atmphys_driver_convection.o   \
	mpas_atmphys_driver_gwdo.o         \
	mpas_atmphys_driver_lsm.o          \
	mpas_atmphys_driver_lsm_noahmp.o   \
	mpas_atmphys_driver_microphysics.o \
	mpas_atmphys_driver_oml.o          \
	mpas_atmphys_driver_pbl.o          \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_seaice.o       \
	mpas_atmphys_driver_sfclayer.o     \
	mpas_atmphys_finalize.o            \
	mpas_atmphys_init.o                \
	mpas_atmphys_init_microphysics.o   \
	mpas_atmphys_interface.o           \
	mpas_atmphys_landuse.o             \
	mpas_atmphys_lsm_noahinit.o        \
	mpas_atmphys_lsm_noahmpinit.o      \
	mpas_atmphys_lsm_noahmpfinalize.o  \
	mpas_atmphys_lsm_shared.o          \
	mpas_atmphys_manager.o             \
	mpas_atmphys_o3climatology.o       \
	mpas_atmphys_packages.o            \
	mpas_atmphys_driver_photochem.o    \
	mpas_atmphys_photochem_ccpp.o      \
	mpas_atmphys_rrtmg_lwinit.o        \
	mpas_atmphys_rrtmg_swinit.o        \
	mpas_atmphys_sfc_diagnostics.o     \
	mpas_atmphys_todynamics.o          \
	mpas_atmphys_update_surface.o      \
	mpas_atmphys_update.o              \
	mpas_atmphys_vars.o

lookup_tables:
	./checkout_data_files.sh

core_physics_mmm: core_physics_init
	(cd physics_mmm; $(MAKE) -f Makefile.mpas all)

core_UGWP: core_physics_init
	(cd physics_noaa/UGWP; $(MAKE) all)

core_physics_wrf: core_physics_init core_physics_mmm core_UGWP
	(cd physics_wrf; $(MAKE) all COREDEF="$(COREDEF)")

core_physics_noahmp:
	(cd physics_noahmp/utility; $(MAKE) all COREDEF="$(COREDEF)")
	(cd physics_noahmp/src; $(MAKE) all COREDEF="$(COREDEF)")
	(cd physics_noahmp/drivers/mpas; $(MAKE) all COREDEF="$(COREDEF)")

# The CCPP photochemistry subset, compiled inside its own checkout so that its "machine" module does
# not collide with the NoahMP module of the same name. See physics_ccpp/Makefile.
OBJS_photochem = \
	physics_ccpp/ccpp-physics/machine.o \
	physics_ccpp/ccpp-physics/physcons.o \
	physics_ccpp/ccpp-physics/funcphys.o \
	physics_ccpp/ccpp-physics/module_ozphys.o \
	physics_ccpp/ccpp-physics/module_h2ophys.o

core_photochem: core_physics_init
	(cd physics_ccpp; $(MAKE) all COREDEF="$(COREDEF)")

core_physics_init: $(OBJS_init)

core_physics: core_physics_wrf core_physics_noahmp core_photochem
	($(MAKE) phys_interface COREDEF="$(COREDEF)")
	ar -ru libphys.a $(OBJS_init) $(OBJS) $(OBJS_photochem)
	($(MAKE) -C ./physics_mmm -f Makefile.mpas physics_mmm_lib)
	($(MAKE) -C ./physics_wrf physics_wrf_lib)
	($(MAKE) -C ./physics_noahmp/drivers/mpas driver_lib)
	($(MAKE) -C ./physics_noahmp/src src_lib)
	($(MAKE) -C ./physics_noahmp/utility utility_lib)

phys_interface: $(OBJS)

# DEPENDENCIES:
mpas_atmphys_camrad_init.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_control.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

# mpas_atmphys_photochem_ccpp.F is the only MPAS source that USEs a CCPP module, and it has its
# own compilation rule below so that only it sees the CCPP module directory.
mpas_atmphys_photochem_ccpp.o: \
	mpas_atmphys_utilities.o

mpas_atmphys_driver_photochem.o: \
	mpas_atmphys_photochem_ccpp.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver.o: \
	mpas_atmphys_driver_photochem.o \
	mpas_atmphys_driver_cloudiness.o \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_driver_gwdo.o \
	mpas_atmphys_driver_lsm.o \
	mpas_atmphys_driver_lsm_noahmp.o \
	mpas_atmphys_driver_pbl.o \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_seaice.o \
	mpas_atmphys_driver_sfclayer.o \
	mpas_atmphys_driver_oml.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_interface.o \
	mpas_atmphys_sfc_diagnostics.o \
	mpas_atmphys_update.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_cloudiness.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_convection.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_gwdo.o: \
	mpas_atmphys_vars.o \
	mpas_atmphys_manager.o

mpas_atmphys_driver_lsm.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_lsm_noahinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_lsm_noahmp.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_microphysics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_init_microphysics.o \
	mpas_atmphys_interface.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_oml.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_pbl.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_radiation_lw.o: \
	mpas_atmphys_camrad_init.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_rrtmg_lwinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_radiation_sw.o: \
	mpas_atmphys_camrad_init.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_rrtmg_swinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_seaice.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_lsm_shared.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_sfclayer.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_finalize.o: \
	mpas_atmphys_lsm_noahmpfinalize.o

mpas_atmphys_init.o: \
	mpas_atmphys_driver_photochem.o \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_driver_lsm.o \
	mpas_atmphys_driver_microphysics.o \
	mpas_atmphys_driver_pbl.o \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_sfclayer.o \
	mpas_atmphys_lsm_noahmpinit.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_o3climatology.o \
	mpas_atmphys_vars.o

mpas_atmphys_interface.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_landuse.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_lsm_noahinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_lsm_noahmpinit.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_lsm_noahmpfinalize.o : \
	mpas_atmphys_vars.o

mpas_atmphys_manager.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_o3climatology.o \
	mpas_atmphys_update_surface.o \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_o3climatology.o: \
	mpas_atmphys_date_time.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_rrtmg_lwinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_rrtmg_swinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_sfc_diagnostics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_todynamics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_update_surface.o: \
	mpas_atmphys_date_time.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_update.o: \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_vars.o

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_wrf; $(MAKE) clean )
	( if [ -d physics_mmm ]; then cd physics_mmm; $(MAKE) -f Makefile.mpas clean; fi; )
	( cd physics_noahmp/drivers/mpas; $(MAKE) clean )
	( cd physics_noahmp/src; $(MAKE) clean )
	( cd physics_noahmp/utility; $(MAKE) clean )
	( if [ -d physics_noaa/UGWP ]; then cd physics_noaa/UGWP; $(MAKE) clean; fi )
	( cd physics_ccpp; $(MAKE) clean )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

# ----------------------------------------------------------------------
# mpas_atmphys_photochem_ccpp.F needs module_ozphys.mod, module_h2ophys.mod and, transitively,
# the CCPP machine.mod. NoahMP defines a module named Machine too (physics_noahmp/utility), so the
# CCPP directory must not be added to the shared rule below: it is put first on the include path of
# this one file only.
# ----------------------------------------------------------------------
PHOTOCHEM_INC = -I./physics_ccpp/ccpp-physics

mpas_atmphys_photochem_ccpp.o: mpas_atmphys_photochem_ccpp.F
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -I../../framework $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(PHOTOCHEM_INC) $(FCINCLUDES) -I.. -I../../framework $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(PHOTOCHEM_INC) $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../framework $(MPAS_ESMF_INC)
endif

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -I../../framework $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I./physics_mmm -I./physics_wrf -I./physics_noahmp -I./physics_noahmp/utility -I./physics_noahmp/drivers/mpas -I./physics_noahmp/src -I./physics_noaa/UGWP -I.. -I../../framework $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I./physics_mmm -I./physics_wrf -I./physics_noahmp -I./physics_noahmp/utility -I./physics_noahmp/drivers/mpas -I./physics_noahmp/src -I./physics_noaa/UGWP -I.. -I../../framework $(MPAS_ESMF_INC)
endif
