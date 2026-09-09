#=================================================================================================================
# Build rules for the CCPP ozone photochemistry parameterization.
#
# The sources compiled here are NOT part of MPAS: they are taken, unmodified, from the CCPP physics
# repository that is checked out by manage_externals into
#
#     src/core_atmosphere/chemistry/physics_ccpp
#
# (see src/core_atmosphere/Externals.cfg).  Only the files needed by module_ozphys are compiled:
#
#     physics/hooks/machine.F         kind definitions used throughout CCPP
#     physics/hooks/physcons.F90      physical constants        (needed by funcphys)
#     physics/tools/funcphys.f90      thermodynamic functions   (module_ozphys uses fpkapx)
#     physics/photochem/module_ozphys.F90   the ozone photochemistry scheme itself
#
# Everything that adapts this code to MPAS lives in ../mpas_chemistry_ozphys.F and the other
# mpas_chemistry_*.F files, so that this directory never needs to be edited when the external
# checkout is updated.
#=================================================================================================================

CCPP_ROOT = ./../physics_ccpp/physics

all: dummy ccpp_ozone

dummy:
	echo "****** compiling ccpp_ozone ******"

OBJS = \
	machine.o        \
	physcons.o       \
	funcphys.o       \
	module_ozphys.o

# DEPENDENCIES:
physcons.o:      machine.o
funcphys.o:      machine.o physcons.o
module_ozphys.o: machine.o funcphys.o

ccpp_ozone: $(OBJS)
	ar -ru ./../libchem.a $(OBJS)

#-----------------------------------------------------------------------------------------------------------------
# Compiler flags.
#
# The CCPP sources are compiled with the MPAS flags, minus two things:
#
#  1) the Fortran standard flag.  funcphys.f90 mixes real kinds in intrinsic calls, which gfortran
#     rejects under -std=f2008 ("GNU Extension: Different type kinds").
#
#  2) the free-form flag, for machine.F only, which is fixed-form Fortran.  The free-form flag used by
#     the selected compiler is removed and the matching fixed-form flag is added instead.  Note that the
#     standard 72-column limit must be kept for this file: its continuation lines carry a trailing '&'
#     beyond column 72 that is meant to be ignored, so widening the line length breaks its declarations.
#     If your compiler is not covered below, pass the correct flag on the command line, e.g.
#
#         make ... CCPP_FIXED_FLAG=-my-fixed-form-flag
#
# MPAS single-precision builds (-DSINGLE_PRECISION) set -DSINGLE_PREC so that the CCPP kind_phys
# matches the MPAS RKIND.  The MPAS/CCPP interface converts explicitly in either case, so this only
# avoids needless copies.
#-----------------------------------------------------------------------------------------------------------------

FFLAGS_CCPP = $(filter-out -std=f2008 -std=f2003 -std=f95 -std=f2018,$(FFLAGS))

ifneq (,$(findstring -ffree-form,$(FFLAGS)))
    CCPP_FREE_FLAG = -ffree-form
    CCPP_FIXED_FLAG = -ffixed-form
else ifneq (,$(findstring -Mfreeform,$(FFLAGS)))
    CCPP_FREE_FLAG = -Mfreeform
    CCPP_FIXED_FLAG = -Mnofreeform
else ifneq (,$(findstring -Mfree,$(FFLAGS)))
    CCPP_FREE_FLAG = -Mfree
    CCPP_FIXED_FLAG = -Mnofree
else ifneq (,$(findstring -qfree,$(FFLAGS)))
    CCPP_FREE_FLAG = -qfree=f90
    CCPP_FIXED_FLAG = -qfixed
else ifneq (,$(findstring -freeform,$(FFLAGS)))
    CCPP_FREE_FLAG = -freeform
    CCPP_FIXED_FLAG = -fixedform
else ifneq (,$(findstring -ffree,$(FFLAGS)))
    CCPP_FREE_FLAG = -ffree
    CCPP_FIXED_FLAG = -ffixed
else ifneq (,$(findstring -f free,$(FFLAGS)))
    CCPP_FREE_FLAG = -f free
    CCPP_FIXED_FLAG = -f fixed
else ifneq (,$(findstring -free,$(FFLAGS)))
    CCPP_FREE_FLAG = -free
    CCPP_FIXED_FLAG = -fixed
endif

FFLAGS_CCPP_FIXED = $(filter-out $(CCPP_FREE_FLAG),$(FFLAGS_CCPP)) $(CCPP_FIXED_FLAG)

ifneq (,$(findstring SINGLE_PRECISION,$(CPPFLAGS)))
    CCPP_PREC_DEF = -DSINGLE_PREC
endif

#-----------------------------------------------------------------------------------------------------------------
# Compilation rules.  One rule per file, since the CCPP sources use three different suffixes and
# machine.F is fixed-form.
#-----------------------------------------------------------------------------------------------------------------

machine.o: $(CCPP_ROOT)/hooks/machine.F
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CCPP_PREC_DEF) $(CPPINCLUDES) $< > machine.f
	$(FC) $(FFLAGS_CCPP_FIXED) -c machine.f $(FCINCLUDES) -I.
else
	$(FC) $(CPPFLAGS) $(CCPP_PREC_DEF) $(FFLAGS_CCPP_FIXED) -c $< $(CPPINCLUDES) $(FCINCLUDES) -I.
endif

physcons.o: $(CCPP_ROOT)/hooks/physcons.F90
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CCPP_PREC_DEF) $(CPPINCLUDES) $< > physcons.f90
	$(FC) $(FFLAGS_CCPP) -c physcons.f90 $(FCINCLUDES) -I.
else
	$(FC) $(CPPFLAGS) $(CCPP_PREC_DEF) $(FFLAGS_CCPP) -c $< $(CPPINCLUDES) $(FCINCLUDES) -I.
endif

funcphys.o: $(CCPP_ROOT)/tools/funcphys.f90
	$(FC) $(FFLAGS_CCPP) -c $< $(FCINCLUDES) -I.

module_ozphys.o: $(CCPP_ROOT)/photochem/module_ozphys.F90
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CCPP_PREC_DEF) $(CPPINCLUDES) $< > module_ozphys.f90
	$(FC) $(FFLAGS_CCPP) -c module_ozphys.f90 $(FCINCLUDES) -I.
else
	$(FC) $(CPPFLAGS) $(CCPP_PREC_DEF) $(FFLAGS_CCPP) -c $< $(CPPINCLUDES) $(FCINCLUDES) -I.
endif

clean:
	$(RM) *.f *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
