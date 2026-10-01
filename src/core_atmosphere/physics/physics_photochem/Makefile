#=================================================================================================================
# The CCPP ozone and stratospheric water vapor photochemistry.
#
# This directory holds only the checkout made by manage_externals:
#
#   ccpp-physics/          the ufs-community/ccpp-physics working tree, thinned by the sparse
#                          specification below. Nothing in it is modified by MPAS.
#   ccpp_photochem.sparse   the list of files checked out (see ../../Externals.cfg).
#
# The MPAS code that drives these schemes lives one directory up, with the rest of the physics:
#
#   ../mpas_atmphys_photochem_ccpp.F    the only MPAS file that USEs the CCPP modules
#   ../mpas_atmphys_driver_photochem.F  init, run and tendency entry points
#
# The CCPP sources are compiled from inside their checkout, so that their .o and .mod files land next to
# the sources, as they do for every other external in physics/. That also keeps the CCPP "machine" module
# out of physics/, where NoahMP defines a module of the same name (physics_noahmp/utility/Machine.F90).
# Only ../mpas_atmphys_photochem_ccpp.o is compiled with this directory on its include path.
#
# The objects are archived into ../libphys.a by the physics Makefile, which lists them in OBJS_photochem.
#=================================================================================================================

.SUFFIXES: .F .o

all: dummy ccpp

dummy:
	echo "****** compiling photochem ******"

CCPP_DIR  = ccpp-physics
CCPP_ROOT = $(CCPP_DIR)/physics

CCPP_OBJS = \
	$(CCPP_DIR)/machine.o \
	$(CCPP_DIR)/physcons.o \
	$(CCPP_DIR)/funcphys.o \
	$(CCPP_DIR)/module_ozphys.o \
	$(CCPP_DIR)/module_h2ophys.o

ccpp: $(CCPP_OBJS)

clean:
	$(RM) *.o *.mod *.f *.f90
	$(RM) $(CCPP_DIR)/*.o $(CCPP_DIR)/*.mod $(CCPP_DIR)/*.f $(CCPP_DIR)/*.f90
	@# Certain systems with Intel compilers generate *.i files
	$(RM) *.i

# ----------------------------------------------------------------------
# Compiler flags for the CCPP sources
#
# The MPAS flags are used, minus two things:
#
#  1) the Fortran standard flag. module_ozphys.F90 uses the convert= open
#     specifier and funcphys.f90 mixes real kinds in intrinsic calls, both
#     rejected by gfortran under -std=f2008 as GNU extensions.
#
#  2) the free-form flag, for machine.F only, which is fixed-form. Every
#     compiler treats a .F file as fixed-form by default, so it is enough
#     to take back out the free-form flag that MPAS adds; no explicit
#     fixed-form flag is needed. -f free (Cray) is two tokens, hence the
#     subst. CCPP_FIXED_FLAG is an escape hatch for a compiler that has
#     to be told explicitly:  make ... CCPP_FIXED_FLAG=-fixed
#
#     The 72-column limit must stay in force for machine.F: its
#     continuation lines carry a trailing '&' beyond column 72 that is
#     meant to be ignored, so do not add a -ffixed-line-length flag here.
#
# MPAS single-precision builds (-DSINGLE_PRECISION) set -DSINGLE_PREC so
# that the CCPP kind_phys matches the MPAS RKIND.
# ----------------------------------------------------------------------

CCPP_FREE_FLAGS = -ffree-form -ffree -free -freeform -Mfree -Mfreeform -qfree -qfree=f90

FFLAGS_CCPP       = $(filter-out -std=f2008 -std=f2003 -std=f95 -std=f2018,$(FFLAGS))
FFLAGS_CCPP_FIXED = $(subst -f free,-f fixed,$(filter-out $(CCPP_FREE_FLAGS),$(FFLAGS_CCPP))) $(CCPP_FIXED_FLAG)

ifneq (,$(findstring SINGLE_PRECISION,$(CPPFLAGS)))
    CCPP_PREC_DEF = -DSINGLE_PREC
endif

CCPP_SRC       = $(subst $(CCPP_DIR)/,,$<)
CCPP_CPP       = $(CPP) $(CPPFLAGS) $(CCPP_PREC_DEF)
CCPP_COMPILE   = $(FC) $(CPPFLAGS) $(CCPP_PREC_DEF) $(FFLAGS_CCPP)       -I.
CCPP_COMPILE_F = $(FC) $(CPPFLAGS) $(CCPP_PREC_DEF) $(FFLAGS_CCPP_FIXED) -I.

# ----------------------------------------------------------------------
# Compilation rules for the CCPP sources
#
# machine.F is fixed-form and gets its own rule; the other four share one.
# ----------------------------------------------------------------------

$(CCPP_DIR)/machine.o: $(CCPP_ROOT)/hooks/machine.F
ifeq "$(GEN_F90)" "true"
	$(CCPP_CPP) $< > $(CCPP_DIR)/machine.f
	(cd $(CCPP_DIR); $(FC) $(FFLAGS_CCPP_FIXED) -I. -c machine.f)
else
	(cd $(CCPP_DIR); $(CCPP_COMPILE_F) -c $(CCPP_SRC))
endif

# The source file is listed first in each line below, because the shared
# recipe refers to it as $<, i.e. the first prerequisite.
$(CCPP_DIR)/physcons.o:       $(CCPP_ROOT)/hooks/physcons.F90            $(CCPP_DIR)/machine.o
$(CCPP_DIR)/funcphys.o:       $(CCPP_ROOT)/tools/funcphys.f90            $(CCPP_DIR)/machine.o $(CCPP_DIR)/physcons.o
$(CCPP_DIR)/module_ozphys.o:  $(CCPP_ROOT)/photochem/module_ozphys.F90   $(CCPP_DIR)/machine.o $(CCPP_DIR)/funcphys.o
$(CCPP_DIR)/module_h2ophys.o: $(CCPP_ROOT)/photochem/module_h2ophys.F90  $(CCPP_DIR)/machine.o

$(CCPP_DIR)/physcons.o $(CCPP_DIR)/funcphys.o $(CCPP_DIR)/module_ozphys.o $(CCPP_DIR)/module_h2ophys.o:
ifeq "$(GEN_F90)" "true"
	$(CCPP_CPP) $< > $(CCPP_DIR)/$(notdir $*).f90
	(cd $(CCPP_DIR); $(FC) $(FFLAGS_CCPP) -I. -c $(notdir $*).f90)
else
	(cd $(CCPP_DIR); $(CCPP_COMPILE) -c $(CCPP_SRC))
endif
