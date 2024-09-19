ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = slvfe
nobase_pkgdata_SCRIPTS = \
	tools/AMBER/gen_structure tools/AMBER/gen_input \
	tools/GROMACS/gen_structure tools/GROMACS/gen_input \
	tools/NAMD/gen_structure tools/NAMD/gen_input \
	tools/LAMMPS/gen_structure tools/LAMMPS/gen_input \
	tools/CHARMM/gen_structure tools/CHARMM/gen_input \
	tools/GENESIS/gen_structure tools/GENESIS/gen_input \
	tools/GENESIS/gen_structure_AMBER tools/GENESIS/gen_structure_CHARMM \
	tools/bootstrap \
	tools/refs_converge

#	tools/GENESIS/gen_structure_GROMACS # <- In GENESIS, GROMACS is unsupported officially so not included in SCRIPTS


install-exec-hook:
	cd $(DESTDIR)$(bindir) && rm -f ermod$(EXEEXT) && $(LN_S) erdst$(EXEEXT) ermod$(EXEEXT)

uninstall-hook:
	cd $(DESTDIR)$(bindir) && rm -f ermod$(EXEEXT)

EXTRA_DIST = parameters_er.example parameters_fe.example
EXTRA_DIST += extensions/*
EXTRA_DIST += trajectory-formats/*
EXTRA_DIST += tools/AMBER/*
EXTRA_DIST += tools/GROMACS/*
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/LAMMPS/*
EXTRA_DIST += tools/CHARMM/*
EXTRA_DIST += tools/GENESIS/*
EXTRA_DIST += tools/mol_dissection
EXTRA_DIST += tools/refs_converge
EXTRA_DIST += tools/bootstrap
CLEANFILES = *.mod

extra_flags = 

# module dependencies
# slvfe part
slvfe_SOURCES = slvfe.F90 sfemain.F90 sfecorrect.F90 engmain.F90
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
slvfe_LINK = $(FC_CANONICAL) $(FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
slvfe.$(OBJEXT): slvfe.F90 sfemain.$(OBJEXT) sfecorrect.$(OBJEXT) engmain.$(OBJEXT)
sfemain.$(OBJEXT): sfemain.F90 engmain.$(OBJEXT)
sfecorrect.$(OBJEXT): sfecorrect.F90 sfemain.$(OBJEXT)
engmain.$(OBJEXT): engmain.F90

# these are necessary to compile slvfe's sources without mpi
$(top_builddir)/slvfe.$(OBJEXT): $(top_srcdir)/slvfe.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfemain.$(OBJEXT): $(top_srcdir)/sfemain.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfecorrect.$(OBJEXT): $(top_srcdir)/sfecorrect.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/engmain.$(OBJEXT): $(top_srcdir)/engmain.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
