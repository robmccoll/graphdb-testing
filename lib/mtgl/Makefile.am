ACLOCAL_AMFLAGS = -I config

SUBDIRS = mtgl test tutorial

EXTRA_DIST = \
	AUTHORS \
	ChangeLog \
	COPYING \
	FIRST_TUTORIAL \
	INSTALL \
	README \
	doc \
	win32 \
	win32_example

all-local:
	@echo ""
	@echo "The default or 'all' target does nothing as MTGL requires nothing to be built"
	@echo "before installation.  Available non-standard targets are"
	@echo ""
	@echo "    buildtests    Make the standard tests."
	@echo "    buildall      Make all the tests."
	@echo "    tutorial      Make the tutorial programs."

.PHONY: tutorial buildtests buildall

tutorial:
	$(MAKE) -C tutorial tutorial

buildtests:
	$(MAKE) -C test buildtests

buildall:
	$(MAKE) -C test buildall

if ON_MTA
mostlyclean-local:
	/bin/rm -f *.ti
	/bin/rm -f *.ii
	/bin/rm -f *.pl
endif
