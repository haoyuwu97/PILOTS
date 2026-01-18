# Minimal Sphinx Makefile.

SPHINXBUILD ?= sphinx-build
SOURCEDIR = source

# Build into the repository's checked-in manual output directory.
BUILDDIR = ../manual

.PHONY: html clean

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR)/*
