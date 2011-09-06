# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example


# trace the chain of included makefiles
makefiles += example_AMRGodunov_execPolytropic

# the base name(s) of the application(s) in this directory
ebase = amrGodunov

# the location of Chombo lib dir
HOSTNAME:=$(shell hostname)
ifeq ($(findstring cees,$(HOSTNAME)),cees)
	CHOMBO_HOME = /home/jkozdon/packages/Chombo/lib
else
ifeq ($(findstring ranger,$(HOSTNAME)),ranger)
	CHOMBO_HOME = /share/home/01643/jeremyko/packages/Chombo/lib
else
	CHOMBO_HOME = /Users/jkozdon/apps/Chombo/lib
endif
endif

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../srcLinElast ../ibcLinElast

# input file for 'run' target
INPUT = test.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# Note: tell Make.rules to not use -q for 'run' targe20044t
_space = # # NOTE: keep the space after the "=" or bad things happen
VERBOSE = $(_space)

# application-specific targets
