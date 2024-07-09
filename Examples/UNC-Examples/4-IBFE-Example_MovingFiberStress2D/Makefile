######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
#
LIBS = $(LIBS) -lnetcdf -lcurl 
#
#KD (MODULE 3803 <-- Feb. 2014 (?) Module Install)
#IBAMR_SRC_DIR = /nas02/apps/ibamr-3803/OFED-1.5.2/IBAMR2
#IBAMR_BUILD_DIR = /nas02/apps/ibamr-3803/ibamr-objs-opt
#
# KD (MODULE 2016-11  <-- Nov 2016 Module Install)
# IBAMR_SRC_DIR = /nas02/apps/ibamr-2016-11/OFED-1.5.2/sfw/ibamr/IBAMR
# IBAMR_BUILD_DIR = /nas02/apps/ibamr-2016-11/OFED-1.5.2/sfw/ibamr/ibamr-objs-opt
#
# DW (MODULE 2018-03  <-- March 2018 Module Install)
IBAMR_SRC_DIR = /nas/longleaf/apps-dogwood/ibamr/2018-03/sfw/ibamr/IBAMR
IBAMR_BUILD_DIR  =/nas/longleaf/apps-dogwood/ibamr/2018-03/sfw/ibamr/ibamr-objs-opt
#

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

## Needed for Xcode to capture compiler errors and warnings.
ifdef XCODE_VERSION_ACTUAL
CXXFLAGS += -fno-color-diagnostics
endif

######################################################################
## Build the application.
##
## NOTE: The following assumes that all .cpp files in the present
##       directory are used to build the executable.

SRC = $(wildcard *.C)
CPPFLAGS += -MD -MP
PDIM = 2
OBJS = $(SRC:%.C=%.o) $(IBAMR_LIB_2D) $(IBTK_LIB_2D)

main2d: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o $@
clean:
	$(RM) main2d
	$(RM) *.o *.lo *.objs *.ii *.int.c 
	$(RM) -r .libs

-include $(SRC:%.C=%.d)
