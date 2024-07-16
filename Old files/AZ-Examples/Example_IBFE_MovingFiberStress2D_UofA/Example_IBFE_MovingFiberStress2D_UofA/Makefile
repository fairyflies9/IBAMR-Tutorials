######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
#
#
#KD (MODULE 3803 <-- Feb. 2014 (?) Module Install)
#IBAMR_SRC_DIR = /nas02/apps/ibamr-3803/OFED-1.5.2/IBAMR2
#IBAMR_BUILD_DIR = /nas02/apps/ibamr-3803/ibamr-objs-opt
#
# KD (MODULE 2016-11  <-- Nov 2016 Module Install)
#
IBAMR_SRC_DIR = /groups/lauram9/ib10/ibamr/IBAMR
IBAMR_BUILD_DIR = /groups/lauram9/ib10/ibamr/ibamr-objs-opt
LIBSNEW = $(LIBS) -lnetcdf -lcurl 

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

######################################################################
## Build the IB tester application.
SOURCES = main.C 
OBJS = main.o


default:
	@echo "make one of: main2d, main3d"

main2d:
	if (test -f stamp-3d); then $(MAKE) clean; fi
	touch stamp-2d
	$(MAKE) PDIM=2 main-2d

main3d:
	if (test -f stamp-2d); then $(MAKE) clean; fi
	touch stamp-3d
	$(MAKE) PDIM=3 main-3d

main-2d: $(IBAMR_LIB_2D) $(IBTK_LIB_2D) $(OBJS) $(SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) \
	$(IBAMR_LIB_2D) $(IBTK_LIB_2D) $(LIBSNEW) -DNDIM=$(PDIM) -o main2d

main-3d: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS) $(SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) \
	$(IBAMR_LIB_3D) $(IBTK_LIB_3D) -DNDIM=$(PDIM) $(LIBSNEW) -o main3d

clean:
	$(RM) main2d main3d core
	$(RM) *.o *.lo *.objs *.ii *.int.c stamp-[23]d
	$(RM) -r .libs
