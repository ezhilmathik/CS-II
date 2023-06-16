#
# Makefile for finite element code TsunAWI 
#
# Please adjust Fortran compiler and flags, path to netcdf
# to your specific computer in Makefile.in

include Makefile.in

##########################################################
# No change necessary below 
#--------------------------------------------------------

.SUFFIXES: .F90 .o 

SRC =  data_types.F90 \
           parameters.F90 \
           mesh.F90 \
           elements.F90 \
           swe.F90 \
	   benchmark.F90 \
           initial_conditions.F90 \
           okada_fault_parameters.F90 \
	   output.F90 \
           output_raster.F90 \
           quakeml.F90 \
           tsunami.F90

FOXY_LIB = ./libs/foxy/libfoxy.a
FOXY_INC = -I./libs/foxy

#---------------------------------------------------

$(EXE): $(SRC)
	make clean
ifeq ($(FOXY_FLAG),-DNO_FOXY)	
	$(FC) $(FCFLAGS) $(WP_FLAG)  -o $(EXE) $(SRC) $(NC_LIB) $(NC_INC) $(FOXY_FLAG)
else
	make libfoxy.a
	$(FC) $(FCFLAGS) $(WP_FLAG)  -o $(EXE) $(SRC) $(NC_LIB) $(NC_INC) $(FOXY_FLAG) $(FOXY_LIB) $(FOXY_INC)
endif

clean :
	rm -f *~ *.o *.mod $(EXE)

allclean:
	make clean
	make -C libs/foxy clean

libfoxy.a:
	make -C libs/foxy
