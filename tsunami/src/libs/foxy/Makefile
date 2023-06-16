#
# Makefile for Fortran-FOSS XML parser FoXy 

include ../../Makefile.in

##########################################################
# No change necessary below 
#--------------------------------------------------------

.SUFFIXES: .f90 .F90 .o 

OBJS =  penf_global_parameters_variables.o \
	penf_b_size.o \
	penf_stringify.o \
	penf.o \
        befor64_pack_data_m.o \
	befor64.o \
	stringifor_string_t.o \
	stringifor.o \
        foxy_xml_tag.o \
	foxy_xml_file.o \
	foxy.o 


#---------------------------------------------------

libfoxy.a: $(OBJS)
	$(AR) libfoxy.a $(OBJS)

clean :
	rm -f *~ *.o *.mod $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $*.F90

.f90.o:
	$(FC) $(FCFLAGS) -c $*.f90
