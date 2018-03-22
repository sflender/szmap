CPP := mpicxx
CPPFLAGS := -O3 -g -fpermissive

SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) \
            -I /opt/local/include/healpix_cxx \
            -I /opt/local/include
#../../software/cfitsio 
#../../software/Healpix_3.31/src/cxx/generic_gcc/include

OBJECTS := $(OBJDIR)/main.o 

#linking to the libraries
xtract : $(OBJECTS) ; $(CPP) $(CPPFLAGS) $(OBJECTS) \
    -L/opt/local/lib\
    -lgsl -lgslcblas \
    -lhealpix_cxx \
    -lcxxsupport \
    -lcfitsio \
    -o szmap

#../../software/Healpix_3.31/src/cxx/generic_gcc/lib 
#    -L../../software/cfitsio/lib 

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp ; $(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o

clean: ; rm $(OBJDIR)/*.o szmap
