FC   = gfortran
OPTS = -O3 #-fopenmp
CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
CPPFLAGS = #-DOPENMP
LINKOPTS = -L./ -lpoisDiag -L/usr/local/lib -lopenblas
INCLUDE  = 


SRC = poisson_direct_sp.F90 \
      poisson_direct_dp.F90 \
      poisDiag.F90 \
      param.F90

OBJS = $(addsuffix .o, $(basename $(SRC)))
FFLAGS  =  $(OPTS)
AR      = ar cru

.SUFFIXES:      .F90 .f90 .o

all : libpoisDiag

test: libpoisDiag
	$(FC) $(FFLAGS) test.F90 -o test $(LINKOPTS) $(INCLUDE)

libpoisDiag : $(OBJS)
	$(AR) libpoisDiag.a $(OBJS)

.F90.o:
	$(CPP) $(CPPFLAGS) $*.F90 > $*.f90
	$(FC) $(FFLAGS) $(OUTPUTINC) -c $*.f90

clean:
	rm -f *.f90 *.o *.a *.mod

poisDiag.o : poisson_direct_dp.o poisson_direct_sp.o param.o
poisson_direct_dp.o : param.o
poisson_direct_sp.o : param.o
