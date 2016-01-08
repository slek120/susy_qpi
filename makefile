FC = gfortran
FFLAGS = --fixed-line-length-none --free-line-length-none -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
CFLAGS = 

susy_qpi : susy_qpi.f90 \
	dcuhre.o invert4x4.o \
	inc/susy_qpi.f90 inc/write_data.f90 inc/QPI.f90 inc/Gmatrix.f90
	$(FC) $(FFLAGS) -o susy_qpi.out susy_qpi.f90 dcuhre.o invert4x4.o

dcuhre.o : src/dcuhre.f
	$(FC) $(FFLAGS) -O -c src/dcuhre.f

invert4x4.o : src/invert4x4.f
	$(FC) $(FFLAGS) -O -c src/invert4x4.f