#FC = gfortran 
#FFLAGS =   -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
#LDFLAGS =  -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
# FC = ifort 
FC = ifx
FFLAGS =   -O2  -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian
LDFLAGS =  -O2   -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian
FMPI = mpiifort
MPIFLAGS =  -O2 -lmpi  -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian

#FC     = pgf77
#FFLAGS = -O3  -mp  -byteswapio -mcmodel=medium -Mlarge_arrays -Mnoframe -Munroll -Knoieee
#LDFLAGS = -O3  -mp   -byteswapio -mcmodel=medium -Mlarge_arrays -Mnoframe -Munroll -Knoieee

OBJ = PMP2mod_tools.o  PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_power.o PMP2linker.o PMP2mod_analyze.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2MGsolver_fR.o PMP2MGsolver_DGP.o PMP2MGsolver_sym.o PMP2MGsolver_kmf.o PMP2MGsolver_csf.o PMP2extradof.o

PMP2main: $(OBJ)  PMP2main.o                 
	$(FC) $(LDFLAGS) -o $@.exe $^                
PMP2analyze: $(OBJ)  PMP2analyze.o                 
	$(FC) $(LDFLAGS) -o $@.exe $^                
PMP2spectrum: PMP2mod_tools.o  PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_power.o PMP2mod_analyze.o  PMP2spectrum.o                 
	$(FC) $(LDFLAGS) -o $@.exe $^                
PMP2init: PMP2init.o
	$(FC) $(LDFLAGS)  -o $@.exe $^
PMP2start: PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2mod_density.o PMP2mod_power.o PMP2start.o
	$(FC) $(LDFLAGS) -o $@.exe $^
PMP2correlation: PMP2mod_tools.o  PMP2mod_random.o PMP2correlation.o
	$(FC) $(LDFLAGS) -o $@.exe  $^
PMP2density: PMP2mod_tools.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2density.o
	$(FC) $(LDFLAGS) -o $@.exe $^
PMP2startMPI.o: PMP2startMPI.f90
	$(FMPI) $(MPIFLAGS) -c PMP2startMPI.f90 
PMP2startMPI: $(OBJ) PMP2startMPI.o
	$(FMPI) $(MPIFLAGS)  -o $@.exe $^
PMP2mainMPI.o: PMP2mainMPI.f90
	$(FMPI) $(MPIFLAGS) -c PMP2mainMPI.f90 
PMP2mainMPI: $(OBJ) PMP2mainMPI.o
	$(FMPI) $(MPIFLAGS)  -o $@.exe $^

# Compile MG program
PMP2MG: PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_power.o PMP2linker.o PMP2mod_analyze.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2MGsolver_fR.o PMP2MGsolver_DGP.o PMP2MGsolver_sym.o PMP2MGsolver_kmf.o PMP2MGsolver_csf.o PMP2extradof.o PMP2main.o
	$(FC) $(LDFLAGS) -o $@.exe $^

.f.o: 
	$(FC) -c $(FFLAGS) $<

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(FC) $(FFLAGS) -w -c $<

clean:
	rm -f *.o *.exe *.mod
