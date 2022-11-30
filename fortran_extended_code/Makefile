.SUFFIXES: .f90 .mod

EXE = extract_dr2_dr4
FC = gfortran
FC_FLAGS = -c
LAPACK_LIB = -llapack

APP_SRC =  \
	defs.f90 consts.f90 fermi.f90 matrices.f90 solver.f90 main.f90

APP_OBJ= \
	defs.o consts.o fermi.o matrices.o solver.o main.o

$(EXE) : $(APP_OBJ)
	$(FC) -o $(EXE) $(APP_OBJ) $(LAPACK_LIB)


.f90.o: $(APP_SRC)
	$(FC) $(FC_FLAGS) $(APP_SRC)       

clean:
	-rm -f *.o *.mod $(EXE)

