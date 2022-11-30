.SUFFIXES: .f90 .mod

FC = f2py
FC_FLAGS = -c
F2PY_FLAGS = -m
LAPACK_LIB = -llapack

#APP_SRC =  \
#	defs.f90 consts.f90 matrices.f90 solver.f90 main.f90

#APP_OBJ= \
#	defs.o consts.o matrices.o solver.o main.o

#$(EXE) : $(APP_OBJ)
#	$(FC) -o $(EXE) $(APP_OBJ) $(LAPACK_LIB)


#.f90.o: $(APP_SRC)
#	$(FC) $(FC_FLAGS) $(APP_SRC)       

APP_SRC = \
	defs.f90 consts.f90 matrices.f90 solver.f90

f2py_module: $(APP_SRC)
	$(FC) $(F2PY_FLAGS) nuclear_moments $(FC_FLAGS) $(APP_SRC) $(LAPACK_LIB)

clean:
	-rm -f *.o *.mod *.so $(EXE)

