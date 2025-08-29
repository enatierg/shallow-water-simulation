LIB_LIST = -lblas -llapack

FC = gfortran
OPT = -g -Wall -fimplicit-none

all: shallow_water

initial_condition.o: initial_condition.f90
	${FC} -c ${OPT} $<

create_matrices.o: create_matrices.f90 helmholtz_solve.o
	${FC} -c ${OPT} $<

helmholtz_solve.o: helmholtz_solve.f90
	${FC} -c ${OPT} $<

timestepping.o: timestepping.f90 helmholtz_solve.o
	${FC} -c ${OPT} ${LIB_LIST} $<

save_fields.o: save_fields.f90
	${FC} -c ${OPT} $<

test_helmholtz.o: test_helmholtz.f90 helmholtz_solve.o
	${FC} -c ${OPT} $<

# Main program
shallow_water.o: shallow_water.f90 initial_condition.o create_matrices.o timestepping.o save_fields.o helmholtz_solve.o test_helmholtz.o
	${FC} -c ${OPT} $<

shallow_water: shallow_water.o initial_condition.o create_matrices.o timestepping.o save_fields.o helmholtz_solve.o test_helmholtz.o
	${FC} -o shallow_water $^ ${LIB_LIST}

clean:
	rm -f *.o core.* shallow_water
