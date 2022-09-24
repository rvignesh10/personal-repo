EXECS=mesh_test fespace_test Integrator_test

MPICXX?=mpicxx

all: ${EXECS}

mesh_test: mesh_test.cpp
	${MPICXX} -g mesh_test.cpp -o mesh_test
fespace_test: fespace_test.cpp
	${MPICXX} -g fespace_test.cpp -o fespace_test
Integrator_test: Integrator_test.cpp
	${MPICXX} Integrator_test.cpp -o Integrator_test 

clean:
	rm -f ${EXECS}
