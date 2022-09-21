EXECS=mesh_test fespace_test
MPICXX?=mpicxx

all: ${EXECS}

mesh_test: mesh_test.cpp
	${MPICXX} -g -o mesh_test mesh_test.cpp
fespace_test: fespace_test.cpp
	${MPICXX} -g -o fespace_test fespace_test.cpp

clean:
	rm -f ${EXECS}
