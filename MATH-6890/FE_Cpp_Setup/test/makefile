EXECS=mesh_test 
MPICXX?=mpicxx

all: ${EXECS}

mesh_test: mesh_test.cpp
	${MPICXX} -o mesh_test mesh_test.cpp

clean:
	rm -f ${EXECS}
