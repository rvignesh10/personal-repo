CXX= mpic++
RUN= mpirun
opt= -O3
CXXFLAGS= -fPIC -O3
EXECS= heat1dmpi
BLK_OPTIONS= -nx=100000 -debug=0 -option=0 -tFinal=0.5 -commOption=0

all: ${EXECS}

heat1dmpi:
	${CXX} -c heat1d_mpi.C
	${CXX} ${CXXFLAGS} -o $@.bin heat1d_mpi.o

cases_blocking:
	${RUN} -n=1 ./heat1dmpi.bin ${BLK_OPTIONS} -matlabFileName=heat1dNx100000Np1.m
	${RUN} -n=2 ./heat1dmpi.bin ${BLK_OPTIONS} -matlabFileName=heat1dNx100000Np2.m
	${RUN} -n=4 ./heat1dmpi.bin ${BLK_OPTIONS} -matlabFileName=heat1dNx100000Np4.m
	${RUN} -n=8 ./heat1dmpi.bin ${BLK_OPTIONS} -matlabFileName=heat1dNx100000Np8.m
	${RUN} -n=16 ./heat1dmpi.bin ${BLK_OPTIONS} -matlabFileName=heat1dNx100000Np16.m

clean: 
	rm *.o
	rm *.bin
