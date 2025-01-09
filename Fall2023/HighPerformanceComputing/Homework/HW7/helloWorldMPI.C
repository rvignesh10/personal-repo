#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int myRank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  printf("Hello world, from rank=%d out of total world size=%d\n", myRank, np);
  MPI_Finalize();
  return 0;
}
