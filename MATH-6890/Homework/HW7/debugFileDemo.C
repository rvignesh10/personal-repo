#define USE_PPP

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <unistd.h>

typedef double Real;

#include <string>
using std::string;
using std::max;

#include "parseCommand.h"

#include <ctime>

inline double getCPU()
{
  #ifndef USE_PPP
    return MPI_Wtime();
  #else
    return (1.0*std::clock())/CLOCKS_PER_SEC;
  #endif
}


Real getMaxValue(Real value, int processor=-1, MPI_Comm comm=MPI_COMM_WORLD)
{
  Real maxValue=value;
  #ifndef USE_PPP
  if (processor=-1)
    MPI_Allreduce(&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, comm);
  else
    MPI_Reduce(&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, processor, comm);
  #endif
  return maxValue;
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if(myRank == 0)
    printf("Usage: debugFileDemo -debug=i\n");

  int debug=1;
  string line;
  bool echo=false;
  for(int i=1; i<argc; i++)
  {
    line=argv[i];
    if(parseCommand(line, "-debug=", debug, echo) ) {}
  }

  FILE *debugFile = NULL;
  if (debug>0)
  {
    char debugFileName[80];
    sprintf(debugFileName, "debug/debugFileDemoNp%dProc%d.debug", np, myRank);
    debugFile=fopen(debugFileName, "w");
  }

  for( int ifile=0; ifile<=1; ifile++ ){
    FILE *file = ifile==0? stdout : debugFile;
    if ( (ifile==0 && myRank==0) || (ifile==1 && debug>0) ){
      fprintf(file, "------------------------------- DebugFileDemo ------------------------------- \n");
      fprintf(file, " np=%d, myRank=%d\n", np, myRank);
    }
  }
  Real cpu0=getCPU();
  usleep( myRank*100000 );

  Real cpu = getCPU() - cpu0;
  if (debug>0){
    fprintf(debugFile, " Time to sleep = %9.2e(s)\n", cpu);
    fflush(debugFile); // flush output to file
  }

  cpu = getMaxValue(cpu);
  if(myRank==0)
    printf(" Max time to sleep = %9.2e(s)\n", cpu);

  if (debugFile)
    fclose(debugFile);

  MPI_Finalize();
  return 0;
}
