#include <stdio.h>
#include <string>
using std::string;
using std::max;

#include <mpi.h>

typedef double Real;
typedef int Integer;

#define DUSE_PPP

#include <ctime>
// ------------------------------------------------------------------------
// Return the current wall-clock time in seconds
// ------------------------------------------------------------------------
inline double getCPU()
{
  #ifndef DUSE_PPP
    return MPI_Wtime();
  #else
    return (1.0*std::clock())/CLOCKS_PER_SEC;
  #endif
}

#include "parseCommand.h"
#include "getLocalIndexBounds.h"

//#include "A++.h"
//typedef doubleSerialArray RealArray;
//typedef intSerialArray IntegerArray;

int main(int argc, char* argv[]){
  
  MPI_Init(&argc, &argv);
  Integer myRank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  
  Integer nx = 100;
  Integer ny = nx;

  enum arrayIndexing{
    scalarIndexing=0,
    arrayIndexing=1,
    cIndexing=2,
    fIndexing=3
  };

  Integer option=cIndexing; // only C-indexing available in code sorry:(

  int root = 0;
  int debug= 1;
  int saveMatlab=1;
  bool echo = root==myRank? true:false;
  string line; 
  for(int i=1; i<argc; i++)
  {
    line = argv[i];
    if( parseCommand(line, "-nx=", nx, echo) ) {}
    else if( parseCommand(line, "-ny=", ny, echo) ) {}
    else if( parseCommand(line, "-root=", root, echo) ) {}
    else if( parseCommand(line, "-option=", option, echo) ) {
      if (option!=cIndexing){abort();}
    }
    else if( parseCommand(line, "-debug=", debug, echo) ) {}
    else if( parseCommand(line, "-saveMatlab=", saveMatlab, echo) ) {}
  }

  if(debug>0)
  {
    fprintf(stdout, "nx=%d\n", nx);
    fflush(stdout);
  }

  MPI_Finalize();

  return 0;
}
