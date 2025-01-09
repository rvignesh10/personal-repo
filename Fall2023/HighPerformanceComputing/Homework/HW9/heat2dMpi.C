// =====================================================================
//
//     HEAT EQUATION IN TWO DIMENSIONS
//
// =====================================================================
#define DUSE_PPP

#include <stdio.h>
#include <mpi.h>

typedef double Real;
typedef int Integer;

#include <string>
using std::string;
using std::max;

#include <float.h>
#include <limits.h>
#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN

#include <ctime>
// -----------------------------------------------------------------------
// Return the current wall-clock time in seconds
// -----------------------------------------------------------------------
inline double getCPU()
{
   #ifndef DyUSE_PPP
        return MPI_Wtime();
    #else
        return (1.0*std::clock())/CLOCKS_PER_SEC;
    #endif
}

// parseCommand
#include "parseCommand.h"
// get local index
#include "getLocalIndexBounds.h"

// polynomial manufactured solution
static const Real c0=.2, c1=.1, c2=.3;
static const Real b0=1., b1=.5, b2=.25;
static const Real a0=1., a1=.3, a2=0.;
#define UTRUE(x,y,t) (b0 + (x)*( b1 + (x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 ) )
#define UTRUET(x,y,t) (b0 + (x)*( b1 +(x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a1 + 2.*(t)*a2 )
#define UTRUEXX(x,y,t) ( 2.*b2 )*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 ) )
#define UTRUEYY(x,y,t) (b0 + (x)*( b1 + (x)*b2 ))*( 2.*c2 )*( a0 + (t)*( a1 + (t)*a2 ) )
#define FORCE(x,y,t) ( UTRUET(x,y,t) - kappa*( UTRUEXX(x,y,t)+UTRUEYY(x,y,t) ) )

#include "A++.h"
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;


int main(int argc, char *argv[])
{
  Index::setBoundsCheck(on);
  enum BoundaryConditionsEnum
  {
    periodic=-1,
    dirichlet=1,
    neumann=2
  };
  enum OptionsEnum
  {
    scalarIndexing=0,
    arrayIndexing =1,
    cIndexing     =2,
    fortranRoutine=3
  };
  int option = scalarIndexing;
  const int numberOfDimensions=2;
  int debug = 1;
  int saveMatlab=2;
  std::string matlabFileName= "heat2dMPI.m";
  
  MPI_Init(&argc, &argv);
  Integer root=0;
  int myRank;
  int *np;
  //MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, np);
  
  Real tFinal= 0.1;
  Integer nx = 100;
  Integer ny = nx;
  string line;
  bool echo = root==myRank? true : false;
  for( int i=1; i<argc; i++){
    line = argv[i];
    if (parseCommand(line, "-nx=", nx, echo)) {}
    else if( parseCommand(line, "-ny=", ny, echo) ) {}
    else if( parseCommand(line, "-root=", root, echo) ) {}
    else if( parseCommand(line, "-debug=", debug, echo) ) {}
    else if( parseCommand(line, "-saveMatlab=", saveMatlab, echo) ) {}
    else if( parseCommand(line, "-matlabFileName=", matlabFileName, echo) ) {}
    else if( parseCommand(line, "-option=", option, echo) ) {}
  }

  if(myRank == root)
  {
    fprintf(stdout, "Usage: heat2d -nx=<i> -ny=<i> -option=[0|1|2|3] -tFinal=<f> -debug=<i> -saveMatlab=[0|1|2] - matlabFile=<s>\n""   option : 0=scalarIndexing, 1=arrayIndexing, 2=cIndexing, 3=fortranRoutine\n" );
  }

  FILE *debugFile = NULL;
  if (debug > 0){
    char debugFileName[80];
    sprintf(debugFileName, "debug/heat2dNp%dProc%d.debug", np, myRank);
    debugFile = fopen(debugFileName, "w");
  }
  
  Real xa=0.0, xb=1.0;
  Real ya=0.0, yb=1.0;
  Real kappa=0.1;
  Real cfl=0.9;
  
  Real dx[2];
  dx[0] = (xb-xa)/nx;
  dx[1] = (yb-ya)/ny;
  
  const Integer numGhost = 1;
  const Integer n1a   = 0;
  const Integer n1b   = nx;
  const Integer nd1a  = n1a-numGhost;
  const Integer nd1b  = n1b+numGhost;
  const Integer nd1   = nd1b-nd1a+1;
  
  Integer ny_l, n2a_l, n2b_l;
  getLocalIndexBounds(myRank, np, ny, ny_l, n2a_l, n2b_l);
  const Integer nd2a_l = n2a_l-numGhost;
  const Integer nd2b_l = n2b_l+numGhost;
  const Integer nd2_l  = nd2b_l-nd2a_l+1;
  
  IntegerArray gridIndexRange( 2, numberOfDimensions );
  IntegerArray dimension( 2, numberOfDimensions );
  IntegerArray boundaryCondition( 2, numberOfDimensions );

  gridIndexRange(0, 0) = n1a; gridIndexRange(1, 0) = n1b;
  gridIndexRange(0, 1) = n2a_l; gridIndexRange(1, 1) = n2b_l;
  
  dimension(0, 0) = nd1a; dimension(1, 0) = nd1b;
  dimension(0, 1) = nd2a_l; dimension(1, 1) = nd2b_l;

  boundaryCondition(0, 0) = dirichlet;
  boundaryCondition(1, 0) = dirichlet;

  // Grid points
  Range Rx(nd1a, nd1b), Ry(nd2a_l, nd2b_l);
  RealArray x(Rx, Ry, 2);

  int i1, i2;
  for( i2=nd2a_l; i2<=nd2b_l; i2++ )
  for( i1=nd1a; i1<=nd1b; i1++ )
  {
    x(i1, i2, 0) = xa + (i1-n1a)*dx[0];
    x(i1, i2, 1) = ya + (i2-n2a_l)*dx[1];
  }    
  
  if (debug > 0){
    fprintf(debugFile, "x(:, :, 0)=\n");
    for( i1=nd1a; i1<=nd1b; i1++ ){
    for( i2=nd2a_l; i2<=nd2b_l; i2++ )
    {
      fprintf(debugFile, "%1.5e\t", x(i1, i2, 0));  
    }
    fprintf(debugFile, "\n");
    }
    fprintf(debugFile, "x(:, :, 1)=\n");
    for( i1=nd1a; i1<=nd1b; i1++ ){
    for( i2=nd2a_l; i2<=nd2b_l; i2++ )
    {
      fprintf(debugFile, "%1.5e\t", x(i1, i2, 1));
    }
      fprintf(debugFile, "\n");
    }
    fflush(debugFile);
  }

  MPI_Finalize();
  return 0;
}


