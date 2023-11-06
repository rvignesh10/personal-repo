#include <stdio.h>
#include <string>
using std::string;
using std::max;
#include <cmath>
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
#include "writeMatlabArray.h"

//#include "A++.h"
//typedef doubleSerialArray RealArray;
//typedef intSerialArray IntegerArray;
// polynomial manufactured solution
static const Real c0=.2, c1=.1, c2=.3;
static const Real b0=1., b1=.5, b2=.25;
static const Real a0=1., a1=.3, a2=0.;
#define UTRUE(x,y,t) (b0 + (x)*( b1 + (x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 ) )
#define UTRUET(x,y,t) (b0 + (x)*( b1 +(x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a1 + 2.*(t)*a2 )
#define UTRUEXX(x,y,t) ( 2.*b2 )*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 ) )
#define UTRUEYY(x,y,t) (b0 + (x)*( b1 + (x)*b2 ))*( 2.*c2 )*( a0 + (t)*( a1 + (t)*a2 ) )
#define FORCE(x,y,t) ( UTRUET(x,y,t) - kappa*( UTRUEXX(x,y,t)+UTRUEYY(x,y,t) ) )

int main(int argc, char* argv[]){
  
  MPI_Init(&argc, &argv);
  Integer myRank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Integer nx = 100;
  Integer ny = nx;

  enum BoundaryConditionsEnum
  {
    periodic  = -1,
    dirichlet = 1,
    neumann   = 2
  };

  enum arrayIndexing{
    scalarIndexing=0,
    arrayIndexing=1,
    cIndexing=2,
    fIndexing=3
  };

  enum sendTag{
    send_top=0,
    send_btm=1
  };
  enum recvTag{
    recv_btmG=0,
    recv_topG=1
  };

  Integer option=cIndexing; // only C-indexing available in code sorry:(
  string optionName="c-Indexing";
  int root = 0;
  int debug= 1;
  int saveMatlab=1;
  string matlabFileName= "heat2dMPI.m";

  Real tFinal=0.1;
  Real kappa=0.1;

  if(myRank == root)
  {
    printf("Usage: heat2d -nx=<i> -ny=<i> -option=[0|1|2|3] -tFinal=<f> -debug=<i> -saveMatlab=[0|1|2] - matlabFile=<s>\n""   option : 0=scalarIndexing, 1=arrayIndexing, 2=cIndexing, 3=fortranRoutine\n");
  }

  bool echo = root==myRank? true:false;
  string line; 
  for(int i=1; i<argc; i++)
  { 
    line = argv[i];
    if( parseCommand(line, "-nx=", nx, echo) ) {ny=nx;}
    else if( parseCommand(line, "-ny=", ny, echo) ) {}
    else if( parseCommand(line, "-root=", root, echo) ) {}
    else if( parseCommand(line, "-option=", option, echo) ) {
      if (option!=cIndexing){abort();}
    }
    else if( parseCommand(line, "-debug=", debug, echo) ) {}
    else if( parseCommand(line, "-saveMatlab=", saveMatlab, echo) ) {}
    else if( parseCommand(line, "-tFinal=", tFinal, echo) ) {}
  }

  FILE *debugFile=NULL;
  if(debug>0)
  {
    char debugFileName[80];
    sprintf(debugFileName, "debug/heat2dNp%dProc%d.debug", np, myRank);
    debugFile = fopen(debugFileName, "w");
  }
  
  Real xa=0.0, xb=1.0;
  Real ya=0.0, yb=1.0;
  Real dx[2]; 
  dx[0] = (xb-xa)/nx;
  dx[1] = (yb-ya)/ny;

  const Integer numGhost=1;
  const Integer n1a=0;
  const Integer n1b=nx;
  const Integer nd1a=n1a-numGhost;
  const Integer nd1b=n1b+numGhost;
  const Integer nd1= nd1b-nd1a+1;

  Integer ny_l, n2a_l, n2b_l;
  getLocalIndexBounds(myRank, np, ny, ny_l, n2a_l, n2b_l);
  const Integer nd2a_l=n2a_l-numGhost;
  const Integer nd2b_l=n2b_l+numGhost;
  const Integer nd2_l= nd2b_l-nd2a_l+1;

  // Local 2D grid
  Real *x_p[2];
  x_p[0] = new Real [nd1*nd2_l];
  x_p[1] = new Real [nd1*nd2_l];
  #define x(i1, i2, i3) x_p[i3][((i1)-nd1a) + nd1*((i2)-nd2a_l)]
  
  for(int j=nd2a_l; j<=nd2b_l; j++)
  for(int i=nd1a; i<=nd1b; i++)
  {
    x(i, j, 0) = xa + (i)*dx[0];
    x(i, j, 1) = ya + (j)*dx[1];
  }

  if(debug>0)
  { 
    fprintf(debugFile, "x-direction discretization:\nnd1a=%d, n1a=%d, n1b=%d, nd1b=%d, nd1=%d\n", nd1a, n1a, n1b, nd1b, nd1);
    fprintf(debugFile, "y-direction (local) discretization:\nnd2a=%d, n2a=%d, n2b=%d, nd2b=%d, nd2=%d\n", nd2a_l, n2a_l, n2b_l, nd2b_l, nd2_l);
    if (debug>1)
    {
      fprintf(debugFile, "x(:, :, 0) = \n");
      for(int j=nd2a_l; j<=nd2b_l; j++){
      for(int i=nd1a; i<=nd1b; i++)
      {
        fprintf(debugFile, "%1.5e\t", x(i, j, 0));
      }
        fprintf(debugFile, "\n");
      } 
      fprintf(debugFile, "\nx(:, :, 1) = \n");
      for(int j=nd2a_l; j<=nd2b_l; j++){
      for(int i=nd1a; i<=nd1b; i++)
      {
        fprintf(debugFile, "%1.5e\t", x(i, j, 1));
      }
        fprintf(debugFile, "\n");
      }    
    }
    fflush(debugFile);
  }

  FILE *rootFile = stdout;
  if(myRank == root)
  {
    fprintf(rootFile, "----- Solve the Heat Equation in two dimensions, np=%d ------\n", np);
    fprintf(rootFile, "      option=%d : %s \n",option, optionName.c_str());
    fprintf(rootFile, "      saveMatlab=%d, matlabFileName=%s \n",saveMatlab,matlabFileName.c_str());
    fprintf(rootFile, "      kappa=%.3g, nx=%d, ny=%d, tFinal=%6.2f\n",kappa,nx,ny,tFinal);
    fflush(rootFile);
  }

  const Integer numberOfDimensions=2;
  Integer gridIndexRange_p[2][numberOfDimensions];
  #define gridIndexRange(side, axis) gridIndexRange_p[side][axis]
  Integer dimension_p[2][numberOfDimensions];
  #define dimension(side, axis) dimension_p[side][axis]
  Integer boundaryCondition_p[2][numberOfDimensions];
  #define boundaryCondition(side, axis) boundaryCondition_p[side][axis]

  gridIndexRange(0, 0)= n1a;   gridIndexRange(1, 0)= n1b;
  gridIndexRange(0, 1)= n2a_l; gridIndexRange(1, 1)= n2b_l;

  dimension(0, 0)= nd1a;   dimension(1, 0)= nd1b;

  boundaryCondition(0, 0)= dirichlet; boundaryCondition(1, 0)= dirichlet; // left right
  boundaryCondition(0, 1)= dirichlet; boundaryCondition(1, 1)= dirichlet; // top bottom

  Real *u_p[2];
  u_p[0] = new Real [nd1 * nd2_l];
  u_p[1] = new Real [nd1 * nd2_l];
  #define uc(id1, id2) u_p[curr][((id1)-nd1a) + nd1*((id2)-nd2a_l)]
  #define un(id1, id2) u_p[next][((id1)-nd1a) + nd1*((id2)-nd2a_l)]
  // initial conditions
  Real t=0.0;
  Integer curr = 0;
  for(int j=nd2a_l; j<=nd2b_l; j++)
  for(int i=nd1a; i<=nd1b; i++)
  uc(i, j) = UTRUE( x(i, j, 0), x(i, j, 1), t );

  if(debug>1)
  { 
    fprintf(debugFile, "\nu(:, :, t=0): \n");
    for(int j=nd2a_l; j<=nd2b_l; j++){
      for(int i=nd1a; i<=nd1b; i++){
        fprintf(debugFile, "%1.5e\t", uc(i, j));
      }
      fprintf(debugFile, "\n");
    }
    fprintf(debugFile,"\n");
    fflush(debugFile);
  }
  // Time-step restriction:
  // Forward Euler: kappa*dt*( 1/dx^2 + 1/dy^2 ) < cfl*.5
  Real cfl=.9;
  Real dt = cfl*(.5/kappa)/( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) );

  int numSteps=ceil(tFinal/dt);
  dt = tFinal/numSteps; // adjust dt to reach the final time
 
  Real rx = kappa*(dt/(dx[0]*dx[0]));
  Real ry = kappa*(dt/(dx[1]*dx[1]));

  int n;
  Real cpu0= getCPU();
  for(n=0; n<numSteps; n++)
  {
    t = n*dt;
    curr = n%2;
    int next= (n+1)%2;

    if(option==cIndexing)
    {
      // do internal time-stepping
      for(int j=n2a_l; j<=n2b_l; j++)
      for(int i=n1a; i<=n1b; i++)
      {
        un(i, j) = uc(i, j) + rx*( uc(i+1, j) - 2.*uc(i, j) + uc(i-1, j) )
                            + ry*( uc(i, j+1) - 2.*uc(i, j) + uc(i, j-1) ) + dt*FORCE( x(i,j,0), x(i,j,1), t );
      }
    }
    else{
      if(myRank==root)
        printf("ERROR: indexing option not supported \n");
      abort();
    }
    // set boundary conditions
    if(np>1)
    {
    if( myRank==0 )
    {
      // set side boundary conditions - left and right
      for(int j=n2a_l; j<=n2b_l; j++)
      {
        un(n1a, j) = UTRUE( x(n1a, j, 0), x(n1a, j, 1), t+dt );
        //un(nd1a, j)= UTRUE( x(nd1a,j, 0), x(nd1a,j, 1), t+dt );
        un(nd1a,j) = 3.*un(n1a, j) - 3.*un(n1a+1, j) + un(n1a+2, j); // extrapolate
        un(n1b, j) = UTRUE( x(n1b, j, 0), x(n1b, j, 1), t+dt );
        //un(nd1b, j)= UTRUE( x(nd1b,j, 0), x(nd1b,j, 1), t+dt );
        un(nd1b,j) = 3.*un(n1b, j) - 3.*un(n1b-1, j) + un(n1b-2, j); // extrapolate
      }
      // send top
      if(debug>3)
      {
        fprintf(debugFile, "sending top from myRank=%d\n", myRank);
        for(int i=nd1a; i<=nd1b; i++)
          fprintf(debugFile, "%9.2e\t", un(i, n2b_l));
        fprintf(debugFile, "\n\n");
        fflush(debugFile);
      }
      MPI_Send(&un(nd1a,n2b_l), nd1, MPI_DOUBLE, myRank+1, send_top, MPI_COMM_WORLD);
      // fill in bottom side
      for(int i=n1a; i<=n1b; i++)
      {
        un(i, n2a_l) = UTRUE( x(i, n2a_l, 0), x(i, n2a_l, 1), t+dt );
        //un(i, nd2a_l)= UTRUE( x(i, nd2a_l,0), x(i, nd2a_l,1), t+dt );
        un(i, nd2a_l)= 3.*un(i, n2a_l) - 3.*un(i, n2a_l+1) + un(i, n2a_l+2);
      }
      // fill in corners
      un(nd1a, nd2a_l) = UTRUE( x(nd1a, nd2a_l, 0), x(nd1a, nd2a_l, 1), t+dt );
      un(nd1b, nd2a_l) = UTRUE( x(nd1b, nd2a_l, 0), x(nd1b, nd2a_l, 1), t+dt );
      // receive the topGhost layer
      MPI_Recv(&un(nd1a, nd2b_l), nd1, MPI_DOUBLE, myRank+1, recv_topG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if(myRank==np-1)
    {
      // set side boundary conditions - left and right
      for(int j=n2a_l; j<=n2b_l; j++)
      {
        un(n1a, j) = UTRUE( x(n1a, j, 0), x(n1a, j, 1), t+dt );
        //un(nd1a, j)= UTRUE( x(nd1a,j, 0), x(nd1a,j, 1), t+dt );
        un(nd1a,j) = 3.*un(n1a, j) - 3.*un(n1a+1, j) + un(n1a+2, j); // extrapolate
        un(n1b, j) = UTRUE( x(n1b, j, 0), x(n1b, j, 1), t+dt );
        //un(nd1b, j)= UTRUE( x(nd1b,j, 0), x(nd1b,j, 1), t+dt );
        un(nd1b,j) = 3.*un(n1b, j) - 3.*un(n1b-1, j) + un(n1b-2, j); // extrapolate
      }
      // send btm
      MPI_Send(&un(nd1a,n2a_l), nd1, MPI_DOUBLE, myRank-1, send_btm, MPI_COMM_WORLD);
      // fill in top side
      for(int i=n1a; i<=n1b; i++)
      {
        un(i, n2b_l) = UTRUE( x(i, n2b_l, 0), x(i, n2b_l, 1), t+dt );
        //un(i, nd2b_l)= UTRUE( x(i, nd2b_l,0), x(i, nd2b_l,1), t+dt );
        un(i, nd2b_l)= 3.*un(i, n2b_l) - 3.*un(i, n2b_l-1) + un(i, n2b_l-2);
      }
      // fill in corners
      un(nd1a, nd2b_l) = UTRUE( x(nd1a, nd2b_l, 0), x(nd1a, nd2b_l, 1), t+dt );
      un(nd1b, nd2b_l) = UTRUE( x(nd1b, nd2b_l, 0), x(nd1b, nd2b_l, 1), t+dt );
      // receive the bottom ghost layer
      MPI_Recv(&un(nd1a, nd2a_l), nd1, MPI_DOUBLE, myRank-1, recv_btmG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
      // set side boundary conditions - left and right
      for(int j=n2a_l; j<=n2b_l; j++)
      {
        un(n1a, j) = UTRUE( x(n1a, j, 0), x(n1a, j, 1), t+dt );
        un(nd1a, j)= UTRUE( x(nd1a,j, 0), x(nd1a,j, 1), t+dt );
        //un(nd1a,j) = 3.*un(n1a, j) - 3.*un(n1a+1, j) + un(n1a+2, j); // extrapolate
        un(n1b, j) = UTRUE( x(n1b, j, 0), x(n1b, j, 1), t+dt );
        un(nd1b, j)= UTRUE( x(nd1b,j, 0), x(nd1b,j, 1), t+dt );
        //un(nd1b,j) = 3.*un(n1b, j) - 3.*un(n1b-1, j) + un(n1b-2, j); // extrapolate
      }
      // send top and bottom
      MPI_Send(&un(nd1a, n2a_l), nd1, MPI_DOUBLE, myRank-1, send_btm, MPI_COMM_WORLD);
      MPI_Send(&un(nd1a, n2b_l), nd1, MPI_DOUBLE, myRank+1, send_top, MPI_COMM_WORLD);
      // recieve topG and bottomG
      MPI_Recv(&un(nd1a, nd2a_l), nd1, MPI_DOUBLE, myRank-1, recv_btmG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&un(nd1a, nd2b_l), nd1, MPI_DOUBLE, myRank+1, recv_topG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    }
    else{
      for(int j=n2a_l;j<=n2b_l;j++)
      {
        un(n1a, j) = UTRUE( x(n1a, j, 0), x(n1a, j, 1), t+dt );
        //un(nd1a, j)= UTRUE( x(nd1a,j, 0), x(nd1a,j, 1), t+dt );                                                                                                            
        un(nd1a,j) = 3.*un(n1a, j) - 3.*un(n1a+1, j) + un(n1a+2, j); // extrapolate                                                                                          
        un(n1b, j) = UTRUE( x(n1b, j, 0), x(n1b, j, 1), t+dt );
        //un(nd1b, j)= UTRUE( x(nd1b,j, 0), x(nd1b,j, 1), t+dt );                                                                                                            
        un(nd1b,j) = 3.*un(n1b, j) - 3.*un(n1b-1, j) + un(n1b-2, j); // extrapolate 
      }
      for(int i=n1a; i<=n1b; i++)
      {
        un(i, n2a_l) = UTRUE( x(i, n2a_l, 0), x(i, n2a_l, 1), t+dt );
        //un(i, nd2a_l)= UTRUE( x(i, nd2a_l,0), x(i, nd2a_l,1), t+dt );                                                                                                      
        un(i, nd2a_l)= 3.*un(i, n2a_l) - 3.*un(i, n2a_l+1) + un(i, n2a_l+2);
        un(i, n2b_l) = UTRUE( x(i, n2b_l, 0), x(i, n2b_l, 1), t+dt );
        //un(i, nd2b_l)= UTRUE( x(i, nd2b_l,0), x(i, nd2b_l,1), t+dt );                                                                                                      
        un(i, nd2b_l)= 3.*un(i, n2b_l) - 3.*un(i, n2b_l-1) + un(i, n2b_l-2);
      }  
      un(nd1a, nd2a_l) = UTRUE( x(nd1a, nd2a_l, 0), x(nd1a, nd2a_l, 1), t+dt );
      un(nd1b, nd2a_l) = UTRUE( x(nd1b, nd2a_l, 0), x(nd1b, nd2a_l, 1), t+dt );
      un(nd1a, nd2b_l) = UTRUE( x(nd1a, nd2b_l, 0), x(nd1a, nd2b_l, 1), t+dt );
      un(nd1b, nd2b_l) = UTRUE( x(nd1b, nd2b_l, 0), x(nd1b, nd2b_l, 1), t+dt );
    }
    // check for error if debug
    if(debug>1)
    {
      // check and write max errors
      Real maxErr = 0.;
      for(int j=nd2a_l; j<=nd2b_l; j++)
      for(int i=nd1a; i<=nd1b; i++)
      {
        Real err = fabs( un(i, j) - UTRUE( x(i, j, 0), x(i, j, 1), t+dt ) );
        maxErr = max(maxErr, err);
      }
      fprintf(debugFile, "step=%d, t=%9.3e, maxErr=%9.2e\n", n+1, t+dt, maxErr);
      fflush(debugFile);
    }
  }
  Real cpuTimeStep  = getCPU()-cpu0;
  Real cpuTimeStep_g= cpuTimeStep;
  MPI_Reduce(&cpuTimeStep, &cpuTimeStep_g, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

  // --------- check final error -------------
  t += dt;
  if( fabs(t-tFinal) > 1e-03*dt/tFinal )
  {
    if(myRank==root)
    {
      fprintf(stdout, "ERROR: AFTER TIME STEPPING: t=%16.8e IS NOT EQUAL TO tFinal=%16.8e\n", t, tFinal);
    }
    abort();
  }
  curr = numSteps%2;
  Real *err_p = new Real [nd1*nd2_l];
  Real maxErr_l=0., maxErr_g;
  #define err(i1, i2) err_p[((i1)-nd1a) + nd1*( (i2)-nd2a_l )]
  for(int j=nd2a_l; j<=nd2b_l; j++)
  for(int i=nd1a; i<=nd1b; i++)
  {
    err(i, j) = uc(i, j) - UTRUE( x(i, j, 0), x(i, j, 1), t );
    maxErr_l= max( maxErr_l, fabs(err(i, j)) );
  }
  maxErr_g = maxErr_l;
  MPI_Reduce(&maxErr_l, &maxErr_g, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
  if(debug>0)
  { 
    if(debug>1)
    {
      fprintf(debugFile, "\nu(:, :, t=%f): \n",tFinal);
      for(int j=nd2a_l; j<=nd2b_l; j++){
        for(int i=nd1a; i<=nd1b; i++){
          fprintf(debugFile, "%1.5e\t", uc(i, j));
        }
        fprintf(debugFile, "\n");
      }
      fprintf(debugFile,"\n");      
    }
    fprintf(debugFile, "maxErr local : %1.9e\n", maxErr_l);
    fprintf(debugFile, "maxErr global: %1.9e\n", maxErr_g);
    fflush(debugFile);
  }
  if(myRank==root)
  {
  fprintf(stdout,"\n maxErr global: %1.9e, max CPUTime: %1.9e\n", maxErr_g, cpuTimeStep_g);
  fflush(stdout);
  }

  Integer NTot=0;
  Integer YsendCount;
  Integer XsendCount = nd1;
  if(np>1)
  {
    YsendCount = (myRank==0 || myRank==np-1)? ny_l+2 : ny_l+1;
  }
  else
    YsendCount = ny_l+1+2*numGhost;
  Integer sendCount = XsendCount*YsendCount;
  Integer *recvCount_p = new Integer [np];
  MPI_Allgather(&sendCount, 1, MPI_INTEGER, recvCount_p, 1, MPI_INTEGER, MPI_COMM_WORLD);
  MPI_Allreduce(&sendCount, &NTot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
  Integer *displacement_p = new Integer [np];
  int sum=0;
  for(int i=0; i<np; i++)
  {
    displacement_p[i] = sum;
    sum += recvCount_p[i];
  }
  Integer Ystart = myRank==0   ? nd2a_l: n2a_l;
  Integer Yend   = myRank==np-1? nd2b_l: n2b_l;
  
  Real *uFinal_p = new Real [NTot];
  Real *errorG_p = new Real [NTot];
  Real *xloc_p[2];
  xloc_p[0] = new Real [NTot];
  xloc_p[1] = new Real [NTot];

  MPI_Gatherv(&uc(nd1a,Ystart), sendCount, MPI_DOUBLE, uFinal_p, 
              recvCount_p, displacement_p, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Gatherv(&err(nd1a, Ystart), sendCount, MPI_DOUBLE, errorG_p, 
              recvCount_p, displacement_p, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Gatherv(&x(nd1a, Ystart, 0), sendCount, MPI_DOUBLE, xloc_p[0], 
              recvCount_p, displacement_p, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Gatherv(&x(nd1a, Ystart, 1), sendCount, MPI_DOUBLE, xloc_p[1], 
              recvCount_p, displacement_p, MPI_DOUBLE, root, MPI_COMM_WORLD);
  
  if(saveMatlab>0 && myRank==root)
  {
    int t_ny, t_n2a, t_n2b;
    int f_ny, f_n2a, f_n2b;
    getLocalIndexBounds(np-1, np, ny, t_ny, t_n2a, t_n2b);
    getLocalIndexBounds(0   , np, ny, f_ny, f_n2a, f_n2b);
    FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
    fprintf(matlabFile, "%% File written by heat1Dmpi.C\n");
    fprintf(matlabFile,"xa=%g; xb=%g; ya=%g; yb=%g; kappa=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",xa,xb,ya,yb,kappa,tFinal,maxErr_g, cpuTimeStep_g);
    fprintf(matlabFile,"n1a=%d; n1b=%d; nd1a=%d; nd1b=%d;\n",n1a,n1b,nd1a,nd1b);
    fprintf(matlabFile,"n2a=%d; n2b=%d; nd2a=%d; nd2b=%d;\n",f_n2a,t_n2b,f_n2a-1,t_n2b+1);
    fprintf(matlabFile,"dx(1)=%14.6e; dx(2)=%14.6e; numGhost=%d;\n",dx[0],dx[1],numGhost);
    fprintf(matlabFile,"option=%d; optionName=\'%s\';\n",option,optionName.c_str());
    dimension(0, 1)= f_n2a-numGhost; dimension(1, 1)= t_n2b+numGhost;
    if( saveMatlab>1 )
    {
      writeMatlabArray(matlabFile,uFinal_p,"uFinal",dimension_p);
      writeMatlabArray(matlabFile,errorG_p,"errorG",dimension_p);
      writeMatlabArray(matlabFile,xloc_p[0],"x{1}",dimension_p);
      writeMatlabArray(matlabFile,xloc_p[1],"x{2}",dimension_p);
    }
    fclose(matlabFile);
    printf("Wrote file [%s]\n", matlabFileName.c_str());
  }

  MPI_Finalize();

  delete [] x_p[0];
  delete [] x_p[1];
  delete [] u_p[0];
  delete [] u_p[1];
  delete [] recvCount_p;
  delete [] displacement_p;
  delete [] err_p;
  delete [] uFinal_p;
  delete [] errorG_p;
  delete [] xloc_p[0];
  delete [] xloc_p[1];

  #undef x
  #undef gridIndexRange
  #undef dimension
  #undef boundaryCondition
  #undef uc
  #undef un
  #undef err

  return 0;

}
