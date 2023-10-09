// ===========================================================================================
//
// POISSON EQUATION IN TWO DIMENSIONS
//          -Delta (u) = f
//
//       Solve with A++ arrays
//
// ===========================================================================================

#include "A++.h"

typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

// include commands tp parse command line arguments
#include "parseCommand.h"

// getCPU() : Return the current wall-clock time in seconds
#include "getCPU.h"

// enum for BC's
enum BoundaryConditionsEnum
{
  periodic=-1,
  dirichlet=1,
  neumann=2
};

// store parameters in this class
class PoissonParameters
{
public:
  Real dx[2];
  IntegerArray gridIndexRange;
  IntegerArray dimension;
  IntegerArray boundaryCondition;
  RealArray x;
  Real tol;
  RealArray uTrue;
  int maxIterations;
  int debug;
  int intervalToCheckResidual;
};


// -----------------------------------------------------------------------------
// Return the max-norm residual
// -----------------------------------------------------------------------------

Real getMaxResidual( RealArray &u, RealArray &f, PoissonParameters &par )
{
  const IntegerArray &gridIndexRange = par.gridIndexRange;
  const Real (&dx)[2]                   = par.dx;             // NOTE: reference to an array

  const Real dxSq = dx[0]*dx[0];
  const Real dySq = dx[1]*dx[1];

  const int n1a = gridIndexRange(0, 0), n1b = gridIndexRange(1, 0);
  const int n2a = gridIndexRange(0, 1), n2b = gridIndexRange(1, 1);

  Real maxRes = 0.;
  for ( int i2=n2a+1; i2<=n2b-1; i2++ )
  for ( int i1=n1a+1; i1<=n1b-1; i1++ )
  {
    Real res = f(i1, i2) + ( ( u(i1+1, i2) -2.*u(i1, i2) + u(i1-1, i2) )/dxSq +
		             ( u(i1, i2+1) -2.*u(i1, i2) + u(i1, i2-1) )/dySq );
    maxRes = max( maxRes, res );
  } 
  return maxRes;
}

// -----------------------------------------------------------------------------
// Return the max-norm error
// -----------------------------------------------------------------------------

Real getMaxError( RealArray &u, RealArray &err, PoissonParameters &par )
{
  const IntegerArray & gridIndexRange = par.gridIndexRange;
  const IntegerArray & dimension = par.dimension;
  const IntegerArray & boundaryCondition = par.boundaryCondition;
  const RealArray & x = par.x;
  const Real & tol = par.tol;
  const RealArray & uTrue = par.uTrue;
  const Real (&dx)[2] = par.dx; // Note: reference to an array
  const int & maxIterations = par.maxIterations;
  const int & debug = par.debug;

  const int n1a = gridIndexRange(0,0), n1b = gridIndexRange(1,0);
  const int n2a = gridIndexRange(0,1), n2b = gridIndexRange(1,1);

  Real maxErr=0.;
  for ( int i2=n2a; i2<=n2b; i2++ )
  for ( int i1=n1a; i1<=n1b; i1++ )
  {
    Real xi = x(i1, i2, 0);
    Real yi = x(i1, i2, 1);
    err(i1, i2) = fabs( u(i1, i2) - uTrue(i1, i2) );
    maxErr = max( err(i1, i2), maxErr );
  }

  if ( n1b-n1a+1 <=10 )
  {
    u.display("u");
    err.display("err");
  }
  return maxErr;
}



#define F(i1, i2) f_p[(i1) + nd1*(i2)]
#define U(i1, i2) u_p[(i1) + nd1*(i2)]
#define UN(i1, i2) un_p[(i1) + nd1*(i2)]

// ----------------------------------------------------------------------------
// Jacobi iteration
// ----------------------------------------------------------------------------

int jacobiIteration( RealArray &u, RealArray &f, PoissonParameters &par )
{
  const IntegerArray & gridIndexRange = par.gridIndexRange;
  const IntegerArray & dimension = par.dimension;
  const IntegerArray & boundaryCondition = par.boundaryCondition;
  const RealArray & x = par.x;
  const Real & tol = par.tol;
  const RealArray & uTrue = par.uTrue;
  const Real (&dx)[2] = par.dx; // Note: reference to an array
  const int & maxIterations = par.maxIterations;
  const int & debug = par.debug;
  const int & intervalToCheckResidual = par.intervalToCheckResidual;

  const int n1a = gridIndexRange(0,0), n1b = gridIndexRange(1,0);
  const int n2a = gridIndexRange(0,1), n2b = gridIndexRange(1,1);
  const int nd1a = dimension(0,0), nd1b = dimension(1,0);
  const int nd2a = dimension(0,1), nd2b = dimension(1,1);
  const int nd1 = nd1b-nd1a+1;

  const Real pi = M_PI;

  Range Rx(nd1a,nd1b), Ry(nd2a,nd2b);
  doubleArray ua[2];
  ua[0].redim(Rx, Ry); ua[0] = 0.;
  ua[1].redim(Rx, Ry); ua[1] = 0.;

  int current = 0;
  ua[current] = u; // initial guess

  int n;
  Real h=dx[0]; // assume dx=dy
  assert( dx[0] == dx[1] );
  Real omega = 1.;

  Real res0 = getMaxResidual( ua[current], f, par );
  Real maxRes = res0, maxResOld = res0;
  Real CR=1.;

  const Real *f_p = f.getDataPointer();

  Real cpu1 = getCPU();
  for ( n=0; n<maxIterations; n++ )
  {
    const int next = (current+1)%2;
    doubleArray &u = ua[current];
    doubleArray &un = ua[next];
    
    const Real *u_p = u.getDataPointer();
    Real *un_p = un.getDataPointer();
    
    // omega-Jacobi iterations
    for ( int i2=n2a+1; i2<=n2b-1; i2++ )
    for ( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
      Real z = .25*( h*h*F(i1, i2) + U(i1+1, i2) + U(i1-1, i2) + U(i1, i2+1) + U(i1, i2-1) );
      UN(i1, i2) = U(i1, i2) + omega*( z - U(i1, i2) );
    }
    current = next;
    
    if ( (n% intervalToCheckResidual == 0) || (n == maxIterations-1) )
    {
      // check for convergence
      maxRes = getMaxResidual( u, f, par );
      CR = pow( maxRes/maxResOld, 1./intervalToCheckResidual );
      maxResOld = maxRes;

     if (false)
	printf("n=%6d, maxRes=%9.3e\n",n,maxRes);
     if (maxRes < tol )
	break; 
    }

  }
  Real cpu = getCPU() - cpu1;

  const int numIterations=n;
  // printf("numIterations=%d, nx=%d, ny=%d, cpu time= %9.2e (s)\n",numIterations,nx,ny,cpu);
  if( maxRes<tol )
  {
    // printf("CONVERGENCE: maxRes=%8.2e < tol=%9.2e.\n",maxRes,tol);
  }
  else
    printf("Jacobi: ERROR: maxRes=%8.2e > tol=%9.2e : NO CONVERGENCE\n",maxRes,tol);

  // --- compute errors ---
  u = ua[current];
  RealArray err(Rx,Ry);

  Real maxErr = getMaxError( u, err, par );
  maxRes = getMaxResidual( u,f,par );

  // Average convergence rate:
  Real aveCR = pow( (maxRes/res0), 1./max(1,numIterations) );
  // Asymptotic convergence rate:
  Real ACR = fabs( 1. - omega*(1.-cos(pi*h)) );
  printf("Jac: omega=%4.2f Its=%6d res=%8.2e err=%8.2e CR=%7.5f aveCR=%7.5f ACR=%7.5f cpu=%7.1es cpu/it=%7.1es\n",
         omega,numIterations,maxRes,maxErr,CR,aveCR,ACR,cpu,cpu/numIterations);

  return numIterations;

}

// ---------------------------------------------
// Gauss-Seidel or SOR
// ---------------------------------------------
int gaussSeidelIteration( RealArray & u, RealArray & f, PoissonParameters & par )
{
  const IntegerArray & gridIndexRange = par.gridIndexRange;
  const IntegerArray & dimension = par.dimension;
  const IntegerArray & boundaryCondition = par.boundaryCondition;
  const RealArray & x = par.x;
  const Real & tol = par.tol;

  const RealArray & uTrue = par.uTrue;
  const Real (&dx)[2] = par.dx; // Note: reference to an array
  const int & maxIterations = par.maxIterations;
  const int & debug = par.debug;
  const int & intervalToCheckResidual = par.intervalToCheckResidual;

  const int n1a = gridIndexRange(0,0), n1b = gridIndexRange(1,0);
  const int n2a = gridIndexRange(0,1), n2b = gridIndexRange(1,1);
  const int nd1a = dimension(0,0), nd1b = dimension(1,0);
  const int nd2a = dimension(0,1), nd2b = dimension(1,1);
  const int nd1 = nd1b-nd1a+1;
  Range Rx(nd1a,nd1b), Ry(nd2a,nd2b);

  const Real pi = M_PI;

  int n;
  Real h=dx[0]; // assume dx=dy
  assert( dx[0]==dx[1] );

  Real omega = 2./(1.+sin(pi*h)); // optimal omega

  Real res0 = getMaxResidual( u,f,par );
  Real maxRes=res0, maxResOld=maxRes, CR=1.;

  const Real *f_p = f.getDataPointer();
  Real *u_p = u.getDataPointer();

  Real cpu1 = getCPU();
  for( n=0; n<maxIterations; n++ )
  {
    // omega-GS Iteration
    for( int i2=n2a+1; i2<=n2b-1; i2++ )
    for( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
      Real z = .25*( h*h*F(i1,i2) + U(i1+1,i2) + U(i1-1,i2) + U(i1,i2+1) + U(i1,i2-1) );
                                    U(i1,i2) = U(i1,i2) + omega*( z - U(i1,i2) );
    }

    if( ( n % intervalToCheckResidual) ==0 || n==(maxIterations-1) )
    {
    // check for convergence
    maxRes = getMaxResidual( u,f,par );

    CR = pow( (maxRes/maxResOld), 1./intervalToCheckResidual );
    maxResOld=maxRes;
    if( false )
    {
      printf("GS: n=%6d, maxRes=%9.3e, CR=%7.4f\n",n,maxRes,CR);
    }

      if( maxRes<tol )
        break;
    }

  }
  Real cpu = getCPU()-cpu1;

  const int numIterations=n;
  if ( maxRes < tol )
  {
  }
  else
    printf("GS: ERROR: maxRes=%8.2e > tol=%9.2e : NO CONVERGENCE\n",maxRes,tol);
  
  // --- compute errors ---
  RealArray err(Rx,Ry);
  Real maxErr = getMaxError( u, err, par );
  maxRes = getMaxResidual( u,f,par );

  // Average convergence rate:
  Real aveCR = pow( (maxRes/res0), 1./max(1,numIterations) );
  // Asymptotic convergence rate:
  Real ACR = omega-1.;
  printf("GS : omega=%4.2f Its=%6d res=%8.2e err=%8.2e CR=%7.5f aveCR=%7.5f ACR=%7.5f cpu=%7.1es cpu/it=%7.1es\n",
		omega,numIterations,maxRes,maxErr,CR,aveCR,ACR,cpu,cpu/numIterations);

  return numIterations;

}

// ---------------------------------------------
// Red-Black Gauss-Seidel iteration
// ---------------------------------------------
// int redBlackIteration( RealArray & u, RealArray & f, Real xa, Real ya, Real dx[], int
// maxIterations, Real tol, RealArray & uTrue )
int redBlackIteration( RealArray & u, RealArray & f, PoissonParameters & par )
{

  const IntegerArray & gridIndexRange = par.gridIndexRange;
  const IntegerArray & dimension = par.dimension;
  const IntegerArray & boundaryCondition = par.boundaryCondition;
  const RealArray & x = par.x;
  const Real & tol = par.tol;
  const RealArray & uTrue = par.uTrue;
  const Real (&dx)[2] = par.dx; // Note: reference to an array
  const int & maxIterations = par.maxIterations;
  const int & debug = par.debug;
  const int & intervalToCheckResidual = par.intervalToCheckResidual;

  const int n1a = gridIndexRange(0,0), n1b = gridIndexRange(1,0);
  const int n2a = gridIndexRange(0,1), n2b = gridIndexRange(1,1);
  const int nd1a = dimension(0,0), nd1b = dimension(1,0);
  const int nd2a = dimension(0,1), nd2b = dimension(1,1);
  const int nd1 = nd1b-nd1a+1;
  Range Rx(nd1a,nd1b), Ry(nd2a,nd2b);

  const Real pi = M_PI;

  doubleArray ua[2];
  ua[0].redim(Rx,Ry); ua[0]=0.;
  ua[1].redim(Rx,Ry); ua[1]=0.;

  int current = 0;
  ua[0] = u; // initial guess
  ua[1] = u; // initial guess

  int n;
  Real h=dx[0]; // assume dx=dy
  assert( dx[0]==dx[1] );

  // this is the optimal omega for ANY ordering (Owl book p32)
  const Real omega = 2./(1.+sin(pi*h));

  Real res0 = getMaxResidual( ua[current],f,par );
  Real maxRes=res0, maxResOld=res0;
  Real CR=1.;

  const Real *f_p = f.getDataPointer();

  Real cpu1 = getCPU();
  for( n=0; n<maxIterations; n++ )
  {
    const int next = (current+1) % 2;
    doubleArray & u = ua[current];
    doubleArray & un = ua[next];

    Real *u_p = u.getDataPointer();
    Real *un_p = un.getDataPointer();

    // --- omega Red-Black ---
    // Question : what is the fastest way to loop over red and black points?
    // (1) Use if statement
    // (2) Use loops with stride of 2
    // Note: loops are faster if LHS array is different from RHS arrays
    for( int i2=n2a+1; i2<=n2b-1; i2++ )
    for( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
      if( (i1+i2)%2 == 0 ) // red points
      {
        Real z = .25*( h*h*F(i1,i2) + U(i1+1,i2) + U(i1-1,i2) + U(i1,i2+1) + U(i1,i2-1) );
                                      UN(i1,i2) = U(i1,i2) + omega*( z - U(i1,i2) );
      }
      else
      {
        UN(i1,i2) = U(i1,i2);
      }
    }

    for( int i2=n2a+1; i2<=n2b-1; i2++ )
    for( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
      if( (i1+i2)%2 == 1 ) // black points
      {
        Real z = .25*( h*h*F(i1,i2) + UN(i1+1,i2) + UN(i1-1,i2) + UN(i1,i2+1) + UN(i1,i2-1) );
                                      U(i1,i2) = UN(i1,i2) + omega*( z - UN(i1,i2) );
      }
      else
      {
        U(i1,i2) = UN(i1,i2);
      }
    }
    
    if ( (n%intervalToCheckResidual == 0) || ( n == maxIterations-1) )
    {
      // check for convergence
      maxRes = getMaxResidual( u,f,par );
      CR = pow( (maxRes/maxResOld), 1./intervalToCheckResidual );
      maxResOld=maxRes;

      if( false )
        printf("RB: n=%6d, maxRes=%9.3e, CR=%7.5f\n",n,maxRes,CR);

      if( maxRes<tol )
        break;
    }
  }
  Real cpu = getCPU()-cpu1;

  const int numIterations=n;
  // printf("numIterations=%d, nx=%d, ny=%d, cpu time= %9.2e (s)\n",numIterations,nx,ny,cpu);
  if( maxRes<tol )
  {
    // printf("CONVERGENCE: maxRes=%8.2e < tol=%9.2e.\n",maxRes,tol);
  }
  else
     printf("Red-Black: ERROR: maxRes=%8.2e > tol=%9.2e : NO CONVERGENCE\n",maxRes,tol);


  // --- compute errors ---
  u = ua[current];
  RealArray err(Rx,Ry);
  Real maxErr = getMaxError( u, err, par );
  maxRes = getMaxResidual( u,f,par );

  // Average convergence rate:
  Real aveCR = pow( (maxRes/res0), 1./max(1,numIterations) );
  // Asymptotic convergence rate:
  Real ACR = omega-1.;
  printf("RB : omega=%4.2f Its=%6d res=%8.2e err=%8.2e CR=%7.5f aveCR=%7.5f ACR=%7.5f cpu=%7.1es cpu/it=%7.1es\n",
               omega,numIterations,maxRes,maxErr,CR,aveCR,ACR,cpu,cpu/numIterations);

  return numIterations;

}



// ---------------------------------------------
// Conjugate Gradient iteration
// ---------------------------------------------
int conjugateGradientIteration( RealArray & u, RealArray & f, PoissonParameters & par )
{
  const IntegerArray & gridIndexRange = par.gridIndexRange;
  const IntegerArray & dimension = par.dimension;
  const IntegerArray & boundaryCondition = par.boundaryCondition;
  const RealArray & x = par.x;
  const Real & tol = par.tol;
  const RealArray & uTrue = par.uTrue;
  const Real (&dx)[2] = par.dx; // Note: reference to an array
  const int & maxIterations = par.maxIterations;
  const int & debug = par.debug;
  const int & intervalToCheckResidual = par.intervalToCheckResidual;
  
  const int n1a = gridIndexRange(0,0), n1b = gridIndexRange(1,0);
  const int n2a = gridIndexRange(0,1), n2b = gridIndexRange(1,1);
  const int nd1a = dimension(0,0), nd1b = dimension(1,0);
  const int nd2a = dimension(0,1), nd2b = dimension(1,1);
  const int nd1 = nd1b-nd1a+1;

  const Real dxSq = dx[0]*dx[0];
  const Real dySq = dx[1]*dx[1];

  Range Rx(nd1a,nd1b), Ry(nd2a,nd2b);
 
  // Temporary arrays
  RealArray z(Rx,Ry), r(Rx,Ry), p(Rx,Ry);

  const Real *f_p = f.getDataPointer();

  Real *u_p = u.getDataPointer();
  Real *z_p = z.getDataPointer();
  Real *r_p = r.getDataPointer();
  Real *p_p = p.getDataPointer();
  #define Z(i1,i2) z_p[(i1) + nd1*(i2)]
  #define R(i1,i2) r_p[(i1) + nd1*(i2)]
  #define P(i1,i2) p_p[(i1) + nd1*(i2)]

  Real cpu1 = getCPU();
  
  u = 0.; // initial guess -- *** FIX ME for non-zero initial guess ***
  r = 0; // initial residual goes here
  for( int i2=n2a+1; i2<=n2b-1; i2++ )
  for( int i1=n1a+1; i1<=n1b-1; i1++ )
  {
    R(i1,i2) = F(i1,i2); // initial residual
    P(i1,i2) = R(i1,i2); // initial search direction
  }

  Real alpha,beta, rNormSquared, rNormSquaredNew, pz;

  // rNormSquared = r^T r
  rNormSquared=0.;
  for( int i2=n2a+1; i2<=n2b-1; i2++ )
  for( int i1=n1a+1; i1<=n1b-1; i1++ )
  {
    rNormSquared += R(i1,i2)*R(i1,i2);
  }

  Real res0 = sqrt(rNormSquared);
  Real maxRes=res0, maxResOld=res0;
  Real CR=1.;
  
  int n, nOld = -1;
  for ( n=0; n<maxIterations; n++ )
  {
    // NOTE: Some loops have been combined for speed
    // z  = A p
    // pz = p^T z
    // Note: A = -Delta
    pz = 0.;
    for( int i2=n2a+1; i2<=n2b-1; i2++ )
    for( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
      Z(i1,i2) = -( (P(i1+1,i2) -2.*P(i1,i2) + P(i1-1,i2))/dxSq
                + (P(i1,i2+1) -2.*P(i1,i2) + P(i1,i2-1))/dySq );
      pz += P(i1,i2)*Z(i1,i2);
    }

    alpha = rNormSquared/pz; // step length
    // Merging these into one loop is faster in serial
    // u += alpha*p
    // r -= alpha*z
    // rNormSquaredNew = r^T r
    rNormSquaredNew = 0.;
    for( int i2=n2a+1; i2<=n2b-1; i2++ )
    for( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
      U(i1,i2) += alpha*P(i1,i2); // new solution
      R(i1,i2) -= alpha*Z(i1,i2); // residual
      rNormSquaredNew += R(i1,i2)*R(i1,i2);
    }

    beta = rNormSquaredNew/rNormSquared; // improvement this step
    rNormSquared=rNormSquaredNew;

    // p = r + beta*p
    for( int i2=n2a+1; i2<=n2b-1; i2++ )
    for( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
      P(i1,i2) = R(i1,i2) + beta*P(i1,i2); // new search direction
    }

    if( ( n % intervalToCheckResidual) ==0 || n==(maxIterations-1) )
    {
      // check for convergence **** WE SHOULD REALLY USE EXISTING 2-norm RESIDUAL ****
      // Do this for now to be consistent with other schemes.
      maxRes = getMaxResidual( u,f,par );
      CR = pow( (maxRes/maxResOld), 1./max(1,n-nOld) );
      maxResOld=maxRes;
      nOld=n;

      if (debug > 0 )
	printf("CG: n=%3d: pz=%9.2e, alpha=%9.2e, beta=%9.2e, maxRes=%8.2e, CR=%7.5f \n",n,pz,alpha,beta,maxRes,CR);
      if ( maxRes < tol )
	break;
    }

  }
  Real cpu = getCPU() - cpu1;
  
  const int numIterations=n;
  if ( maxRes < tol )
  {
  }
  else
    printf("CG: ERROR: maxRes=%8.2e > tol=%9.2e : NO CONVERGENCE\n",maxRes,tol);

  // --- compute errors ---
  RealArray err(Rx,Ry);
  Real maxErr = getMaxError( u, err, par );
  maxRes = getMaxResidual( u,f,par );

  // Average convergence rate:
  Real aveCR = pow( (maxRes/res0), 1./max(1,numIterations) );
  // Asymptotic convergence rate:
  const Real pi = M_PI;
  Real ACR = 1. - pi*dx[0]; // ** CHECK ME **
  printf("CG :            Its=%6d res=%8.2e err=%8.2e CR=%7.5f aveCR=%7.5f ACR=%7.5f cpu=%7.1es cpu/it=%7.1es\n",
		  numIterations,maxRes,maxErr,CR,aveCR,ACR,cpu,cpu/numIterations);

  return numIterations;

}


int main( int argc, char *argv[] )
{
  printf("Usage: poisson -nx=<i> -tol=<f> -maxIterations=<i> -debug=<i>\n"
		  "    nx = number of grid points in x and y directions\n"
		  "    tol = convergence tolerance\n"
		  "    maxIterations = max number of iterations\n");
  const Real pi = M_PI; // 4.*atan2(1., 1.);

  // Parameters are stored here:
  PoissonParameters par;

  // Make references to parameters for clarity
  IntegerArray & gridIndexRange = par.gridIndexRange;
  IntegerArray & dimension = par.dimension;
  IntegerArray & boundaryCondition = par.boundaryCondition;
  RealArray & x = par.x;
  Real & tol = par.tol;
  RealArray & uTrue = par.uTrue;
  Real (&dx)[2] = par.dx; // Note: reference to an array
  int & maxIterations = par.maxIterations;
  int & debug = par.debug;
  int & intervalToCheckResidual = par.intervalToCheckResidual;

  intervalToCheckResidual=50; // check residual every this many iterations

  const int numberOfDimensions=2;
  Real xa=0., xb=1.; // domain is [xa,xb] X [ya,yb]
  Real ya=0., yb=1.;

  debug = 0;
  maxIterations = 1000;
  tol=1.e-3;
  int nx=100, ny=nx;

  string line;
  for( int i=1; i<argc; i++ )
  {
    line=argv[i];
    // printf("Input: argv[%d] = [%s]\n",i,line.c_str());
    if( parseCommand( line,"-nx=",nx) )
    {
      ny=nx;
    }
    else if( parseCommand( line,"-debug=",debug) ){}
    else if( parseCommand( line,"-maxIterations=",maxIterations) ){}
    else if( parseCommand( line, "-tol=",tol) ){}

  }
  printf("-----------------------------------------------------------------\n");
  printf("--------- Solve the Poisson Equation in two dimensions ----------\n");
  printf("              - Delta (u) = f                                    \n");
  printf("              Dirichlet BCs                                      \n");
  printf("   nx=%d, ny=%d, maxIterations=%d, tol=%9.2e\n",nx,ny,maxIterations,tol);
  printf("-----------------------------------------------------------------\n");

  const Real kx=1., ky=1.;
  const Real kxp=kx*pi;
  const Real kyp=ky*pi;

  // This true solution is too easy for CG since u is a multiple of f:
  // #define UTRUE(x,y) sin( kxp*(x) )*sin( kyp*(y) )
  // f = - Delta (u)
  // #define FORCE(x,y) (kxp*kxp +kyp*kyp)*UTRUE(x,y)

  #define UTRUE(x,y) (x)*(1.-(x)) * (y)*(1.-(y))
  // f = - Delta (u)
  #define FORCE(x,y) 2.*( (x)*(1.-(x)) + (y)*(1.-(y)) )

  const int numGhost=0;
  const int n1a = numGhost;
  const int n1b = n1a + nx;
  const int nd1a = n1a-numGhost;
  const int nd1b = n1b+numGhost;
  const int nd1 = nd1b-nd1a+1;

  const int n2a = numGhost;
  const int n2b = n2a + ny;
  const int nd2a = n2a-numGhost;
  const int nd2b = n2b+numGhost;
  const int nd2 = nd2b-nd2a+1;

  gridIndexRange.redim(2,numberOfDimensions);
  dimension.redim(2,numberOfDimensions);
  boundaryCondition.redim(2,numberOfDimensions);

  gridIndexRange(0,0)=n1a; gridIndexRange(1,0)=n1b;
  gridIndexRange(0,1)=n2a; gridIndexRange(1,1)=n2b;

  dimension(0,0)=nd1a; dimension(1,0)=nd1b;
  dimension(0,1)=nd2a; dimension(1,1)=nd2b;

  boundaryCondition(0,0)=dirichlet; // left
  boundaryCondition(1,0)=dirichlet; // right

  boundaryCondition(0,1)=dirichlet; // bottom
  boundaryCondition(1,1)=dirichlet; // top

  // Grid points
  Range Rx(nd1a,nd1b), Ry(nd2a,nd2b);
  x.redim(Rx,Ry,2);

  // Real dx[2];
  dx[0] = (xb-xa)/nx;
  dx[1] = (yb-ya)/ny;

  int i1,i2;
  for( i2=nd2a; i2<=nd2b; i2++ )
  for( i1=nd1a; i1<=nd1b; i1++ )
  {
    x(i1,i2,0) = xa + (i1-n1a)*dx[0];
    x(i1,i2,1) = ya + (i2-n2a)*dx[1];
  }


  RealArray f(Rx,Ry);
  uTrue.redim(Rx,Ry);

  RealArray u(Rx,Ry);

  Real xi,yi;
  
  for( int i2=nd2a; i2<=nd2b; i2++ )
  for( int i1=nd1a; i1<=nd1b; i1++ )
  {
    xi = x(i1,i2,0);
    yi = x(i1,i2,1);
    f(i1,i2)= FORCE(xi,yi);
    uTrue(i1,i2) = UTRUE(xi,yi);
  }

  // ===== JACOBI =====
  u=0.; // initial guess
  jacobiIteration( u, f, par );

  // ===== Gauss-Seidel =====
  u=0.; // initial guess
  gaussSeidelIteration( u, f, par );
 
  // ===== Red-Black Gauss-Seidel =====
  u=0.; // initial guess
  redBlackIteration( u, f, par );

  // ==== Conjugate Gradient =========
  u=0.; // initial guess
  conjugateGradientIteration( u,f,par );

  return 0;

}



