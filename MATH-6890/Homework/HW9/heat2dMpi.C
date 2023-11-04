// =====================================================================
//
//     HEAT EQUATION IN TWO DIMENSIONS
//
// =====================================================================
#define USE_PPP

#include <mpi.h>

typedef double Real;
typedef int Integer;

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
    #ifndef USE_PPP
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

