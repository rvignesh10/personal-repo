//#ifndef GET_CPU_H
//#define GET_CPU_H

//#include <ctime>

// ---------------------------------------------
// Return the current wall-clock time in seconds
// ---------------------------------------------
//inline double getCPU() {
//    return ( 1.0*std::clock() )/CLOCKS_PER_SEC ;
//}

#ifndef GET_CPU_H
#define GET_CPU_H

#include <ctime>
// ---------------------------------------------
// Return the current wall-clock time in seconds
// ---------------------------------------------
inline double getCPU()
{
  #if defined(_OPENMP)
    return omp_get_wtime();
  #else
    return ( 1.0*std::clock() )/CLOCKS_PER_SEC ;
  #endif
}

#endif

//#endif
