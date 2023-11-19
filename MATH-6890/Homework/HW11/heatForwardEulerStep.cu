/* multiply.cu */

#include <cuda.h>
#include <cuda_runtime.h>

const double b0 = 1.0, b1 = 0.5, b2 = 0.25;
const double a0 = 1.0, a1 = 0.3;

#define UTRUE(x, t) ( b0 + (x) * (b1 + (x) * b2) ) * ( a0 + (t) * a1 )
#define UTRUEX(x, t) ( b1 + 2. * (x) * b2 ) * ( a0 + (t) * ( a1 ) )
#define UTRUET(x, t) (b0 + (x) * ( b1 + (x) * b2 )) * ( a1 )
#define UTRUEXX(x, t) ( 2.*b2 )*( a0 + (t)*( a1 ) )
#define FORCE(x, t) ( UTRUET(x, t) - kappa*UTRUEXX(x, t) )

__global__ void __setSpatialMeshNodes1D__(double *dev_x, double *dev_xa, double *dev_dx, int *dev_nd1L, int *dev_Nt)
{
    const int idx = threadIdx.x + blockIdx.x * (*dev_Nt);
    if( idx < *dev_nd1L )
        dev_x[idx] = *dev_xa + (idx)*(*dev_dx);
}

__global__ void __heat1dSetInitialCondition__ (double *dev_uc, double *dev_x, int *dev_Nt, int *dev_nd1L) 
{
    const int idx = threadIdx.x + blockIdx.x * (*dev_Nt);
    if( idx < *dev_nd1L )
        dev_uc[idx] = UTRUE( dev_x[idx], 0. );
}

__global__ void __heat1dForwardEulerUpdate__ (double *dev_uc, double *dev_un, double *dev_x, double *dev_rx, double *dev_t, int *dev_Nt, int *dev_n1aL, int *dev_n1bL, int *dev_numGhost)
{
    // const int idx = threadIdx.x + blockIdx.x * blockDim.x;
    const int idx = threadIdx.x + blockIdx.x * (*dev_Nt);
    // check if idx is within the internal range
    if( (idx>= *dev_n1aL + *dev_numGhost) || (idx <= *dev_n1bL + *dev_numGhost) )
        dev_un[idx] = (*dev_rx)*(dev_uc[idx+1] -2.*dev_uc[idx] + dev_uc[idx-1]) + FORCE( dev_x[idx], *dev_t ) ;
}

__global__ void __heat1dErrorCalc__(double *dev_err, double *dev_uc, double *dev_x, double *dev_t, int *dev_nd1L, int *dev_Nt)
{
    const int idx = threadIdx.x + blockIdx.x * (*dev_Nt);
    if(idx < *dev_nd1L)
        dev_err[idx] = dev_uc[idx] - UTRUE( dev_x[idx], *dev_t );
}

extern "C" void setSpatialMeshNodes1D( double *dev_x, double *dev_xa, double *dev_dx, int *dev_nd1L, int Nb, int Nt, int *dev_Nt )
{
    __setSpatialMeshNodes1D__<<<Nb, Nt>>>(dev_x, dev_xa, dev_dx, dev_nd1L, dev_Nt);
    safecall(cudaThreadSynchronize());
    safecall(cudaGetLastError());
}

extern "C" void heat1dSetInitialCondition(double *dev_uc, double *dev_x, int Nb, int Nt, int *dev_Nt, int *dev_nd1L)
{   
    __heat1dSetInitialCondition__<<<Nb, Nt>>> (dev_uc, dev_x, dev_Nt, dev_nd1L);
    safecall(cudaThreadSynchronize());
    safecall(cudaGetLastError());
}

extern "C" void heat1dForwardEulerUpdate(double *dev_uc, double *dev_un, double *dev_x, double *dev_rx, double *dev_t, 
                                         int Nb, int Nt, int *dev_Nt, int *dev_n1aL, int *dev_n1bL, int *dev_numGhost, float *gpuStepTime)
{   
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);
    __heat1dForwardEulerUpdate__ <<<Nb, Nt>>> ( dev_uc, dev_un, dev_x, dev_rx, dev_t, dev_Nt, dev_n1aL, dev_n1bL, dev_numGhost );
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(gpuStepTime, start, stop);

    safecall(cudaThreadSynchronize());
    safecall(cudaGetLastError());
}

extern "C" void heat1dErrorCalc(double *dev_err, double *dev_uc, double *dev_x, double *dev_t, int Nb, int Nt, int *dev_Nt, int *dev_nd1L)
{
    __heat1dErrorCalc__<<<Nb, Nt>>>(dev_err, dev_uc, dev_x, dev_t, dev_nd1L, dev_Nt);
    safecall(cudaThreadSynchronize());
    safecall(cudaGetLastError());
}

/* --------------------------------------------------------------------------------------------------- */
// wrapper functions for allocating and free-ing cuda memory 
/* --------------------------------------------------------------------------------------------------- */
extern "C" void AllocateCudaMemory( double *dev_var, int n_size ) { cudaMalloc((void **)&dev_var, (n_size)*sizeof(double)); }
extern "C" void AllocateCudaMemory( int *dev_var, int n_size ) { cudaMalloc((void **)&dev_var, (n_size)*sizeof(int)); }

extern "C" void FreeCudaMemory( double *dev_var ) { cudaFree(dev_var); }
extern "C" void FreeCudaMemory( int *dev_var ) { cudaFree(dev_var); }

/* --------------------------------------------------------------------------------------------------- */
// wrapper functions for copying memory
/* --------------------------------------------------------------------------------------------------- */
extern "C" void MemcpyHostToDev( double *dev_var, double *hst_var, int n_size ) { cudaMemcpy( dev_var, hst_var, (n_size)*sizeof(double), cudaMemcpyHostToDevice ); }
extern "C" void MemcpyHostToDev( int *dev_var, int *hst_var, int n_size ) { cudaMemcpy( dev_var, hst_var, (n_size)*sizeof(int), cudaMemcpyHostToDevice ); }

extern "C" void MemcpyDevToHost( double *hst_var, double *dev_var, int n_size ) { cudaMemcpy( hst_var, dev_var, (n_size)*sizeof(double), cudaMemcpyDeviceToHost ); }
extern "C" void MemcpyDevToHost( int *hst_var, int *dev_var, int n_size ) { cudaMemcpy( hst_var, dev_var, (n_size)*sizeof(int), cudaMemcpyDeviceToHost ); }
