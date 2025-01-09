#include <stdio.h>

#include <ctime>
// ------------------------------------------------------------------------
// Return the current wall-clock time in seconds
// ------------------------------------------------------------------------
inline double getCPU(){
    return ( 1.0 * std::clock() )/CLOCKS_PER_SEC ;
}

#include "parseCommand.h"

__global__ void add( int *a_d, int *b_d, int *c_d, int *n_d, int *nt_d ){
  int idx= threadIdx.x + blockIdx.x*(*nt_d);
  if(idx < *n_d)
    c_d[idx] = a_d[idx] + b_d[idx];
}

int main( int argc, char *argv[] ) {

  int n=1000;

  std::string line;
  bool echo=true;
  for( int i=1; i<argc; i++ )
  {
    line=argv[i];
    if( parseCommand( line, "-n=", n, echo ) ){}
  }

  int *a_p = new int [n];
  int *b_p = new int [n];
  int *c_p = new int [n];

  double cpu0 = getCPU();
  for( int i=0; i<n; i++ ) {
    a_p[i] = -i; b_p[i] = i*i;
    c_p[i] = a_p[i] + b_p[i];
  }
  double cpuTime= getCPU() - cpu0;
  int *a_d;
  int *b_d;
  int *c_d;
  int *n_d;
  int *nt_d;
  cudaMalloc((void **)&c_d, n*sizeof(int));
  cudaMalloc((void **)&a_d, n*sizeof(int));
  cudaMalloc((void **)&b_d, n*sizeof(int));
  cudaMalloc((void **)&n_d, sizeof(int));
  cudaMalloc((void **)&nt_d,sizeof(int));
  cudaMemcpy(a_d, a_p, n*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(b_d, b_p, n*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(n_d, &n ,   sizeof(int), cudaMemcpyHostToDevice);
  //cudaMemcpy(nt_d,&Nt,   sizeof(int), cudaMemcpyHostToDevice);
  int Nt;
  for( int k=0; k<=10; k++ ){
    Nt=(int)(pow(2, k));
    int Nb = ceil( (1.*n)/Nt );
    cudaMemcpy(nt_d,&Nt,   sizeof(int), cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);
    add<<<Nb, Nt>>>(a_d, b_d, c_d, n_d, nt_d);
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    float gpuTime;
    cudaEventElapsedTime(&gpuTime, start, stop);

    cudaMemcpy(c_p, c_d, n*sizeof(int), cudaMemcpyDeviceToHost);

    double maxErr=0.0;
    for(int i=0; i<n; i++){
      double t = fabs(c_p[i]- a_p[i]-b_p[i]);
      maxErr = t>=maxErr? t:maxErr;
    }
    // printf("maxErr=%f\n",maxErr); 
   printf("%10d \t %10d \t %10d \t %10d \t %25.15e \t %25.15e \t %f\n", Nb, Nt, Nb*Nt, n, maxErr, gpuTime, cpuTime*1000/gpuTime);   
  }
  cudaFree(a_d); cudaFree(b_d);
  cudaFree(c_d); cudaFree(n_d); cudaFree(nt_d);
  delete [] a_p;
  delete [] b_p;
  delete [] c_p;
  
  return 0;
}