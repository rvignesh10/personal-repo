#include <stdio.h>
#include "parseCommand.h"
typedef double Real;

__global__ void add( Real *a_d, Real *b_d, Real *c_d, int n, int Nt ){
	   int index = threadIdx.x + blockIdx.x*Nt;
	   if (index < n){
	      a_d[index] = -index;
	      b_d[index] = index*index;
	      c_d[index] = a_d[index] + b_d[index];
	   }
}

int main( int argc, char *argv[] ) {

    int n=1e+04;;
    int Nt=1;
    string line;
    bool echo = true;
    for( int i=1; i<argc; i++){
      line = argv[i];
      if (parseCommand(line, "-n=", n, echo)) {}
      else if(parseCommand(line, "-Nt=", Nt, echo)) {}
    }
    Real *a_p = new Real [n];
    Real *b_p = new Real [n];
    Real *c_p = new Real [n];

    // serial version
    for( int i=0; i<n; i++ ){
        a_p[i] = -i;
	b_p[i] = i*i;
	c_p[i] = a_p[i] + b_p[i];
    }

    // cuda version
    int Nb = ceil( (1.*n)/Nt);
    Real *a_d;
    Real *b_d;
    Real *c_d;
    int  *n_d;

    // cuda allocate memory
    cudaMalloc((void**)&a_d, n*sizeof(Real));
    cudaMalloc((void**)&b_d, n*sizeof(Real));
    cudaMalloc((void**)&c_d, n*sizeof(Real));
    cudaMalloc((void**)&n_d, sizeof(int));

    // cuda copy memory
    cudaMemcpy(a_d, a_p, n*sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(b_d, b_p, n*sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(c_d, c_p, n*sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(n_d, &n, sizeof(int), cudaMemcpyHostToDevice);

    add<<<Nb, Nt>>>(a_d, b_d, c_d, n, Nt);

    cudaMemcpy(a_p, a_d, n*sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(b_p, b_d, n*sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_p, c_d, n*sizeof(Real), cudaMemcpyDeviceToHost);

    return 0;
}