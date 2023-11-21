/*  ------------------------- Solving Heat Equation in 1D ------------------------- */

// ------------------------------------------------------------------------
// ways to execute this code 
// ./heat1d -Nx #discretizationSize -matlabFileName #fileName.m
// examples:
// ./heat1d -Nx 40 -matlabFileName heat1dBcDDNx20.m
// ------------------------------------------------------------------------

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <float.h>
#include <assert.h>

// define a new type "Real" which is equivalent to double
typedef double Real;
#define MAX_NUM_THREADS 128

#include <string>
using std::string;
using std::max;

#include <ctime>


// some necessary constants
Real kappa  = 0.1;
const Real kx = 3.0;
const Real kxpi = kx * M_PI;
const Real kappaPiSq = kappa * kxpi * kxpi;

const Real b0 = 1.0, b1 = 0.5, b2 = 0.25;
const Real a0 = 1.0, a1 = 0.3;

// defining different type of solutions available in the code
#define TRIG_DD 1
#define TRIG_NN 2
#define POLY_DD 3
#define POLY_NN 4

#ifndef SOLUTION
    #define SOLUTION TRIG_DD
    // #define SOLUTION TRIG_NN
    // #define SOLUTION POLY_DD
    // #define SOLUTION POLY_NN
#endif

// assign boundary conditions for different problems
#if SOLUTION == TRIG_DD
    const char solutionName[] = "trueDD";

    #define UTRUE(x, t) sin( kxpi * (x) ) * exp( -kappaPiSq * (t) )
    #define UTRUEX(x, t) kxpi * cos( kxpi * (x) ) * exp( -kappaPiSq * (t) )
    #define UTRUET(x, t) -kappaPiSq * UTRUE(x, t)
    #define FORCE(x, t) (0.)

#elif SOLUTION == TRIG_NN
    const char solutionName[] = "trueNN";

    #define UTRUE(x, t) cos( kxpi * (x) ) * exp( -kappaPiSq * (t) )
    #define UTRUEX(x, t)  -kxpi * sin( kxpi * (x) ) * exp( -kappaPiSq * (t) )
    #define UTRUET(x, t) -kappaPiSq * UTRUE(x, t)
    #define FORCE(x, t) (0.)

#elif SOLUTION == POLY_DD || SOLUTION == POLY_NN
    #if SOLUTION == POLY_DD
        const char solutionName[] = "polyDD";
    #else
        const char solutionName[] = "polyNN";
    #endif

    #define UTRUE(x, t) ( b0 + (x) * (b1 + (x) * b2) ) * ( a0 + (t) * a1 )
    #define UTRUEX(x, t) ( b1 + 2. * (x) * b2 ) * ( a0 + (t) * ( a1 ) )
    #define UTRUET(x, t) (b0 + (x) * ( b1 + (x) * b2 )) * ( a1 )
    #define UTRUEXX(x, t) ( 2.*b2 )*( a0 + (t)*( a1 ) )
    #define FORCE(x, t) ( UTRUET(x, t) - kappa*UTRUEXX(x, t) )
#else
    const char solutionName[] = "wrongChoice";
#endif

/* --------------------------------------------------------------------- */
//                              cuda kernels
/* --------------------------------------------------------------------- */

__global__ void mesh1d(double *x_d, int *nd1_d, int *nd1a_d, int *nd1b_d, double *xa_d, double *dx_d, int *nt_d){
    const int idx = threadIdx.x + blockIdx.x * (*nt_d) + *nd1a_d ;
    #define x(i1)  x_d[i1 - *nd1a_d]
    if (idx >= *nd1a_d && idx <= *nd1b_d)
        x(idx) = *xa_d + idx*(*dx_d);
    #undef x
}

__global__ void setInitialCondition(double *u_d, int *nd1_d, int *nd1a_d, int *nd1b_d, double *x_d, double *t_d, int *nt_d){
    const int idx = threadIdx.x + blockIdx.x * (*nt_d) + *nd1a_d ;
    #define u(i1) u_d[i1 - *nd1a_d]
    #define x(i1) x_d[i1 - *nd1a_d]

    if(idx >= *nd1a_d && idx <= nd1b_d)
        u(idx) = UTRUE( x(idx), *t_d );

    #undef x
    #undef u
}

__global__ void heat1dForwardEulerTimeStep( double *uc_d, double *un_d, double *x_d,
                                            int *nd1_d, int *nd1a_d, int *nd1b_d, 
                                            double *rx_d, double *dt_d, double *t_d, int *nt_d ){
    const int idx = threadIdx.x + blockIdx.x * (*nt_d) + *nd1a_d ;

    #define uc(i1) uc_d[i1- *nd1a_d]
    #define un(i1) un_d[i1- *nd1a_d]
    #define x(i1)  x_d[i1 - *nd1a_d]

    if(idx > *nd1a_d && idx < *nd1b_d)
    {
        un(idx) = uc(idx) + (*rx_d) * ( uc(idx-1) + 2.*uc(idx) + uc(idx+1) ) + 
                  (*dt_d) * FORCE( x(idx), *t_d );
    }

    #undef uc
    #undef un
    #undef x
}


// ------------------------------------------------------------------------
// Return the current wall-clock time in seconds
// ------------------------------------------------------------------------
inline double getCPU(){
    return ( 1.0 * std::clock() )/CLOCKS_PER_SEC ;
}

#include "parseCommand.h"

// ------------------------------------------------------------------------
// Function to save a vector to a matlab file
// matlabFile  (input) : save vector to this file
// u_p         (input) : array of vector values
// name        (input) : name for array
// (nd1a:nd1b) (input) : array dimensions
// ------------------------------------------------------------------------
int writeMatlabVector(FILE *matlabFile, Real *u_p, const char *name, int nd1a, int nd1b){
    #define u(i) u_p[i-nd1a]

    const int numPerLine=8; // number of entries per line
    // Save the vector as:
    // name = [ num num num num num ...
    // num num num num num ];
    fprintf(matlabFile,"%s=[",name);
    for( int i=nd1a; i<=nd1b; i++ ) {
        fprintf(matlabFile,"%20.15e ",u(i));
        if( (i-nd1a) % numPerLine == numPerLine-1 )
            fprintf(matlabFile,"...\n"); // continuation line
    }
    fprintf(matlabFile,"];\n");

    return 0;
}


int main(int argc, char *argv[]){

    // parameters and names used if it is not provided 
    int Nx = 20;           // set spatial discretization size
    std::string matlabFileName = "heat1d_gpu.m";
    int debug= 0;
    int saveMatlab=0;
    Real tFinal = 0.2;
    Real cfl    = 0.9;

    std::string line;
    bool echo = true;
    for (int i=1; i<argc; i++)
    {
        line = argv[i];
        if( parseCommand(line, "-nx=", Nx, echo) ) {}
        else if( parseCommand(line, "-debug=", debug, echo) ) {}
        else if( parseCommand(line, "-tFinal=", tFinal, echo) ) {}
        else if( parseCommand(line, "-saveMatlab=", saveMatlab, echo) ) {}
        else if( parseCommand(line, "-matlabFileName=", matlabFileName, echo) ) {}
    }
    
    Real xa = 0.0, xb = 1.0 ;
    Real *xa_d;
    checkCudaErrors( cudaMalloc((void **)&xa_d, sizeof(Real)) ); cudaMemcpy(xa_d, &xa, sizeof(Real), cudaMemcpyHostToDevice);

    // ============= Grid and indexing==============
    //            xa                             xb
    //         G---X---+---+---+---+-- ... ---+---X---G
    //             0   1   2                      Nx
    //             ja                             jb
    //        nd1a                                   nd1b
    // C index: 0 1 2 3 ...

    Real dx = (xb - xa)/Nx;
    Real *dx_d;
    checkCudaErrors( cudaMalloc((void **)&dx_d, sizeof(Real)) ); cudaMemcpy(dx_d, &dx, sizeof(Real), cudaMemcpyHostToDevice);
    int numGhost = 1;
    int n1a       = 0;
    int n1b       = Nx;
    int nd1a     = n1a - numGhost;
    int nd1b     = n1b + numGhost;
    int nd1      = nd1b - nd1a + 1; // total number of grid points including ghost nodes

    int nt = MAX_NUM_THREADS;
    int nb = ceil( (1.*nd1)/nt );
    int *nt_d;
    checkCudaErrors( cudaMalloc((void **)&nt_d, sizeof(int)) );  cudaMemcpy(nt_d, &nt, sizeof(int), cudaMemcpyHostToDevice);

    int *n1a_d, *n1b_d, *nd1a_d, *nd1b_d, *nd1_d;
    checkCudaErrors( cudaMalloc((void **)&n1a_d, sizeof(int)) );  cudaMemcpy(n1a_d, &n1a, sizeof(int), cudaMemcpyHostToDevice);
    checkCudaErrors( cudaMalloc((void **)&n1b_d, sizeof(int)) );  cudaMemcpy(n1b_d, &n1b, sizeof(int), cudaMemcpyHostToDevice);
    checkCudaErrors( cudaMalloc((void **)&nd1a_d, sizeof(int)) ); cudaMemcpy(nd1a_d,&nd1a,sizeof(int), cudaMemcpyHostToDevice);
    checkCudaErrors( cudaMalloc((void **)&nd1b_d, sizeof(int)) ); cudaMemcpy(nd1b_d,&nd1b,sizeof(int), cudaMemcpyHostToDevice);
    checkCudaErrors( cudaMalloc((void **)&nd1_d, sizeof(int)) );  cudaMemcpy(nd1_d, &nd1, sizeof(int), cudaMemcpyHostToDevice);
    // creating a 1D array of grid points
    Real *x_p = new Real [nd1];
    # define x(i) x_p[i-nd1a] 
    Real *x_d;
    cudaCheckErrors( cudaMalloc((void **)&x_d, nd1*sizeof(Real)) ); 
    mesh1d<<<nb, nt>>>(x_d, nd1_d, nd1a_d, nd1b_d, xa_d, dx_d, nt_d);
    cudaMemcpy(x_p, x_d, nd1*sizeof(Real), cudaMemcpyDeviceToHost);

    const int numSides = 2;
    const int dirichlet = 1, neumann = 2;
    const int numberOfDimensions = 1;
    // initialize boundary conditions as a 2D matrix with rows for left/right side and columns for dimensions
    int *boundaryCondition_p = new int [numSides * numberOfDimensions];
    #define boundaryCondition(side, axis) boundaryCondition_p[(side) + numSides*(axis)]

    // assign boundary conditions for different problems
    #if SOLUTION == TRIG_DD
        // true solution for dirichlet BC's
        boundaryCondition(0, 0) = dirichlet;
        boundaryCondition(1, 0) = dirichlet;
    
    #elif SOLUTION == TRIG_NN
        // true solution for neumann BC's
        boundaryCondition(0, 0) = neumann;
        boundaryCondition(1, 0) = neumann;

    #elif SOLUTION == POLY_DD || SOLUTION == POLY_NN
        #if SOLUTION == POLY_DD
            boundaryCondition(0, 0) = dirichlet;
            boundaryCondition(1, 0) = dirichlet;

        #else
            boundaryCondition(0, 0) = neumann;
            boundaryCondition(1, 0) = neumann;
        
        #endif
    #else
        std::cerr << "ERROR: unknown choice of solution/case to solve \n";
        abort();
    #endif

    Real *u_p[2]; // two arrays used for storing current and next time step solution vectors
    u_p[0] = new Real [nd1];
    u_p[1] = new Real [nd1];
    
    #define uc(i) u_p[curr][i-nd1a]
    #define un(i) u_p[next][i-nd1a]

    // initial condition set up
    Real t = 0.0;
    Real *t_d;
    checkCudaErrors( cudaMalloc((void **)&t_d, sizeof(Real)) ); cudaMemcpy(t_d, &t, sizeof(Real), cudaMemcpyHostToDevice);
    int curr = 0;
    Real *u_d[2];
    checkCudaErrors( cudaMalloc(void **)&u_d[0], nd1*sizeof(Real) );
    checkCudaErrors( cudaMalloc(void **)&u_d[1], nd1*sizeof(Real) );
    Real *uc_h = &uc(nd1a);
    Real *un_h = nullptr;
    Real *uc_d = u_d+curr;
    Real *un_d = nullptr;
    setInitialCondition(uc_d, nd1_d, nd1a_d, nd1b_d, x_d, t_d, nt_d);
    cudaMemcpy(uc_h, uc_d, nd1*sizeof(Real), cudaMemcpyDeviceToHost);

    /* Time-step restrictions */
    const Real dx2     = dx * dx;
    Real dt            = cfl * 0.5 * dx2 / kappa;
    const int numSteps = ceil( tFinal/dt );
    dt                 = tFinal/numSteps;
    Real *dt_d;
    checkCudaErrors( cudaMalloc((void **)&dt_d, sizeof(Real)) ); cudaMemcpy(dt_d, &dt, sizeof(Real), cudaMemcpyHostToDevice);
    const Real rx      = kappa * dt / dx2;
    Real *rx_d;
    checkCudaErrors( cudaMalloc((void **)&rx_d, sizeof(Real)) ); cudaMemcpy(rx_d, &rx, sizeof(Real), cudaMemcpyHostToDevice);

    printf("------------------- Solve the heat equation in 1D solution=%s --------------------- \n",solutionName);
    printf("  numGhost=%d, n1a=%d, n1b=%d, nd1a=%d, nd1b=%d\n",numGhost,ja,jb,nd1a,nd1b);
    printf("  numSteps=%d, Nx=%d, kappa=%g, tFinal=%g, boundaryCondition(0,0)=%d, boundaryCondition(1,0)=%d\n",numSteps,Nx,kappa,tFinal,boundaryCondition(0,0),boundaryCondition(1,0));
    printf(" ----------------- Using Forward Euler Time Stepping -------------------------- \n");

    /* -------- TIME-STEPPING LOOP --------- */
    Real cpu0 = getCPU();
    for (int n=0; n<numSteps; n++) {
        const int curr = n % 2;
        const int next = (n + 1) % 2;
        t = n * dt; // current time
        cudaMemcpy(t_d, &t, sizeof(Real), cudaMemcpyHostToDevice);
        
        uc_h = &uc(nd1a);
        un_h = &un(nd1a);

        uc_d = u_d+curr;
        un_d = u_d+next;
        cudaMemcpy(uc_d, uc_h, nd1*sizeof(Real), cudaMemcpyHostToDevice);
        heat1dForwardEulerTimeStep<<<nb, nt>>>(uc_d, un_d, x_d, nd1_d, nd1a_d, nd1b_d, rx_d, dt_d, t_d, nt_d);
        cudaMemcpy(un_h, un_d, nd1*sizeof(Real), cudaMemcpyDeviceToHost);
        // set boundary conditions on host
        for (int side=0; side<=1; side++) {
            const int i  = side == 0 ? ja : jb;
            const int is = 1 - 2*side;
            if (boundaryCondition(side, 0) == dirichlet){
                un(i)    = UTRUE( x(i), t+dt );
                un(i-is) = 3.0 * un(i) - 3.0 * un(i+is) + un(i+2*is); // extrapolate ghost
            }
            else {
                // neumann
                un(i-is) = un(i+is) - 2.*is*dx*UTRUEX(x(i),t+dt);
            }
        }

    }
    Real cpuTimeStep = getCPU()-cpu0;

    /*  check errors */
    t += dt; // tFinal
    if ( fabs(t - tFinal) > 1e-3 * dt/tFinal ){
        printf("ERROR AFTER TIME STEPPING: t=%16.8e IS NOT EQUAL to tFinal=%16.8e\n", t, tFinal);
    }

    Real *error_p = new Real [nd1];
    #define error(i) error_p[i-nd1a]

    curr = numSteps % 2;
    Real maxErr = 0.0;
    for (int i=nd1a; i<=nd1b; i++) {
        error(i) = uc(i) - UTRUE( x(i), t );
        maxErr = max( maxErr, abs(error(i)) );
    }
    printf("numSteps=%4d, Nx=%3d, maxErr=%9.2e, cpu=%9.2e(s)\n",numSteps,Nx,maxErr,cpuTimeStep);

// --- Write a file for plotting in matlab ---
    FILE *matlabFile = fopen(matlabFileName.c_str(),"w");
    fprintf(matlabFile,"%% File written by heat1d.C\n");
    fprintf(matlabFile,"xa=%g; xb=%g; kappa=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",xa,xb,kappa,tFinal,maxErr,cpuTimeStep);
    fprintf(matlabFile,"Nx=%d; dx=%14.6e; numGhost=%d; n1a=%d; n1b=%d; nd1a=%d; nd1b=%d;\n",Nx,dx,numGhost,ja,jb,nd1a,nd1b);
    fprintf(matlabFile,"solutionName=\'%s\';\n",solutionName);

    if (saveMatlab > 1)
    {
        writeMatlabVector( matlabFile, x_p, "x", nd1a, nd1b );
        writeMatlabVector( matlabFile, u_p[curr], "u", nd1a, nd1b );
        writeMatlabVector( matlabFile, error_p, "err", nd1a, nd1b );
    }

    fclose(matlabFile);
    printf("Wrote file %s\n\n",matlabFileName.c_str());

    delete [] x_p;
    delete [] boundaryCondition_p;
    delete [] u_p[0];
    delete [] u_p[1];
    delete [] error_p;

    cudaFree(xa_d); 
    cudaFree(nt_d); 
    cudaFree(n1a_d); cudaFree(n1b_d); cudaFree(nd1a_d); cudaFree(nd1b_d); cudaFree(nd1_d);
    cudaFree(x_d); 
    cudaFree(u_d[0]); cudaFree(u_d[1]);
    cudaFree(t_d); cudaFree(dt_d); cudaFree(rx_d); 

    return 0;
}