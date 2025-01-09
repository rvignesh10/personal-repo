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

#include <string>
using std::string;
using std::max;

#include <ctime>
// ------------------------------------------------------------------------
// Return the current wall-clock time in seconds
// ------------------------------------------------------------------------
inline double getCPU(){
    return ( 1.0 * std::clock() )/CLOCKS_PER_SEC ;
}

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
    // defining max arguments
    #define MAX_ARG 5

    printf("Usage: ./heat1d -Nx [Nx] -matlabFileName [matlabFileName.m] \n");
    printf("       Nx = number of grid cells. \n");
    printf("       matlabFileName.m : save result to this file \n");

    // defining different type of solutions available in the code
    #define TRIG_DD 1
    #define TRIG_NN 2
    #define POLY_DD 3
    #define POLY_NN 4

    // parameters and names used if it is not provided 
    int Nx = 20;           // set spatial discretization size
    std::string matlabFileName = "default_out.m";

    #ifndef SOLUTION
        // #define SOLUTION TRIG_DD
        #define SOLUTION TRIG_NN
        // #define SOLUTION POLY_DD
        // #define SOLUTION POLY_NN
    #endif


    // parsing input parameters
    if (argc > 1){
        for (int i=1; i<argc; ) {
            string option = argv[i];
            if (i+1 == argc){
                std::cerr << "Not enough arguments provided \n";
                exit(EXIT_FAILURE);
            }
            else if (option == "-Nx"){
                Nx = atoi(argv[i+1]);
            }
            else if (option == "-matlabFileName"){
                matlabFileName = argv[i+1];
            }
            i += 2;
        }
    }
    
    const Real pi = M_PI;

    Real xa = 0.0, xb = 1.0 ;
    Real kappa  = 0.1;
    Real tFinal = 0.2;
    Real cfl    = 0.9;

    // ============= Grid and indexing==============
    //            xa                             xb
    //         G---X---+---+---+---+-- ... ---+---X---G
    //             0   1   2                      Nx
    //             ja                             jb
    //        nd1a                                   nd1b
    // C index: 0 1 2 3 ...

    Real dx = (xb - xa)/Nx;
    const int numGhost = 1;
    const int ja       = 0;
    const int jb       = Nx;
    const int nd1a     = ja - numGhost;
    const int nd1b     = jb + numGhost;
    const int nd1      = nd1b - nd1a + 1; // total number of grid points including ghost nodes

    // creating a 1D array of grid points
    Real *x_p = new Real [nd1];
    # define x(i) x_p[i-nd1a] 

    for (int i=nd1a; i<=nd1b; i++) {
        x(i) = xa + (i - ja) * dx;
    }

    const int numSides = 2;
    const int dirichlet = 1, neumann = 2;
    const int numberOfDimensions = 1;
    // initialize boundary conditions as a 2D matrix with rows for left/right side and columns for dimensions
    int *boundaryCondition_p = new int [numSides * numberOfDimensions];
    #define boundaryCondition(side, axis) boundaryCondition_p[(side) + numSides*(axis)]

    // some necessary constants
    const Real kx = 3.0;
    const Real kxpi = kx * pi;
    const Real kappaPiSq = kappa * kxpi * kxpi;

    const Real b0 = 1.0, b1 = 0.5, b2 = 0.25;
    const Real a0 = 1.0, a1 = 0.3;

    // assign boundary conditions for different problems
    #if SOLUTION == TRIG_DD
        // true solution for dirichlet BC's
        boundaryCondition(0, 0) = dirichlet;
        boundaryCondition(1, 0) = dirichlet;

        const char solutionName[] = "trueDD";

        #define UTRUE(x, t) sin( kxpi * (x) ) * exp( -kappaPiSq * (t) )
        #define UTRUEX(x, t) kxpi * cos( kxpi * (x) ) * exp( -kappaPiSq * (t) )
        #define UTRUET(x, t) -kappaPiSq * UTRUE(x, t)
        #define FORCE(x, t) (0.)
    
    #elif SOLUTION == TRIG_NN
        // true solution for neumann BC's
        boundaryCondition(0, 0) = neumann;
        boundaryCondition(1, 0) = neumann;

        const char solutionName[] = "trueNN";

        #define UTRUE(x, t) cos( kxpi * (x) ) * exp( -kappaPiSq * (t) )
        #define UTRUEX(x, t)  -kxpi * sin( kxpi * (x) ) * exp( -kappaPiSq * (t) )
        #define UTRUET(x, t) -kappaPiSq * UTRUE(x, t)
        #define FORCE(x, t) (0.)

    #elif SOLUTION == POLY_DD || SOLUTION == POLY_NN
        #if SOLUTION == POLY_DD
            boundaryCondition(0, 0) = dirichlet;
            boundaryCondition(1, 0) = dirichlet;

            const char solutionName[] = "polyDD";
        
        #else
            boundaryCondition(0, 0) = neumann;
            boundaryCondition(1, 0) = neumann;

            const char solutionName[] = "polyNN";
        
        #endif
    
        #define UTRUE(x, t) ( b0 + (x) * (b1 + (x) * b2) ) * ( a0 + (t) * a1 )
        #define UTRUEX(x, t) ( b1 + 2. * (x) * b2 ) * ( a0 + (t) * ( a1 ) )
        #define UTRUET(x, t) (b0 + (x) * ( b1 + (x) * b2 )) * ( a1 )
        #define UTRUEXX(x, t) ( 2.*b2 )*( a0 + (t)*( a1 ) )
        #define FORCE(x, t) ( UTRUET(x, t) - kappa*UTRUEXX(x, t) )
    
    #else
        const char solutionName[] = "wrongChoice";
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
    int curr = 0;
    for (int i=nd1a; i<=nd1b; i++) {
        uc(i) = UTRUE(x(i), t);
    }

    /* Time-step restrictions */
    const Real dx2     = dx * dx;
    Real dt            = cfl * 0.5 * dx2 / kappa;
    const int numSteps = ceil( tFinal/dt );
    dt                 = tFinal/numSteps;
    const Real rx      = kappa * dt / dx2;
    
    printf("------------------- Solve the heat equation in 1D solution=%s --------------------- \n",solutionName);
    printf("  numGhost=%d, n1a=%d, n1b=%d, nd1a=%d, nd1b=%d\n",numGhost,ja,jb,nd1a,nd1b);
    printf("  numSteps=%d, Nx=%d, kappa=%g, tFinal=%g, boundaryCondition(0,0)=%d, boundaryCondition(1,0)=%d\n",numSteps,Nx,kappa,tFinal,boundaryCondition(0,0),boundaryCondition(1,0));

    #define FORWARD_EULER  1
    #define RK_2 2
    // make a choice for explicit time marching scheme
    #ifndef EXPLICIT_SOLVER
        //#define EXPLICIT_SOLVER FORWARD_EULER
        #define EXPLICIT_SOLVER RK_2
    #endif

    #if EXPLICIT_SOLVER == FORWARD_EULER
        printf(" ----------------- Using Forward Euler Time Stepping -------------------------- \n");
    #else
         printf(" ----------------- Using Explicit RK2 Time Stepping -------------------------- \n");
    #endif

    /* -------- TIME-STEPPING LOOP --------- */
    Real cpu0 = getCPU();
    for (int n=0; n<numSteps; n++) {
        const int curr = n % 2;
        const int next = (n + 1) % 2;
        t = n * dt; // current time

        #if EXPLICIT_SOLVER == FORWARD_EULER
            for (int j=ja; j<=jb; j++) {
                un(j) = uc(j) + rx * ( uc(j+1) - 2*uc(j) + uc(j-1) ) + dt * FORCE( x(j), t );
            }
        #elif EXPLICIT_SOLVER == RK_2
            // write RK2 solver
            Real *w_p = new Real[nd1];
            #define w(i) w_p[i-nd1a] 

            for (int i=ja; i<=jb; i++) {
                w(i) = uc(i) + (rx/2.0) * (uc(i+1) - 2.*uc(i) + uc(i-1)) + 0.5 * dt * FORCE( x(i), t );
            }
            // set boundary conditions
            for (int side=0; side<=1; side++) {
                const int i  = side == 0 ? ja : jb;
                const int is = 1 - 2*side;
                if (boundaryCondition(side, 0) == dirichlet){
                    w(i)    = UTRUE( x(i), t+ 0.5*dt );
                    w(i-is) = 3.0 * un(i) - 3.0 * un(i+is) + un(i+2*is); // extrapolate ghost
                }
                else {
                    // neumann
                    w(i-is) = w(i+is) - 2.*is*dx*UTRUEX(x(i),t+ 0.5*dt);
                }
            }

            for (int i=ja; i<=jb; i++) {
                un(i) = uc(i) + rx * (w(i+1) - 2.*w(i) + w(i-1)) + dt * FORCE( x(i), t+0.5*dt );
            }            

            delete [] w_p;
        #else
            std::cerr << "ERROR: wrong time integrator chosen \n";
            abort();
        #endif
        // set boundary conditions
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

    writeMatlabVector( matlabFile, x_p, "x", nd1a, nd1b );
    writeMatlabVector( matlabFile, u_p[curr], "u", nd1a, nd1b );
    writeMatlabVector( matlabFile, error_p, "err", nd1a, nd1b );

    fclose(matlabFile);
    printf("Wrote file %s\n\n",matlabFileName.c_str());

    delete [] x_p;
    delete [] boundaryCondition_p;
    delete [] u_p[0];
    delete [] u_p[1];
    delete [] error_p;
    return 0;
}
