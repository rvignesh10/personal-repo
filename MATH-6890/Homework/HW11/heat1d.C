#include <stdio.h>
#include <iostream>
#include <math.h>
#include <float.h>
#include <assert.h>

#define USE_PPP
#define NUM_THREADS 128

#include <string>
using std::string;
using std::max;

#include <ctime>
// ------------------------------------------------------------------------
// Return the current wall-clock time in seconds
// ------------------------------------------------------------------------
inline double getCPU()
{
  #ifndef USE_PPP
    return MPI_Wtime();
  #else
    return (1.0*std::clock())/CLOCKS_PER_SEC;
  #endif
}

// MPI header file
#include <mpi.h>

// ------------------------------------------------------------------------
// wrapper function definitions
// ------------------------------------------------------------------------
void heat1dSetInitialCondition(double *dev_uc, double *dev_x, int Nb, int Nt, int *dev_Nt, int *dev_nd1L);
void heat1dForwardEulerUpdate(double *dev_uc, double *dev_un, double *dev_x, double *dev_rx, double *dev_t, 
                                         int Nb, int Nt, int *dev_Nt, int *dev_n1aL, int *dev_n1bL, int *dev_numGhost, float *gpuStepTime);
void setSpatialMeshNodes1D( double *dev_x, double *dev_xa, double *dev_dx, int *dev_nd1L, int Nb, int Nt, int *dev_Nt );
void heat1dErrorCalc(double *dev_err, double *dev_uc, double *dev_x, double *dev_t, int Nb, int Nt, int *dev_Nt, int *dev_nd1L);
void AllocateCudaMemory( double *dev_var, int n_size );
void AllocateCudaMemory( int *dev_var, int n_size );
void FreeCudaMemory( double *dev_var );
void FreeCudaMemory( int *dev_var );
void MemcpyHostToDev( double *dev_var, double *hst_var, int n_size );
void MemcpyHostToDev( int *dev_var, int *hst_var, int n_size );
void MemcpyDevToHost( double *hst_var, double *dev_var, int n_size );
void MemcpyDevToHost( int *hst_var, int *dev_var, int n_size );

// parseCommand
#include "parseCommand.h"
// getting LocalIndexBounds
#include "getLocalIndexBounds.h"

#define POLY_DD 0
#define POLY_NN 1

typedef double Real;
typedef int Integer;

enum bcs_send_comm {
    bcs_Lsend = 0,
    bcs_Rsend = 1
};

enum bcs_recv_comm {
    bcs_Lrecv = 1,
    bcs_Rrecv = 0
};

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

int main( int argc, char *argv[] ) {
    MPI_Init(&argc, &argv);
    Integer myRank, np;
    Integer root = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Request sendLeft, recvLeft;
    MPI_Request sendRight, recvRight;
    Integer mpi_op;
    if (myRank == 0)
    {
        printf("Usage: heat1d -nx=<i> -option=[0|1] -tFinal=<f> -root=<i> -debug=<i> -saveMatlab=[0|1|2] - matlabFile=<s> -commOption=<i>\n""       option : 0=polyDD, 1=polyNN\n");
    }

    Real tFinal = .5;
    Integer nx  = 100;
    Integer saveMatlab=2;
    Integer commOption=0;
    Integer option= 0;
    Integer debug = 1;
    string matlabFileName= "heat1dMPI.m";

    string line;
    bool echo = root==myRank? true : false;
    for (int i=1; i<argc; i++)
    {
        line = argv[i];
        if( parseCommand(line, "-nx=", nx, echo) ) {}
        else if( parseCommand(line, "-option=", option, echo) ) {}
        else if( parseCommand(line, "-root=", root, echo) ) {}
        else if( parseCommand(line, "-debug=", debug, echo) ) {}
        else if( parseCommand(line, "-tFinal=", tFinal, echo) ) {}
        else if( parseCommand(line, "-saveMatlab=", saveMatlab, echo) ) {}
        else if( parseCommand(line, "-matlabFileName=", matlabFileName, echo) ) {}
    }

    FILE *debugFile = NULL;
    if (debug > 0)
    {
        char debugFileName[80];
        sprintf(debugFileName, "debug/heat1dNp%dProc%d.debug", np, myRank);
        debugFile=fopen(debugFileName, "w");
    }

    Real xa = 0.0, xb = 1.0;
    Real kappa= 0.1;
    Real cfl  = 0.9;

    Real dx = (xb - xa)/nx;
    // allocate and copy memory from host to device
    Real *dev_xa, *dev_dx;
    AllocateCudaMemory(dev_xa, 1); MemcpyHostToDev(dev_xa, &xa, 1);
    AllocateCudaMemory(dev_dx, 1); MemcpyHostToDev(dev_dx, &dx, 1);

    Integer numGhost = 1;
    int *dev_numGhost;
    AllocateCudaMemory(dev_numGhost, 1); MemcpyHostToDev(dev_numGhost, &numGhost, 1);
    Integer nx_l, n1a_l, n1b_l;
    getLocalIndexBounds(myRank, np, nx, nx_l, n1a_l, n1b_l);
    Integer nd1a_l = n1a_l - numGhost;
    Integer nd1b_l = n1b_l + numGhost;
    Integer nd1_l = nd1b_l - nd1a_l + 1;

    // allocate and copy memory from host to device
    int *dev_n1aL, *dev_n1bL, *dev_nd1L;
    AllocateCudaMemory(dev_n1aL, 1); MemcpyHostToDev(dev_n1aL, &n1a_l, 1);
    AllocateCudaMemory(dev_n1bL, 1); MemcpyHostToDev(dev_n1bL, &n1b_l, 1);
    AllocateCudaMemory(dev_nd1L, 1); MemcpyHostToDev(dev_nd1L, &nd1_l, 1);

    // Local spatial grid
    Real *x_l = new Real [nd1_l];
    #define x(id) x_l[(id) - nd1a_l] 
    
    // find cuda blocks and threads
    int Nt = NUM_THREADS;
    int Nb = ceil( (1.*(nd1_l))/Nt );
    int *dev_Nt;
    AllocateCudaMemory(dev_Nt, 1); MemcpyHostToDev(dev_Nt, &Nt, 1);
    // allocate and copy memory from host to device for spatial location
    Real *dev_x;
    AllocateCudaMemory(dev_x, nd1_l);
    setSpatialMeshNodes1D(dev_x, dev_xa, dev_dx, dev_nd1L, Nb, Nt, dev_Nt);
    MemcpyDevToHost(x_l, dev_x, nd1_l);

    if (debug>0)
    {
        fprintf(debugFile, "Local discretization: \n\tnx: %d\n\tn1a: %d\n\tn1b: %d\n\tnd1: %d\n", 
        nx_l, n1a_l, n1b_l, nd1_l);
        fprintf(debugFile, " X -- X -- X -- ....... -- X -- X\n%d    %d         .......    %d    %d \n", nd1a_l, n1a_l, n1b_l, nd1b_l);
        for (int i=nd1a_l; i<=nd1b_l; i++)
            fprintf(debugFile, "%1.1f\t", x(i));
        fprintf(debugFile, "\n");
    }

    const Integer numSides  = 2;
    const Integer dirichlet = 1, neumann = 2;
    const Integer numberOfDimensions = 1;
    Integer *boundaryCondition_l = new Integer [numSides * numberOfDimensions];
    #define boundaryCondition(side, axis) boundaryCondition_l[(side) + numSides*(axis)]

    string solutionName;
    if (option == POLY_DD)
    {
        boundaryCondition(0, 0) = dirichlet;
        boundaryCondition(1, 0) = dirichlet;
        solutionName = "polyDD";
    }
    else if( option == POLY_NN )
    {
        boundaryCondition(0, 0) = neumann;
        boundaryCondition(1, 0) = neumann;
        solutionName= "polyNN";
    }
    else{
        solutionName= "ERRORS";
        if(myRank == 0)
            fprintf(stdout, "ERROR: wrong option chosen... aborting \n");
        abort();
    }
    const Real b0 = 1.0, b1 = 0.5, b2 = 0.25;
    const Real a0 = 1.0, a1 = 0.3;

    #define UTRUE(x, t) ( b0 + (x) * (b1 + (x) * b2) ) * ( a0 + (t) * a1 )
    #define UTRUEX(x, t) ( b1 + 2. * (x) * b2 ) * ( a0 + (t) * ( a1 ) )
    #define UTRUET(x, t) (b0 + (x) * ( b1 + (x) * b2 )) * ( a1 )
    #define UTRUEXX(x, t) ( 2.*b2 )*( a0 + (t)*( a1 ) )
    #define FORCE(x, t) ( UTRUET(x, t) - kappa*UTRUEXX(x, t) )

    Real *u_l[2]; // two local arrays used for storing current and next time-step solution vector
    u_l[0] = new Real [nd1_l];
    u_l[1] = new Real [nd1_l];

    #define uc(id) u_l[curr][(id)-nd1a_l]
    #define un(id) u_l[next][(id)-nd1a_l]

    double *dev_uc, *dev_un;
    AllocateCudaMemory(dev_uc, nd1_l);
    AllocateCudaMemory(dev_un, nd1_l);

    // initial conditions set-up
    Real t = 0.0;
    // allocate memory and copy value on gpu
    double *dev_t;
    AllocateCudaMemory(dev_t, 1); MemcpyHostToDev(dev_t, &t, 1);
    Integer curr = 0;
    heat1dSetInitialCondition(dev_uc, dev_x, Nb, Nt, dev_Nt, dev_nd1L);
    MemcpyDevToHost(&uc(nd1a_l), dev_uc, nd1_l);

    if(debug>0)
    {   
        fprintf(debugFile, "u_l(x, 0): [\t");
        for(int i=nd1a_l; i<=nd1b_l; i++)
        fprintf(debugFile, "%5.9e\t", uc(i));
        fprintf(debugFile, "\t]\n");
        fflush(debugFile);
    }

    /* Time-step restrictions */
    const Real dx2     = dx * dx;
    Real dt            = cfl * 0.5 * dx2 / kappa;
    const int numSteps = ceil( tFinal/dt );
    dt                 = tFinal/numSteps;
    Real rx      = kappa * dt / dx2;
    double *dev_rx;
    AllocateCudaMemory(dev_rx, 1); MemcpyHostToDev(dev_rx, &rx, 1);

    if(myRank == 0)
    {
        printf("------------------- Solve the heat equation in 1D in parallel with solution=%s --------------------- \n", solutionName.c_str());        
    }

    /* ---- Time stepping loop ---- */
    Real cpu0 = getCPU();
    float gpuTime = 0.0;
    for (int n=0; n<numSteps; n++)
    {
        const Integer curr = n%2;
        const Integer next = (n+1)%2;
        t = n*dt; // current time

        float gpuStepTime;
        // copy current value to device from host
        MemcpyHostToDev(dev_t, &t, 1);
        MemcpyHostToDev(dev_uc, &uc(nd1a_l), nd1_l);
        // perform time stepping in GPU
        heat1dForwardEulerUpdate(dev_uc, dev_un, dev_x, dev_rx, dev_t, Nb, Nt, dev_Nt, dev_n1aL, dev_n1bL, dev_numGhost, &gpuStepTime);
        gpuTime += gpuStepTime;
        // copy values of un back to host
        MemcpyDevToHost(&un(nd1a_l), dev_un, nd1_l);

        // set boundary conditions
        if (np > 1)
        {
            if(myRank == 0)
            {   
                // right side
                MPI_Send(&un(n1b_l), 1, MPI_DOUBLE, myRank+1, bcs_Rsend, MPI_COMM_WORLD);
                // left side
                if(boundaryCondition(0, 0) == dirichlet)
                {
                    un(n1a_l) = UTRUE( x(n1a_l), t+dt );
                    un(nd1a_l)= 3.0*un(n1a_l) - 3.0*un(n1a_l+1) + un(n1a_l+2);
                }
                else
                {
                    // neumann
                    un(n1a_l) = UTRUE( x(n1a_l), t+dt );
                    un(nd1a_l)= un(n1a_l+1) - 2.*dx*UTRUEX( x(n1a_l), t+dt );
                }
                MPI_Recv(&un(nd1b_l), 1, MPI_DOUBLE, myRank+1, bcs_Rrecv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else if(myRank == np-1 )
            {   
                // left side
                MPI_Send(&un(n1a_l), 1, MPI_DOUBLE, myRank-1, bcs_Lsend, MPI_COMM_WORLD);
                // right side
                if( boundaryCondition(1, 0) == dirichlet )
                {
                    un(n1b_l) = UTRUE( x(n1b_l), t+dt );
                    un(nd1b_l)= 3.0*un(n1b_l) - 3.0*un(n1b_l-1) + un(n1b_l-2);
                }
                else
                {   
                    // neumann
                    un(n1b_l) = UTRUE( x(n1b_l), t+dt );
                    un(nd1b_l)= un(n1b_l-1) + 2.*dx*UTRUEX( x(n1b_l), t+dt );
                }
                MPI_Recv(&un(nd1a_l), 1, MPI_DOUBLE, myRank-1, bcs_Lrecv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {   
                // left and right side
                MPI_Send(&un(n1a_l), 1, MPI_DOUBLE, myRank-1, bcs_Lsend, MPI_COMM_WORLD);
                MPI_Send(&un(n1b_l), 1, MPI_DOUBLE, myRank+1, bcs_Rsend, MPI_COMM_WORLD);

                MPI_Recv(&un(nd1a_l), 1, MPI_DOUBLE, myRank-1, bcs_Lrecv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&un(nd1b_l), 1, MPI_DOUBLE, myRank+1, bcs_Rrecv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
            }
        }
        else
        {
            for (int side=0; side<=1; side++) {
                const int i  = side == 0 ? n1a_l : n1b_l;
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

        if(debug>1)
        {
            // check and write max errors
            Real maxErr = 0.;
            for( int i=nd1a_l; i<=nd1b_l; i++ )
            {
                Real err = fabs( un(i) - UTRUE( x(i), t+dt ) );
                maxErr = max( maxErr, err );
            }
            fprintf(debugFile, "step= %d, t= %9.3e, maxErr= %9.2e\n", n+1, t+dt, maxErr);
            fflush(debugFile);
        }
    }
    Real cpuTimeStep = getCPU()-cpu0;
    Real cpuTimeStep_g = cpuTimeStep;
    Real gpuTime_g = gpuTime;
    MPI_Reduce(&gpuTime, &gpuTime_g, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    mpi_op = MPI_Reduce(&cpuTimeStep, &cpuTimeStep_g, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    if (mpi_op != MPI_SUCCESS)
        abort();
    MPI_Barrier(MPI_COMM_WORLD);
    // ---------- check final error -------------
    t += dt;
    if( fabs(t-tFinal) > 1e-03*dt/tFinal )
    {
        if(myRank == 0)
        {
            fprintf(stdout, "ERROR: AFTER TIME STEPPING: t=%16.8e IS NOT EQUAL TO tFinal=%16.8e\n", t, tFinal);
        }
    }

    curr = numSteps%2;
    Integer NTot = 0;
    Integer sendCount;
    if(np > 1)
        sendCount = (myRank==0 || myRank==np-1) ? nx_l+2 : nx_l+1;
    else
        sendCount = nx_l+1+2*numGhost;
    Integer *recvCount_p = new Integer [np];
    MPI_Allgather(&sendCount, 1, MPI_INTEGER, recvCount_p, 1, MPI_INTEGER, MPI_COMM_WORLD);
    MPI_Allreduce(&sendCount, &NTot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    Integer *displacement_p = new Integer [np];
    Integer sum=0;
    for(int i=0; i<np; i++)
    {
        displacement_p[i] = sum;
        sum += recvCount_p[i];
    }

    if(debug > 0)
    {   
        fprintf(debugFile, "uLocal: [ \t");
        for(int i=nd1a_l; i<=nd1b_l; i++)
        fprintf(debugFile, "%20.15e\t", uc(i));
        fprintf(debugFile, "]\n");
        fprintf(debugFile, "recvCount: ");
        for(int i=0; i<np; i++)
        {
            fprintf(debugFile, "%d \t", recvCount_p[i]);
        }
        fprintf(debugFile, "\n");
        fprintf(debugFile, "displacement: ");
        for(int i=0; i<np; i++)
        {
            fprintf(debugFile, "%d \t", displacement_p[i]);
        }
        fprintf(debugFile, "\n");
        fprintf(debugFile, "NTot: %d\n", NTot);
        fflush(debugFile);
    }



    Real *err_p = new Real [nd1_l];
    #define err(i) err_p[(i) - nd1a_l]

    Real *dev_err;
    AllocateCudaMemory(dev_err, nd1_l);
    MemcpyHostToDev(dev_uc, &uc(nd1a_l), nd1_l);
    MemcpyHostToDev(dev_t, &t, 1);
    heat1dErrorCalc(dev_err, dev_uc, dev_x, dev_t, Nb, Nt, dev_Nt, dev_nd1L);
    MemcpyDevToHost(err_p, dev_err, nd1_l);

    Real maxErr_l = 0., maxErr_g;
    for(int i=nd1a_l; i<=nd1b_l; i++)
    {
        maxErr_l= max(maxErr_l, abs(err(i)));
    }
    maxErr_g = maxErr_l;
    MPI_Reduce(&maxErr_l, &maxErr_g, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    
    if (debug>0)
    {
        fprintf(debugFile, "maxErr local : %1.9e\n", maxErr_l);
        fprintf(debugFile, "maxErr global: %1.9e\n", maxErr_g);
        fflush(debugFile);
    }

    Real *sendArray_p = new Real [sendCount];
    Real *sendError_p = new Real [sendCount];
    Real *sendxloc_p  = new Real [sendCount];
    Integer start = myRank==0   ? nd1a_l : n1a_l;
    Integer end   = myRank==np-1? nd1b_l : n1b_l;
    Integer k=0;
    for(int i=start; i<=end; i++)
    {
        sendArray_p[k] = uc(i);
        sendError_p[k] = err(i);
        sendxloc_p[k]  = x(i);
        k++;
    }

    Real *uFinal_p = new Real [NTot];
    Real *errorG_p = new Real [NTot];
    Real *xloc_p = new Real [NTot];

    MPI_Gatherv(sendArray_p, sendCount, MPI_DOUBLE, 
                uFinal_p, recvCount_p, displacement_p, MPI_DOUBLE, root, MPI_COMM_WORLD);

    MPI_Gatherv(sendError_p, sendCount, MPI_DOUBLE, 
                errorG_p, recvCount_p, displacement_p, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    MPI_Gatherv(sendxloc_p, sendCount, MPI_DOUBLE, 
                xloc_p, recvCount_p, displacement_p, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    if(debug>0)
    {
        fprintf(debugFile, "uFinal: [\t");
        for(int i=0; i<NTot; i++)
        fprintf(debugFile, "%1.9e\t", uFinal_p[i]);
        fprintf(debugFile, "]\n");
        fflush(debugFile);
    }

    if(saveMatlab==2 && myRank == root)
    {   
        int t_nx, t_n1a, t_n1b;
        getLocalIndexBounds(np-1, np, nx, t_nx, t_n1a, t_n1b);
        FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
        fprintf(matlabFile, "%% File written by heat1Dmpi.C\n");
        fprintf(matlabFile,"xa=%g; xb=%g; kappa=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",xa,xb,kappa,tFinal,maxErr_g, cpuTimeStep_g);
        fprintf(matlabFile,"nx=%d; dx=%14.6e; numGhost=%d; n1a=%d; n1b=%d; nd1a=%d; nd1b=%d;\n",nx,dx,numGhost,n1a_l,t_n1b,nd1a_l,t_n1b+numGhost);
        fprintf(matlabFile,"solutionName=\'%s\';\n",solutionName.c_str());

        writeMatlabVector( matlabFile, xloc_p, "x", nd1a_l, t_n1b+numGhost );
        writeMatlabVector( matlabFile, uFinal_p, "u", nd1a_l, t_n1b+numGhost );
        writeMatlabVector( matlabFile, errorG_p, "err", nd1a_l, t_n1b+numGhost );
        fflush(matlabFile);
        fclose(matlabFile);
        fprintf(stdout, "Wrote file %s\n\n", matlabFileName.c_str());
    }

    if (debugFile)
        fclose(debugFile);

    MPI_Finalize();

    // delete pointers
    delete [] x_l;
    delete [] boundaryCondition_l;
    delete [] u_l[0];
    delete [] u_l[1];
    delete [] recvCount_p;
    delete [] displacement_p;
    delete [] err_p;
    delete [] sendArray_p;
    delete [] sendError_p;
    delete [] sendxloc_p;
    delete [] uFinal_p;
    delete [] errorG_p;
    delete [] xloc_p;

    // Free cuda memory
    FreeCudaMemory(dev_xa); FreeCudaMemory(dev_dx);
    FreeCudaMemory(dev_n1aL); FreeCudaMemory(dev_n1bL); FreeCudaMemory(dev_nd1L); FreeCudaMemory(dev_numGhost);
    FreeCudaMemory(dev_Nt); FreeCudaMemory(dev_x);
    FreeCudaMemory(dev_uc); FreeCudaMemory(dev_un);
    FreeCudaMemory(dev_rx); FreeCudaMemory(dev_t);
    FreeCudaMemory(dev_err);

    // undef macros
    #undef x
    #undef boundaryCondition
    #undef UTRUE
    #undef UTRUEX
    #undef UTRUET
    #undef UTRUEXX
    #undef FORCE
    #undef uc
    #undef un
    
    return 0;
}