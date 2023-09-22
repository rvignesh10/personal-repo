// ==================================================
// Routines for solving tridiagonal linear systems
// ==================================================
#include "A++.h"
#include <assert.h>

typedef double Real;
typedef doubleSerialArray RealArray;

// define macros for array references

// =====================================================================
// Factror a tridiagonal matrix -- the factorization is stored in Ax
// =====================================================================

int factorTridiagonalMatrix ( RealArray &Ax ) {
    const int iax = Ax.getBase(1), ibx = Ax.getBound(1);

    Real *Ax_p = Ax.getDataPointer();
    #define AX(i, j) Ax_p[(i) + 3*(j)]

    // Ax.display("Ax before factor");
    
    // Factor: (no pivoting)
    //
    // [b0 c0                    ]
    // [a1 b1 c1                 ]
    // [   a2 b2 c2              ]
    // [      a3 b3 c3           ]
    // [            ....         ]
    // [                am bm cm ]
    // [                   an bn ]    

    // Note: AX is a tri-diagonal matrix /
    // AX(0, :) - [0  a1 a2 ... am an]' vector
    // AX(1, :) - [b0 b1 b2 ... bm bn]' vector
    // AX(2, :) - [c0 c1 c2 ... cm  0]' vector

    for (int i1=iax+1; i1<=ibx; i1++) {
        Real d = -AX(0, i1)/AX(1, i1-1); // -a[i1]/b[i1-1]
        AX(1, i1) += d*AX(2, i1-1);
        AX(0, i1) = d;
    }
    // Ax.display("Ax after factor");

    return 0;
    #undef AX
}

int factorTridiagonalMatrix(Real *Ax_p, const int &num_row, const int &n1a, const int &n1b) {
    assert(num_row == 3);
    #define AX(i, j) Ax_p[(i) + 3*((j)-n1a)]
    // Factor: (no pivoting)
    //
    // [b0 c0                    ]
    // [a1 b1 c1                 ]
    // [   a2 b2 c2              ]
    // [      a3 b3 c3           ]
    // [            ....         ]
    // [                am bm cm ]
    // [                   an bn ]    

    // Note: AX is a tri-diagonal matrix /
    // AX(0, :) - [0  a1 a2 ... am an]' vector
    // AX(1, :) - [b0 b1 b2 ... bm bn]' vector
    // AX(2, :) - [c0 c1 c2 ... cm  0]' vector

    for (int i1=n1a+1; i1<=n1b; i1++) {
        Real d = -AX(0, i1)/AX(1, i1-1); // -a[i1]/b[i1-1]
        AX(1, i1) += d*AX(2, i1-1);
        AX(0, i1) = d;                   // save d here
    }

    return 0;
    #undef AX
}

// =====================================================================
// Solve the tridiagonal matrix problem given the factored matrix Ax
// =====================================================================

int solveTridiagonal( RealArray & Ax, RealArray & rhs )
{
    const int iax=Ax.getBase(1), ibx=Ax.getBound(1);

    const Real *Ax_p = Ax.getDataPointer();
    #define AX(i, j) Ax_p[(i) + 3*(j)]
    Real *rhs_p = rhs.getDataPointer();
    #define RHS(i) rhs_p[(i)]

    // --- forward elimination ---
    for( int i1=iax+1; i1<=ibx; i1++ )
    {
        RHS(i1) += AX(0,i1)*RHS(i1-1);
    }

    // --- back-substitution ---
    // [b0 c0                   ][x ]   [x0]
    // [   b1 c1                ][x ]   [x1]
    // [      b2 c2             ][x ]   [x2]
    // [         b3 c3          ][x ] = [x3]
    // [            ...         ][x ]   [ ]
    // [                  bm cm ][xm]   [rm]
    // [                     cn ][xn]   [rn]
    RHS(ibx) = RHS(ibx)/AX(1,ibx);
    for( int i1=ibx-1; i1>=iax; i1-- )
    {
        RHS(i1) = (RHS(i1) - AX(2,i1)*RHS(i1+1) )/AX(1,i1);
    }

    return 0;
    #undef AX
    #undef RHS
}

int solveTridiagonal(Real *Ax_p, Real *rhs_p, const int &num_row, const int &n1a, const int &n1b){
    assert(num_row == 3);
    #define AX(i, j) Ax_p[(i) + 3*((j)-n1a)]
    #define RHS(i) rhs_p[(i)-n1a]

    // --- forward elimination ---
    for( int i1=n1a+1; i1<=n1b; i1++ ) {
        RHS(i1) += AX(0,i1)*RHS(i1-1);
    }

    // --- back-substitution ---
    // [b0 c0                   ][x ]   [x0]
    // [   b1 c1                ][x ]   [x1]
    // [      b2 c2             ][x ]   [x2]
    // [         b3 c3          ][x ] = [x3]
    // [            ...         ][x ]   [ ]
    // [                  bm cm ][xm]   [rm]
    // [                     cn ][xn]   [rn]
    RHS(n1b) = RHS(n1b)/AX(1,n1b);
    for( int i1=n1b-1; i1>=n1a; i1-- ) {
        RHS(i1) = (RHS(i1) - AX(2,i1)*RHS(i1+1) )/AX(1,i1);
    }
    
    return 0;
    #undef AX
    #undef RHS
}