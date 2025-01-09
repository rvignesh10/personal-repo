#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H "tridiagonal.h"

// ==================================================
// Routines for solving tridiagonal linear systems
// ==================================================
#include "A++.h"
#include <assert.h>

typedef double Real;
typedef doubleSerialArray RealArray;

int factorTridiagonalMatrix(RealArray &Ax);
int factorTridiagonalMatrix(Real *Ax_p, const int &num_row, const int &n1a, const int &n1b);

int solveTridiagonal(RealArray &Ax, RealArray &rhs);
int solveTridiagonal(Real *Ax_p, Real *rhs_p, const int &num_row, const int &n1a, const int &n1b);

#endif