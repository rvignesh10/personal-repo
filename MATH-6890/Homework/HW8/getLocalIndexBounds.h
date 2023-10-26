#ifndef GET_LOCAL_INDEX_BOUNDS_H
#define GET_LOCAL_INDEX_BOUNDS_H

// =======================================================================
// Distribute an index dimension across npx processors:
// rank (input) : rank of process, 0<= rank < npx
// npx (input) : distribute across this many processors
// nx (input) : total number of grid cells on global grid
// nx_l, n1a_l, n1b_l (output) : local values
//
// Note: nx+1 = total grid points across all processors
// ========================================================================
int getLocalIndexBounds( const int rank, const int npx, const int nx,
                        int & nx_l, int & n1a_l, int & n1b_l )
{
// ---- Example : distribute grid points ------
// Global:
//  X--X--X--X--X--X--X--X--X--X--X
//  0  1  2  3  4  5  6  7  8  9  10 nx=10
// n1a                           n1b
// Local: (rank=0, np=2)
//    X--X--X--X--X--X
//    0  1  2  3  4  5 nx_l=5
//   n1a_l          n1b_l
// Local: (rank=1, np=2)
//    X--X--X--X--X
//    6  7  8  9  10 nx_l=4
//  n1a_l         n1b_l
//
    nx_l = (nx+1)/npx -1; // nx_l+1 = (nx_1+1)/np
    n1a_l = (nx_l+1)*rank; // local starting index

    // There may be extra points if (nx+1) is not a multiple of np
    int extra = nx+1- npx*(nx_l+1);
    if( rank<extra )
    { // add one extra point to proc's on the left side
        nx_l += 1;
        n1a_l += rank;
    }
    else
    {
        n1a_l += extra;
    }
    n1b_l = n1a_l+nx_l; // local end index

    return 0;
}

#endif