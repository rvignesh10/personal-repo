#ifndef WRITE_MATLAB_VECTOR_H
#define WRITE_MATLAB_VECTOR_H
// ---------------------------------------------------------------------------------------
// Function to save a vector to a matlab file.
// matlabFile (input) : save vector to this file
// u_p (input) : array of vector values
// name (input) : name for array
// (nd1a:nd1b) (input) : array dimensions
// ---------------------------------------------------------------------------------------
typedef double Real;
int writeMatlabVector( FILE *matlabFile, Real *u_p, const char *name, int nd1a, int nd1b ) {
    #define u(i) u_p[i-nd1a]

    const int numPerLine=8; // number of entries per line
    // Save the vector as:
    // name = [ num num num num num ...
    //          num num num num num ];
    fprintf(matlabFile,"%s=[",name);
    for( int i=nd1a; i<=nd1b; i++ ) {
        fprintf(matlabFile,"%20.15e ",u(i));
        if( (i-nd1a) % numPerLine == numPerLine-1 )
            fprintf(matlabFile,"...\n"); // continuation line
    }
    fprintf(matlabFile,"];\n");

    return 0;
    #undef u
}

// ------------------------------------------------------------------------------------
// Save a vector to a matlab file.
// ------------------------------------------------------------------------------------
int writeMatlabVector( FILE *matlabFile, RealArray & u, const char *name, int nd1a, int nd1b ) {
    const int numPerLine=8; // number of entries per line
    fprintf(matlabFile,"%s=[",name);
    for( int i=nd1a; i<=nd1b; i++ ) {
        fprintf(matlabFile,"%20.15e ",u(i));
        if( (i-nd1a) % numPerLine == numPerLine-1 )
            fprintf(matlabFile,"...\n"); // continuation line
    }
    fprintf(matlabFile,"];\n");

    return 0;
}
#endif