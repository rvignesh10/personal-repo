#ifndef WRITE_MATLAB_ARRAY_H
#define WRITE_MATLAB_ARRAY_H


int writeMatlabArray(FILE *matlabFile, double *u_p, std::string name, int dimension_p[][2]){
    #define dimension(side, axis) dimension_p[side][axis]
    const int nd1a = dimension(0,0);
    const int nd1b = dimension(1,0);
    const int nd2a = dimension(0,1);
    const int nd2b = dimension(1,1);
    const int nd1 = nd1b-nd1a+1;
    const int nd2 = nd2b-nd2a+1;
    #define u(i1, i2) u_p[ ((i1)-nd1a) + nd1*((i2)-nd2a) ]
    fprintf(matlabFile, "%s=[",name.c_str());
    for(int j=nd2a; j<=nd2b; j++)
    {
        for(int i=nd1a; i<=nd1b; i++)
        {   
            if(i==nd1b)
                fprintf(matlabFile,"%25.15e;", u(i, j));
            else
                fprintf(matlabFile,"%25.15e,", u(i,j));
        }
        fprintf(matlabFile,"\n");
    }
    fprintf(matlabFile,"];\n");

    #undef dimension
    #undef u
    
    return 0;
}

#endif