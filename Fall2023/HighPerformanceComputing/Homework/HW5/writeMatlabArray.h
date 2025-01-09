#ifndef WRITE_MATLAB_ARRAY_H
#define WRITE_MATLAB_ARRAY_H

// ------------------------------------------------------------------------------------
// Save an array to a matlab file.
// ------------------------------------------------------------------------------------
int writeMatlabArray( FILE *matlabFile, RealArray & u, const char *name,
		      int numberOfComponents, IntegerArray & dimension, int nThreads )
{
  const int nd1a = dimension(0,0);
  const int nd1b = dimension(1,0);
  const int nd2a = dimension(0,1);
  const int nd2b = dimension(1,1);
  const int nd1 = nd1b-nd1a+1;
  const int nd2 = nd2b-nd2a+1;

  const int numPerLine=8; // number of entries per line
  if( numberOfComponents==1 )
  {
    fprintf(matlabFile,"%s=zeros(%d,%d);\n",name,nd1,nd2);

    // Real *ucp_p = u.getDataPointer();
    //#define ucp( id1, id2) ucp_p[ (id1) + nd1*(id2) ]
    
    int count=0;
    int i1, i2;
    #pragma omp parallel num_threads(nThreads) shared(u, name, count, matlabFile) private(i1, i2)
    {
      #pragma omp for collapse(2)
      for (i2=nd1a; i2<=nd1b; i2++)
      for (i1=nd2a; i1<=nd2b; i1++)
      {
	fprintf(matlabFile, "%s(%3d,%3d)=%12.5e; ",name,i1-nd1a+1,i2-nd2a+1,u(i1,i2));
	if( count % numPerLine == numPerLine-1 )
	  fprintf(matlabFile,"\n"); // new line
	count++;
      }
    }
    //for( int i2=nd1a; i2<=nd1b; i2++ )
    //for( int i1=nd1a; i1<=nd1b; i1++ )
    //{
    //  fprintf(matlabFile,"%s(%3d,%3d)=%12.5e; ",name,i1-nd1a+1,i2-nd2a+1,u(i1,i2));
    //  if( count % numPerLine == numPerLine-1 )
    //	fprintf(matlabFile,"\n"); // new line
    //  count++;
    //}
    fprintf(matlabFile,"\n");
  }
  else
  {
    fprintf(matlabFile,"%s=zeros(%d,%d,%d);\n",name,nd1,nd2,numberOfComponents);

    int count=0;
    int i1, i2, m;
 
    #pragma omp parallel num_threads(nThreads) shared(u, numberOfComponents, matlabFile, name, count) private(m, i1, i2)
    {
      #pragma omp for collapse(3)
      for (m=0; m<numberOfComponents; m++)
      for (i2=nd1a; i2<=nd1b; i2++)
      for (i1=nd2a; i1<=nd2b; i1++ )
      {
	fprintf(matlabFile,"%s(%3d,%3d,%d)=%12.5e; ",name,i1-nd1a+1,i2-nd2a+1,m+1,u(i1,i2,m));
	if( count % numPerLine == numPerLine-1 )
          fprintf(matlabFile,"\n"); // new line                                                                                                                                        
        count++;
      }
    }
    
    //for( int m=0; m<numberOfComponents; m++ )
    //for( int i2=nd1a; i2<=nd1b; i2++ )
    //for( int i1=nd1a; i1<=nd1b; i1++ )
    //{
    //  fprintf(matlabFile,"%s(%3d,%3d,%d)=%12.5e; ",name,i1-nd1a+1,i2-nd2a+1,m+1,u(i1,i2,m));
    //  if( count % numPerLine == numPerLine-1 )
    //	fprintf(matlabFile,"\n"); // new line
    //  count++;
    //}
    fprintf(matlabFile,"\n");

  }
  return 0;
}

#endif
