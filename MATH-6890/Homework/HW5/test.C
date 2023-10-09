#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

int main( int argc, char *argv[] )
{
  // Convert the arguments into input array
  int N = 20;
  if (argc > 1)
    N = atoi(argv[1]);
  int *a_p[2];
  a_p[0] = new int [N];
  a_p[1] = new int [N];
#define a(i1, i2) a_p[(i1)][(i2)]
  int i, j;
  int ia = 0;
  int sum=0;

#pragma omp parallel  num_threads(4) shared(a_p, ia, N) private(i,j) reduction(+:sum)
  { sum = 0;
#pragma omp for collapse(2)
    for (j=ia; j<N; j++)
    for (i=ia; i<2; i++)
    {
      //printf("i: %d\n", i);
      //printf("a_p[%d] is set by thread %d \n", i, omp_get_thread_num());
      a(i, j) = i*j+10;
      sum += a(i, j)%2;
    }
    printf("local sum at thread %d is:  %d\n", omp_get_thread_num(), sum);
  }
  //for ( i=0; i<N; i++)
  //{
  //  printf("a[%d]: %d\n", i, a_p[i]);
  //}
  printf("global sum is: %d\n", sum);
  delete [] a_p[0];
  delete [] a_p[1];
  return 0;
}
