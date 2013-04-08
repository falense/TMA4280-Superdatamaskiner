/*
  C-program to solve the two-dimensional Poisson equation on 
  a unit square using one-dimensional eigenvalue decompositions
  and fast sine transforms

  einar m. ronquist
  ntnu, october 2000
  revised, october 2001
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <omp.h>

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

Real f_funk(int x, int y, Real h, Real pi){
	
	Real ix, iy;
	ix = (x+1)*h;
	iy = (y+1)*h;
	//double r = exp(x)*sin(2*PI*x)*sin(PI*y);
	Real r = 5*pi*pi*sin(pi*ix)*sin(2*pi*iy);
	return r;
}
Real g_funk(int x, int y, Real h, Real pi){
	Real ix, iy;
	ix = (x+1)*h;
	iy = (y+1)*h;
	Real r = sin(pi*ix)*sin(2*pi*iy);
	return r;
}
Real absolute(Real i){
	if (i < 0.0) return -i;
	else return i;
}
Real wallTime(){
	struct timeval t;
	double msec = 0;
	gettimeofday(&t,0x0);
	msec = t.tv_usec / 1000.0 + t.tv_sec*1000.0;
	
	return msec;

}

main(int argc, char **argv )
{
  Real *diag, **b, **bt, **z;
  Real pi, h, error;
  int i, j, n, m, nn;

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

 if( argc < 2 ) {
    printf("need a problem size\n");
	return 1;
  }
  
  
  n  = atoi(argv[1]);
  m  = n-1;
  nn = 4*n;

  diag = createRealArray (m);
  b    = createReal2DArray (m,m);
  bt   = createReal2DArray (m,m);
  
  double startTime, endTime;
  startTime = omp_get_wtime();
  z = (Real**)malloc(sizeof(Real*)*omp_get_max_threads());
  
#pragma omp parallel
{
  //printf("Thread %d\n",(int)omp_get_thread_num());//omp_get_num_threads());
  int index = omp_get_thread_num();
  z[index]    = createRealArray (nn);
}
  
  
  h    = 1./(Real)n;
  pi   = 4.*atan(1.);
  #pragma omp parallel for schedule(static)
  for (i=0; i < m; i++) {
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  }
  #pragma omp parallel for schedule(static)
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = h*h*f_funk(j,i,h,pi);
    }
  }
  /*
  End init
  
  Step 1
  bt = (Qt)b(Q)
  */
  #pragma omp parallel for schedule(static)
  for (j=0; j < m; j++) {
    fst_(b[j], &n, z[omp_get_thread_num()], &nn);
  }
  
  transpose (bt,b,m);
  
  #pragma omp parallel for schedule(static)
  for (i=0; i < m; i++) {
    fstinv_(bt[i], &n, z[omp_get_thread_num()], &nn);
  }
  /*
  Step 2
  
  bt[j][i] = bt[j][i]/(diag[i]+diag[j])
  
  */
  #pragma omp parallel for schedule(static)
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] /= (diag[i]+diag[j]);
    }
  }
  /*
	End Step 2
	
	
	
	Step 3
	b = (Qt)bt(Q)
*/
  #pragma omp parallel for schedule(static)
  for (i=0; i < m; i++) {
    fst_(bt[i], &n, z[omp_get_thread_num()], &nn);
  }
  
  
  transpose (b,bt,m);
  
  #pragma omp parallel for schedule(static)
  for (j=0; j < m; j++) {
    fstinv_(b[j], &n, z[omp_get_thread_num()], &nn);
  }
/*
	End step 3
	
	Calculate deviation:
	dev = max(b)
	
*/
  endTime = omp_get_wtime();
  
  error = 0.0;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
		Real diff = b[j][i] - g_funk(j,i,h,pi);
		if (error < absolute(diff)) error = absolute(diff);//abs(diff);
		//if (m < 16)printf("\t%f",diff);
    }
	//if (m < 16)printf("\n");
  }
  printf (" error = %f \n",error);
  
  
  printf (" Time spent calculating %f msec \n",(endTime-startTime)*1000);
  
  return 0;
}

void transpose (Real **bt, Real **b, int m)
{
  int i, j;
  
  #pragma omp parallel for schedule(static)
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = b[i][j];
    }
  }
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}
