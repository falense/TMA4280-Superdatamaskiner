#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>


typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

int size, rank, partSize, startIndex, endIndex, maxPartSize = 0;
int *rankToStartIndex, *rankToEndIndex;

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

void Partition_problem(int problemSize){

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	//if (rank == 0) printf("Attempting to split the work across %d threads\n",size);
	partSize = (problemSize)/size;
	int remainder = (problemSize)-partSize*size;
	
	//if (rank == 0) printf("partSize =  %d\n",partSize);
	
	rankToStartIndex = (int*)malloc(sizeof(int)*size);
	rankToEndIndex = (int*)malloc(sizeof(int)*size);
	int i;
	for (i = 0; i < size; i++){
		if (size -1  != i){
			if(remainder > i){
			
				rankToStartIndex[i] = partSize*i+i; 
				rankToEndIndex[i] = partSize*(i+1)+i+1;
			}
			else{
			
				rankToStartIndex[i] = partSize*i+remainder; 
				rankToEndIndex[i] = partSize*(i+1)+remainder;
			}
		}
		else{
			rankToStartIndex[i] = partSize*i+remainder; 
			rankToEndIndex[i] = problemSize;
		}
		if (rankToEndIndex[i]-rankToStartIndex[i] > maxPartSize){
			maxPartSize = rankToEndIndex[i]-rankToStartIndex[i];
		}
		//if (rank == 0)printf("Rank %d, (start,end) = (%d,%d), range = %d\n",i,rankToStartIndex[i],rankToEndIndex[i],rankToEndIndex[i]-rankToStartIndex[i]);
	}
	startIndex = rankToStartIndex[rank];
	endIndex = rankToEndIndex[rank];
	
	//printf("Thread %d: Working part %d -> %d, maxPartSize = %d\n",rank,startIndex,endIndex,maxPartSize);
	
}

void printMatrix(Real **b, int m){
	if (m > 16) return;
	int i, j,curRank;
	if (rank == 0){
		for (i = startIndex; i < endIndex; i++){
			for (j = 0; j < m; j++){
				printf("\t %10.6f",b[i][j]);
			}
			printf("\n");
		}
		Real *buffer = (Real*)malloc(sizeof(Real)*m);
		for( curRank = 1; curRank < size; curRank++){
			for (i = rankToStartIndex[curRank]; i < rankToEndIndex[curRank]; i++){
				MPI_Recv(buffer, m, MPI_DOUBLE, curRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (j = 0; j < m; j++){
					printf("\t %10.6f",buffer[j]);
					buffer[j] = 0.0;
				}
				
				printf("\n");
			}
		}
		free(buffer);
		printf("\n");
	}
	else{
		for (i = startIndex; i < endIndex; i++){
			
			Real *buffer = (Real*)malloc(sizeof(Real)*m);
			memset(buffer,0,m*sizeof(Real));

			MPI_Send(b[i], m, MPI_DOUBLE, 0, 0,MPI_COMM_WORLD);
			free(buffer);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

main(int argc, char **argv )
{
	Real *diag, **b, **bt, *z;
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


	// Distribute the columns or rows between the different threads
	MPI_Init(&argc, &argv);
	Partition_problem(m);
	MPI_Barrier(MPI_COMM_WORLD);


	diag = createRealArray (m);
	b    = createReal2DArray (m,m);
	bt   = createReal2DArray (m,m);

	double startTime, endTime;
	startTime = MPI_Wtime();



	z    = createRealArray (nn);


	h    = 1./(Real)n;
	pi   = 4.*atan(1.);
	
	for (i=0; i < m; i++) {
		diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
	}
	
	for (j=startIndex; j < endIndex; j++) {
		for (i=0; i < m; i++) {
			b[j][i] = h*h*f_funk(j,i,h,pi);
		}
	}
	/*
End init

Step 1
bt = (Qt)b(Q)
*/
	for (j=startIndex; j < endIndex; j++) {
		fst_(b[j], &n, z, &nn);
	}

	transpose (bt,b,m);

	for (i=startIndex; i < endIndex; i++) {
		fstinv_(bt[i], &n, z, &nn);
	}
	/*
Step 2

bt[j][i] = bt[j][i]/(diag[i]+diag[j])

*/
	for (j=startIndex; j < endIndex; j++) {
		for (i=0; i < m; i++) {
			bt[j][i] /= (diag[i]+diag[j]);
		}
	}
	/*
	End Step 2
	
	Step 3
	b = (Qt)bt(Q)
*/
	for (i=startIndex; i < endIndex; i++) {
		fst_(bt[i], &n, z, &nn);
	}

	transpose (b,bt,m);
	
	for (j=startIndex; j < endIndex; j++) {
		fstinv_(b[j], &n, z, &nn);
	}
	/*
	End step 3
	
	Calculate deviation:
	dev = max(b)
	
*/
	endTime = MPI_Wtime();
	error = 0.0;
	for (j=startIndex; j < endIndex; j++) {
		for (i=0; i < m; i++) {
			Real diff = b[j][i] - g_funk(j,i,h,pi);
			if (error < absolute(diff)) error = absolute(diff);//abs(diff);
			
		}
	}
	if(rank != 0){
		MPI_Send(&error, 1, MPI_DOUBLE, 0, 0,MPI_COMM_WORLD);
	}
	else{
		for (i = 1; i < size; i++){
			double otherError = 0;
			MPI_Recv(&otherError, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (otherError > error) error = otherError;
			//printf("%f\n",otherError);
		}
		printf (" error = %f \n",error);
	}




	if (rank == 0)printf (" Time spent calculating %f msec \n",(endTime-startTime)*1000);

	MPI_Finalize();
	return 0;
}

void transpose (Real **bt, Real **b, int m)
{
	int i, j,dstRank,srcRank;

	int sendrecvcount = maxPartSize*maxPartSize;
	//if (rank == 0) printf("Entering transpose operation. Send recv buf size is %d * %d\n",sendrecvcount,size);
	Real * sendbuf = (Real*)malloc(sizeof(Real)*sendrecvcount*size);
	Real * recvbuf = (Real*)malloc(sizeof(Real)*sendrecvcount*size);
	int state = 0;
	//printf("Thread %d: state %d, maxBufIndex = %d \n",rank,state++,sendrecvcount*size);
	for (dstRank = 0; dstRank < size; dstRank++){
		for (j=0; j < endIndex-startIndex; j++){
			for (i=0; i < rankToEndIndex[dstRank]-rankToStartIndex[dstRank]; i++) {
				int index = dstRank*sendrecvcount+j*maxPartSize+i;
				sendbuf[index] = b[j+startIndex][i+rankToStartIndex[dstRank]];
			}
		}
	}

	//printf("Thread %d: state %d test\n",rank,state++);
	//int MPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
	int *sendcounts = (int*)malloc(sizeof(int)*size);
	int *displacements = (int*)malloc(sizeof(int)*size);
	for ( dstRank = 0; dstRank < size; dstRank++){
		sendcounts[dstRank] = (rankToEndIndex[dstRank]-rankToStartIndex[dstRank])*(rankToEndIndex[rank]-rankToStartIndex[rank]);
		displacements[dstRank] = dstRank*sendrecvcount;
		//if (rank == 0) printf("Send info rank %d: (%d,%d)\n",rank,sendcounts[dstRank],displacements[dstRank]);
	}

	//printf("Thread %d: state %d \n",rank,state++);
	MPI_Alltoallv(sendbuf, sendcounts, displacements, MPI_DOUBLE, recvbuf, sendcounts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);


	//printf("Thread %d: Passed the mpi call. State %d \n",rank,state++);
	for (srcRank = 0; srcRank < size; srcRank++){
		for (j=0; j < endIndex-startIndex; j++){
			for (i=0; i < rankToEndIndex[srcRank]-rankToStartIndex[srcRank]; i++) {
				int index = srcRank*sendrecvcount+i*maxPartSize+j;
				//printf("(i,j) = (%d,%d), index = %d, srcRank = %d\n",i,j,index,srcRank);
				bt[j+startIndex][i+rankToStartIndex[srcRank]] = recvbuf[index];
				//if (rank == 1) printf("Setting (%d,%d) to %f \n",j+startIndex,i+rankToStartIndex[srcRank],recvbuf[index]);
			}
		}
	}
	free(sendbuf);
	free(recvbuf);
	MPI_Barrier(MPI_COMM_WORLD);
	//printMatrix(b,m);
	//printMatrix(bt,m);
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

	for (i=startIndex; i < endIndex; i++) {
		a[i] = (Real  *)malloc(n2*sizeof(Real));
		memset(a[i],0,n2*sizeof(Real));
	}
	return (a);
}
