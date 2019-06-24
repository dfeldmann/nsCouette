#include"TayC.h"

static cublasHandle_t cublasHandle;
static double2 alpha[1];
static double alphad[1];
static int* infoArray;
static double *dirM, *invM;
double *Ah[NT*NZ];
double *Bh[NT*NZ];
double** dirMp,**invMp;
int* pivotArray;

void cublasCheck(cublasStatus_t error, const char* function )
{
  if(error !=  CUBLAS_STATUS_SUCCESS)
  {
    printf("\n error  %s : %d \n", function, error);
    exit(1);
  }
    
  return;
}  

void setCublas(void){

    cublasCheck(cublasCreate(&cublasHandle),"Cre");


    CHECK_CUDART(cudaMalloc((void**)&infoArray,NT*NZ*sizeof(int)));

    alpha[0].x=1.0;
    alpha[0].y=0.0;

    alphad[0]=1.0;

    //Set to invert infinity matrix
        
    CHECK_CUDART(cudaMalloc((void **)&dirM,8*8*NZ*NT*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void **)&invM,8*8*NZ*NT*sizeof(double)));
    
    CHECK_CUDART(cudaMalloc((void **)&dirMp,NZ*NT*sizeof(double*)));
    CHECK_CUDART(cudaMalloc((void **)&invMp,NZ*NT*sizeof(double*)));
    
    CHECK_CUDART(cudaMalloc((void **)&pivotArray,NZ*NT*sizeof(int)));

    
    for(int i=0;i<NT*NZ;i++){
        Ah[i] = &dirM[i*8*8];
        Bh[i] = &invM[i*8*8];
    }

    CHECK_CUDART(cudaMemcpy(dirMp,Ah,NT*NZ*sizeof(double*),cudaMemcpyHostToDevice));
    CHECK_CUDART(cudaMemcpy(invMp,Bh,NT*NZ*sizeof(double*),cudaMemcpyHostToDevice));
    
    return;
    }

void transpose_A(double2* u_2,double2* u_1){

	//Transpuesta de [i,k,j][NZ,NT,NR] a -----> [j,i,k][NR,NZ,NT]

	cublasCheck(cublasZgeam(cublasHandle,CUBLAS_OP_T,CUBLAS_OP_T,NZ*NT,NR,
                            alpha,(const double2*)u_1,NR,0,0,NR,u_2,NZ*NT),"Tr");
	return;


}

void transpose_infforward(double* u_2, double* u_1,int mr){

	//Transpuesta de [j,i,k][mr,NZ,NT] a -----> [i,k,j][NZ,NT,mr]

	cublasCheck(cublasDgeam(cublasHandle,CUBLAS_OP_T,CUBLAS_OP_T,mr,NZ*NT,
                            alphad,(const double*)u_1,NZ*NT,0,0,NZ*NT,u_2,mr),"Tr");
	return;

}

void transpose_infbackward(double* u_2,double* u_1,int mr){

	//Transpuesta de [i,k,j][NZ,NT,mr] a -----> [j,i,k][mr,NZ,NT]

	cublasCheck(cublasDgeam(cublasHandle,CUBLAS_OP_T,CUBLAS_OP_T,NZ*NT,mr,
                            alphad,(const double*)u_1,mr,0,0,mr,u_2,NZ*NT),"Tr");
	return;


}

void invert_infmatrix(double* src){


    transpose_infforward(dirM,src,8*8);
    
    cublasCheck(cublasDmatinvBatched(cublasHandle,8,(const double**)dirMp,8,
                   (double**)invMp,8,infoArray,NT*NZ),"Tr");
     
    transpose_infbackward(src,invM,8*8);


     return;
}
    
