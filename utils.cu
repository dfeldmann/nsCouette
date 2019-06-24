#include "TayC.h"

static __global__ void normalize_kernel(double2* u, double norm, size_t elements)
{
    
    int h = blockIdx.x * blockDim.x + threadIdx.x;

    if (h<elements)
    {
          
      double2 uh=u[h];
      
      uh.x=uh.x/norm;
      uh.y=uh.y/norm;
        
      u[h]=uh;
  }
  
}

static __global__ void decoupleForward_kernel(double2* ur,double2* ut)
{
    
   int h = blockIdx.x * blockDim.x + threadIdx.x;
  
   if (h<NR*NT*NZ)
   {
    
    double2 ur_l,ut_l;
    double2 up_l,um_l;
    
    ur_l=ur[h];
    ut_l=ut[h];
    
    up_l.x=  ur_l.x - ut_l.y;
    up_l.y=  ur_l.y + ut_l.x;
    
    um_l.x=  ur_l.x + ut_l.y;
    um_l.y=  ur_l.y - ut_l.x;
    
    ur[h]=up_l;
    ut[h]=um_l;
        
  }
  
}

static __global__ void decoupleBackward_kernel(double2* up,double2* um)
{
    
   int h = blockIdx.x * blockDim.x + threadIdx.x;
  
   if (h<NR*NT*NZ)
   {
    
    double2 ur_l,ut_l;
    double2 up_l,um_l;

    up_l=up[h];
    um_l=um[h];
    
    ur_l.x=  0.5*(up_l.x + um_l.x);
    ur_l.y=  0.5*(up_l.y + um_l.y);
    
    ut_l.x= -0.5*(um_l.y - up_l.y);
    ut_l.y=  0.5*(um_l.x - up_l.x);
    
    up[h]=ur_l;
    um[h]=ut_l;
    
  }
}

static __global__ void fieldAbs_kernel(double2* ur_diff,double2* ur)
{
    
   int h = blockIdx.x * blockDim.x + threadIdx.x;
  
   if (h<NR*NT*NZ)
   {
    
    double2 urd_l,urr_l;
    
    urr_l=ur[h];
    
    urd_l.x= abs(urr_l.x);
    urd_l.y= abs(urr_l.y);
    
    
    ur_diff[h]=urd_l;
        
  }
  
}

static __global__ void diffAbs_kernel(double2* ur_diff,double2* ur)
{
    
   int h = blockIdx.x * blockDim.x + threadIdx.x;
  
   if (h<NR*NT*NZ)
   {
    
    double2 urd_l,urr_l;
    
    urd_l=ur_diff[h];
    urr_l=ur[h];
    
    urd_l.x= abs(urd_l.x - urr_l.x);
    urd_l.y= abs(urd_l.y - urr_l.y);
    
    
    ur_diff[h]=urd_l;
        
  }
  
}

static __global__ void sumnonlinear_kernel(double2* n1,double2* n2, double d_implicit)
{
    
   int h = blockIdx.x * blockDim.x + threadIdx.x;
  
   if (h<NR*NT*NZ)
   {
    
    double2 u1,u2;
    

    u1=n1[h];
    u2=n2[h];
    
    u1.x=  d_implicit*u1.x + (1.0 - d_implicit)*u2.x;
    u1.y=  d_implicit*u1.y + (1.0 - d_implicit)*u2.y;
    
    n1[h]=u1;
    
  }
  
}


static __global__ void setzero_kernel(double2* u)
{
  
    int h = blockIdx.x * blockDim.x + threadIdx.x;
    int k = h%NZ;
    int j = h/(NT*NZ);
    int i = (h-j*NT*NZ)/NZ;
    
    if (j<NR && i<NT && k<NZ)
    {
        
      double2 zero={0.0,0.0};

      if((abs(i)+abs(k))==0){
        u[h]=zero;
      }

    }
}

static __global__ void base_kernel(double* ur, double* ut, double* uz, 
                                   double* mesh, double eta)
{
  
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (j<NR && k<2*NT*NZ)
    {
      
      
      int h=j*(2*NT*NZ)+k;
      
      double rad = mesh[j];
      double r_i = mesh[0];
      double r_o = mesh[NR-1];
      
      double ur_h, ut_h, uz_h;
      double zero=0.0;
      
      ur_h=zero;
      uz_h=zero;

      double mu=(REYNOLDS_OUTER)/(REYNOLDS_INNER);
      
      double C_1 = (REYNOLDS_INNER/r_o - REYNOLDS_OUTER/r_i)/(r_i/r_o - r_o/r_i);
      double C_2 = (REYNOLDS_INNER*r_o - REYNOLDS_OUTER*r_i)/(r_o/r_i - r_i/r_o);
      
      ut_h = C_1*rad + C_2/rad;
      
      ur[h]=ur_h;
      ut[h]=ut_h;
      uz[h]=uz_h;

//       if(j==NR-1) printf("\nout=%e",ut_h);
//       if(j==0) printf("\nin=%e",ut_h);
    }
}

static __global__ void perturb_kernel(double2* ur, double2* ut, double2* uz, 
                                      double* mesh,double epsilon)
{
  
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (j<NR && k<NT*NZ)
    {
      
      int h=j*NT*NZ+k;
      
      int k0=k%NZ;
      int t0=k/NZ;
      
      double rad = mesh[j];
      
      if(k0==1 && t0==1){
        ut[h].x=epsilon*(1.0-cos(PI2*(rad - mesh[0])/(mesh[NR-1]-mesh[0])));
        ur[h].x=epsilon*(1.0-cos(PI2*(rad - mesh[0])/(mesh[NR-1]-mesh[0])));

        //printf("\nut[%d]=%e",j,ut[h].x);
        ut[h].y=0.0;
        ur[h].y=0.0;
        
      }
    }
}

void normalize(double2* u,double norm,size_t elements){
    
    dim3 grid,block;
    block.x = 128;
    grid.x = (elements + block.x - 1)/block.x;

    normalize_kernel<<<grid,block>>>(u, norm, elements);

    return;

}

void decoupleForward(vfield u){
    
    dim3 grid,block;
    block.x = 128;
    grid.x = (NR*NT*NZ + block.x - 1)/block.x;

    decoupleForward_kernel<<<grid,block>>>(u.r,u.t);
    
    return;
    
}

void decoupleBackward(vfield u){
    
    dim3 grid,block;
    block.x = 128;
    grid.x = (NR*NT*NZ + block.x - 1)/block.x;

    decoupleBackward_kernel<<<grid,block>>>(u.r,u.t);

    return;
    
}

void fieldAbs(double2* u_diff, double2* u){
    
    dim3 grid,block;
    block.x = 128;
    grid.x = (NR*NT*NZ + block.x - 1)/block.x;

    fieldAbs_kernel<<<grid,block>>>(u_diff,u);


    return;
    
}

void diffAbs(double2* u_diff, double2* u){
    
    dim3 grid,block;
    block.x = 128;
    grid.x = (NR*NT*NZ + block.x - 1)/block.x;

    diffAbs_kernel<<<grid,block>>>(u_diff,u);


    return;
    
}

void updateNonlinear(vfield n1, vfield n2, double d_implicit)
{
    
    dim3 grid,block;
    block.x = 128;
    grid.x = (NR*NT*NZ + block.x - 1)/block.x;

    sumnonlinear_kernel<<<grid,block>>>(n1.t,n2.t,d_implicit);
    sumnonlinear_kernel<<<grid,block>>>(n1.r,n2.r,d_implicit);
    sumnonlinear_kernel<<<grid,block>>>(n1.z,n2.z,d_implicit);

    return;
    
}

void copyVfield(vfield u1, vfield u2){

    size_t size=NR*NT*NZ*sizeof(double2);

    CHECK_CUDART(cudaMemcpy(u2.r,u1.r,size,cudaMemcpyDeviceToDevice));
    CHECK_CUDART(cudaMemcpy(u2.t,u1.t,size,cudaMemcpyDeviceToDevice));
    CHECK_CUDART(cudaMemcpy(u2.z,u1.z,size,cudaMemcpyDeviceToDevice));


}

void copyBuffer(double2* u1, double2* u2){

    size_t size=NR*NT*NZ*sizeof(double2);

    CHECK_CUDART(cudaMemcpy(u2,u1,size,cudaMemcpyDeviceToDevice));

}

void onesGPU(double2* u1, size_t elements){
    
    
    size_t size=elements*sizeof(double2);

    double* ones=(double*)malloc(size);
    
    for(int i=0;i<2*elements;i++){
        ones[i]=1.0;
    }
    
    CHECK_CUDART(cudaMemcpy(u1,ones,size,cudaMemcpyHostToDevice));

    free(ones);
}

void zerosGPU(double2* u1, size_t elements){
    
    
    size_t size=elements*sizeof(double2);

    double* zeros=(double*)malloc(size);
    
    for(int i=0;i<2*elements;i++){
        zeros[i]=0.0;
    }
    
    CHECK_CUDART(cudaMemcpy(u1,zeros,size,cudaMemcpyHostToDevice));

    free(zeros);
}

void setZeromode(double2* u)
{

  //Sets mode k=0 to zero
  
  dim3 grid,block;
  block.x = 128;
  grid.x = (NR*NT*NZ + block.x - 1)/block.x;

  setzero_kernel<<<grid,block>>>(u);
  
}

void initField(vfield u, double* mesh)
{

  dim3 grid, block;
  block.x = 128;
  grid.x = (2*NT*NZ + block.x - 1)/block.x;
  grid.y = NR;
  
  double eta=mesh[0]/mesh[NR-1];
  
  double* mesh_d;
  CHECK_CUDART(cudaMalloc(&mesh_d,NR*sizeof(double)));
  CHECK_CUDART(cudaMemcpy(mesh_d, mesh, NR*sizeof(double),cudaMemcpyHostToDevice));
  
  base_kernel<<<grid,block>>>((double*)u.r, (double*)u.t, 
                              (double*)u.z, mesh_d, eta);
 
    
  fftForward_reduced(u.r);
  fftForward_reduced(u.t);
  fftForward_reduced(u.z);
  
  //Introduce a perturbation
  double epsilon=20.0;
  
  perturb_kernel<<<grid,block>>>((double2*)u.r, (double2*)u.t, 
                                 (double2*)u.z, mesh_d,epsilon);
  

  CHECK_CUDART(cudaFree(mesh_d));
  
}

