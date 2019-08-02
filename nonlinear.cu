/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file is part of nsCouette -- A high-performance code for direct         !
! numerical simulations of turbulent Taylor-Couette flow                       !
!                                                                              !
! Copyright (C) 2019 Marc Avila, Bjoern Hof, Jose Manuel Lopez, Markus Rampp,  !
!                    Liang Shi, Alberto Vela-Martin, Daniel Feldmann.          !
!                                                                              !
! nsCouette is free software: you can redistribute it and/or modify it under   !
! the terms of the GNU General Public License as published by the Free         !
! Software Foundation, either version 3 of the License, or (at your option)    !
! any later version.                                                           !
!                                                                              !
! nsCouette is distributed in the hope that it will be useful, but WITHOUT ANY !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    !
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more        !
! details.                                                                     !
!                                                                              !
! You should have received a copy of the GNU General Public License along with !
! nsCouette. If not, see <http://www.gnu.org/licenses/>.                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#include"TayC.h"

static __global__ void deriv_t_kernel(double2* u,double2* du)
{
  
   int h = blockIdx.x * blockDim.x + threadIdx.x;
   int k = h%NZP;
   int j = h/(NTP*NZP);
   int i = (h-j*NTP*NZP)/NZP;
  
  if (j<NR/NSTEPS_CONV && k<NZP && i<NTP)
  {
      
  double kz;
  double kt;

  kt=i<NTP/2 ? (double)i : (double)i-(double)NTP;
  kz=(double)k;

  //Fraction
  kz=(PI2/LZ)*kz;
  kt=(PI2/LT)*kt; 
  
  double2 uh=u[h];
  double2 duh;
  
  duh.x=  -kt*uh.y;
  duh.y=   kt*uh.x;
    
  du[h]=duh;

  }
}

static __global__ void deriv_z_kernel(double2* u,double2* du)
{
  
  int h = blockIdx.x * blockDim.x + threadIdx.x; 
  int k = h%NZP;
  int j = h/(NTP*NZP);
  int i = (h-j*NTP*NZP)/NZP;

  if (j<NR/NSTEPS_CONV && k<NZP && i<NTP)
  {

  double kz;
  double kt;

  kt=i<NTP/2 ? (double)i : (double)i-(double)NTP;
  kz=(double)k;

  //Fraction
  kz=(PI2/LZ)*kz;
  kt=(PI2/LT)*kt; 
  
  double2 uh=u[h];
  double2 duh;

  duh.x=  -kz*uh.y;
  duh.y=   kz*uh.x;

  du[h]=duh;
    
  
  }
}

static __global__ void padForward_kernel(double2* v,double2* u)
{

    int h = blockIdx.x * blockDim.x + threadIdx.x;
    int k = h%NZP;
    int j = h/(NTP*NZP);
    int i = (h-j*NTP*NZP)/NZP;
    
    if(h<NR/NSTEPS_CONV*NTP*NZP)
    {
    
        double factor=1.0;
        
        double2 aux;
        aux.x=0.0;
        aux.y=0.0;

        if((i<NT/2 || i>NT-1) && k<NZ){

            // X indices		
            int ip=i<NT/2 ? i : i-NT/2;
            
            // Z indices
            int kp=k;


            int h =j*NT*NZ+ip*NZ+kp;
            
            
            aux=u[h];
            
            aux.x*=factor;
            aux.y*=factor;

        }	
    
    int hp=j*NTP*NZP+i*NZP+k;
    v[hp]=aux;		
    }
    
}

static __global__ void padBackward_kernel(double2* u,double2* v)
{

    int h = blockIdx.x * blockDim.x + threadIdx.x;
    int k = h%NZ;
    int j = h/(NT*NZ);
    int i = (h-j*NT*NZ)/NZ;
  
    if(j<NR/NSTEPS_CONV & i<NT & k<NZ)
    {

        double factor=1.0;

        double2 aux;

        // X indices		
        int ip=i<NT/2 ? i : NTP + (i-NT);

        // Z indices
        int kp=k;

        int hp=j*NTP*NZP+ip*NZP+kp;
        int h =j*NT*NZ+i*NZ+k;
            
        aux=v[hp];
        aux.x*=factor;
        aux.y*=factor;

        if(i==NT/2 || k==NZ-1){

            aux.x=0.0;
            aux.y=0.0;
            
        }
            
        u[h]=aux;
        
    }
}

static __global__ void calcnonlinear_kernel(double* ur,double* ut,double* uz,
                                            double* drur,double* drut,double* druz,
                                            double* dtur,double* dtut,double* dtuz,
                                            double* dzur,double* dzut,double* dzuz,
                                            double* mesh,int i)
{
  
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  if (j<NR/NSTEPS_CONV && k<2*NTP*NZP)
  {
  
  int h=j*(2*NTP*NZP)+k;
    
  double ur_h=ur[h];
  double ut_h=ut[h];
  double uz_h=uz[h];

  double drur_h=drur[h];
  double drut_h=drut[h];
  double druz_h=druz[h];
  
  double dtur_h=dtur[h];
  double dtut_h=dtut[h];
  double dtuz_h=dtuz[h];

  double dzur_h=dzur[h];
  double dzut_h=dzut[h];
  double dzuz_h=dzuz[h];
  
  double rad=mesh[j+i*NR/NSTEPS_CONV];
    
  double udu_r, udu_t, udu_z;
  
  //Normalize
    
  //udu_r(:,i)  = (u_r(:,i)*drdu_r(:,i)  + u_thOr(:,i)*dthdu_r(:,i)&
  //+ u_z(:,i)*dzdu_r(:,i)  - u_th(:,i)*u_thOr(:,i))*cfac(:) 
  
  udu_r = ur_h*drur_h + (ut_h/rad)*dtur_h + uz_h*dzur_h - ut_h*(ut_h/rad);
  
  //udu_th(:,i) = (u_r(:,i)*drdu_th(:,i) + u_thOr(:,i)*dthdu_th(:,i)&
  //+ u_z(:,i)*dzdu_th(:,i) + u_r(:,i)*u_thOr(:,i))*cfac(:)
  
  udu_t = ur_h*drut_h + (ut_h/rad)*dtut_h + uz_h*dzut_h + ur_h*(ut_h/rad);
  
  //udu_z(:,i)  = (u_r(:,i)*drdu_z(:,i)  + u_thOr(:,i)*dthdu_z(:,i)&
  //+ u_z(:,i)*dzdu_z(:,i))*cfac(:)
  
  udu_z = ur_h*druz_h + (ut_h/rad)*dtuz_h + uz_h*dzuz_h;
  
  //Write
  drur[h]=udu_r;
  drut[h]=udu_t;
  druz[h]=udu_z;
  }
}

static dim3 threadsPerBlock;
static dim3 blocksPerGrid;

static vfield u_pad, dur_pad, dut_pad, duz_pad;
static double *mesh_d;
static double *mesh_h;

cublasHandle_t cublasHandle;

void setNonlinear(double* mesh){
 
    //Setting non-linear terms buffers
    
    size_t size_pad=NR/NSTEPS_CONV*NTP*NZP*sizeof(double2);
        
    CHECK_CUDART(cudaMalloc(&u_pad.r,size_pad));
    CHECK_CUDART(cudaMalloc(&u_pad.t,size_pad));
    CHECK_CUDART(cudaMalloc(&u_pad.z,size_pad));

    CHECK_CUDART(cudaMalloc(&dur_pad.r,size_pad));
    CHECK_CUDART(cudaMalloc(&dur_pad.t,size_pad));
    CHECK_CUDART(cudaMalloc(&dur_pad.z,size_pad));
    
    CHECK_CUDART(cudaMalloc(&dut_pad.r,size_pad));
    CHECK_CUDART(cudaMalloc(&dut_pad.t,size_pad));
    CHECK_CUDART(cudaMalloc(&dut_pad.z,size_pad));
    
    CHECK_CUDART(cudaMalloc(&duz_pad.r,size_pad));
    CHECK_CUDART(cudaMalloc(&duz_pad.t,size_pad));
    CHECK_CUDART(cudaMalloc(&duz_pad.z,size_pad));
    
    CHECK_CUDART(cudaMalloc(&mesh_d,NR*sizeof(double)));
    CHECK_CUDART(cudaMemcpy(mesh_d,mesh,NR*sizeof(double),cudaMemcpyHostToDevice)); 

    mesh_h=(double*)malloc(NR*sizeof(double));
    for(int j=0;j<NR;j++){
      mesh_h[j]=mesh[j];
//       printf("\n%e",mesh_h[j]);
    }
    
    //Used to calclate the maximum velocities
    cublasCheck(cublasCreate(&cublasHandle),"Cre");

}

void padForward(double2* aux,double2* u)
{

	dim3 grid,block;
	block.x = 128;
	grid.x = (NR/NSTEPS_CONV*NTP*NZP + block.x - 1)/block.x;

	padForward_kernel<<<grid,block>>>(aux,u);
	
	return;

}

void padBackward(double2* u,double2* aux)
{
	dim3 grid,block;
	block.x = 128;
	grid.x = (NR/NSTEPS_CONV*NT*NZ + block.x - 1)/block.x;

	padBackward_kernel<<<grid,block>>>(u,aux);
	
	return;

}


void calcDut(vfield u, vfield du)
{
    
	dim3 grid,block;
	block.x = 128;
	grid.x = (NR/NSTEPS_CONV*NTP*NZP + block.x - 1)/block.x;
    
    deriv_t_kernel<<<grid,block>>>(u.t,du.t);
    deriv_t_kernel<<<grid,block>>>(u.r,du.r);
    deriv_t_kernel<<<grid,block>>>(u.z,du.z);

    return;
}

void calcDuz(vfield u, vfield du)
{
    
	dim3 grid,block;
	block.x = 128;
	grid.x = (NR/NSTEPS_CONV*NTP*NZP + block.x - 1)/block.x;
    
    deriv_z_kernel<<<grid,block>>>(u.t,du.t);
    deriv_z_kernel<<<grid,block>>>(u.r,du.r);
    deriv_z_kernel<<<grid,block>>>(u.z,du.z);

    return;

}

void calcNonlinear(vfield u,vfield dur,vfield duz,vfield dut,double* mesh_d,int i)
{
      
   dim3 grid, block;
   block.x = 128;
   grid.x = (2*NTP*NZP + block.x - 1)/block.x;
   grid.y = NR/NSTEPS_CONV;
   
  //Calcs ux and uz out of dd_v and w_y 
  calcnonlinear_kernel<<<grid,block>>>
  ((double*)u.r ,(double*)u.t ,(double*)u.z 
  ,(double*)dur.r ,(double*)dur.t ,(double*)dur.z 
  ,(double*)dut.r ,(double*)dut.t ,(double*)dut.z 
  ,(double*)duz.r ,(double*)duz.t ,(double*)duz.z 
  ,mesh_d,i);
    
}

void calcMaxv(vfield u){
  
 
  int size_l=2*NTP*NR*NZP;
  int index;
  
  double Ut_max;
  double Uz_max;
  double Ur_max;
  
  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)u.t,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&Ut_max,(double*)u.t + index-1, sizeof(double), cudaMemcpyDeviceToHost));

  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)u.z,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&Uz_max,(double*)u.z + index-1, sizeof(double), cudaMemcpyDeviceToHost));

  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)u.r,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&Ur_max,(double*)u.r + index-1, sizeof(double), cudaMemcpyDeviceToHost));

  //Calc largest cfl
  Cfl_th=(double)LT/(double)NT/abs(Ut_max);
  Cfl_z =(double)LZ/(double)NZ/abs(Uz_max);
  
  int ir=(index-1)/(2*NTP*NZP);

//   printf("\nMaxVel(t,z,r)=%e,%e,%e",abs(Ut_max),abs(Uz_max),abs(Ur_max));
//   printf("\nIndex_r=%d,%e",ir,mesh_h[ir]);

  double Dr;
  
  if(ir==0){
    Dr=mesh_h[ir] - mesh_h[ir+1];
  }else if(ir==NR-1){
    Dr=mesh_h[ir] - mesh_h[ir-1];
  }else{
    Dr=(mesh_h[ir+1] - mesh_h[ir-1])/2.0;
  }
   
  Cfl_r=abs(Dr)/abs(Uz_max);

//   printf("\nCfl(r,t,z)=%e,%e,%e",Cfl_r,Cfl_th,Cfl_z);

}

//Non-linear terms computed in fractions

void nonLinear(vfield u, vfield dur,int flag,double time){
  
  
  copyVfield(u,dur);
  
  deriv_R_9E(dur.t);
  deriv_R_9E(dur.z);
  deriv_R_9E(dur.r);
 
  //Print torque
  if(flag==1) readTorq(u.t,dur.t,time);

  //Perform the convolution in steps to save memory
  for(int i=0;i<NSTEPS_CONV;i++){
  
    
        size_t offset=i*NR/NSTEPS_CONV*NT*NZ;
        
        padForward(u_pad.r,u.r+offset);
        padForward(u_pad.t,u.t+offset);
        padForward(u_pad.z,u.z+offset);
        

        calcDuz(u_pad,duz_pad);
        calcDut(u_pad,dut_pad); 

        fftBackward(u_pad.t);
        fftBackward(u_pad.z);
        fftBackward(u_pad.r);
        
        calcMaxv(u_pad);
        
        //You can save all these transforms 
        //by deriving in physical space: future work
        padForward(dur_pad.r, dur.r+offset);
        padForward(dur_pad.t, dur.t+offset);
        padForward(dur_pad.z, dur.z+offset);

        fftBackward(dur_pad.r);
        fftBackward(dur_pad.t);
        fftBackward(dur_pad.z);
                    
        fftBackward(duz_pad.r);
        fftBackward(duz_pad.t);
        fftBackward(duz_pad.z);
            
        fftBackward(dut_pad.r);
        fftBackward(dut_pad.t);
        fftBackward(dut_pad.z);

        calcNonlinear(u_pad,dur_pad,duz_pad,dut_pad,mesh_d,i);
        
        fftForward(dur_pad.t);
        fftForward(dur_pad.z);
        fftForward(dur_pad.r);
        
        padBackward(dur.t+offset,dur_pad.t);
        padBackward(dur.z+offset,dur_pad.z);
        padBackward(dur.r+offset,dur_pad.r);
        
        
    }
     
    //Forcing in the z direction
    //forcing(dur.z);
        
}
