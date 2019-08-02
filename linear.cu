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

#include "TayC.h"

static __global__ void laplacian_kernel(double2* rp, double2* rm, double2* rz,
                                        double2* up, double2* um, double2* uz,
                                        double2* dup,double2* dum, double2* duz,
                                        double2* ddup, double2* ddum, double2* dduz,
                                        double* mesh, double d_implicit, double dt)
{
  
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int h=j*NT*NZ+k;
  
  if (j<NR && k<NT*NZ)
  {
      
        double kz;
        double kt;

        double idt=1.0/dt;

        int i=k/NZ;
        int l=k%NZ;
        
        kt=i<NT/2 ? (double)i : (double)i-(double)NT;
        kz=(double)l;

        //Fraction
        kz=(PI2/LZ)*kz;
        kt=(PI2/LT)*kt; 

        double2 rph=rp[h];
        double2 rmh=rm[h];
        double2 rzh=rz[h];

        double2 uph=up[h];
        double2 umh=um[h];
        double2 uzh=uz[h];

        double2 duph=dup[h];
        double2 dumh=dum[h];
        double2 duzh=duz[h];

        double2 dduph=ddup[h];
        double2 ddumh=ddum[h];
        double2 dduzh=dduz[h];

        double rad=mesh[j];
        
        double2 lap_p, lap_m, lap_z;

        lap_p.x= (1.0-d_implicit)*(dduph.x + duph.x/rad - (kt*kt/(rad*rad) + kz*kz)*uph.x);
        lap_p.y= (1.0-d_implicit)*(dduph.y + duph.y/rad - (kt*kt/(rad*rad) + kz*kz)*uph.y);

        lap_m.x= (1.0-d_implicit)*(ddumh.x + dumh.x/rad - (kt*kt/(rad*rad) + kz*kz)*umh.x);
        lap_m.y= (1.0-d_implicit)*(ddumh.y + dumh.y/rad - (kt*kt/(rad*rad) + kz*kz)*umh.y);

        lap_z.x= (1.0-d_implicit)*(dduzh.x + duzh.x/rad - (kt*kt/(rad*rad) + kz*kz)*uzh.x);
        lap_z.y= (1.0-d_implicit)*(dduzh.y + duzh.y/rad - (kt*kt/(rad*rad) + kz*kz)*uzh.y);

        lap_p.x=lap_p.x + idt*uph.x + (1.0-d_implicit)*(-1.0 - 2.0*kt)/(rad*rad)*uph.x;
        lap_p.y=lap_p.y + idt*uph.y + (1.0-d_implicit)*(-1.0 - 2.0*kt)/(rad*rad)*uph.y;

        lap_m.x=lap_m.x + idt*umh.x + (1.0-d_implicit)*(-1.0 + 2.0*kt)/(rad*rad)*umh.x;
        lap_m.y=lap_m.y + idt*umh.y + (1.0-d_implicit)*(-1.0 + 2.0*kt)/(rad*rad)*umh.y;

        lap_z.x=lap_z.x + idt*uzh.x;
        lap_z.y=lap_z.y + idt*uzh.y;
        
        rph.x=rph.x + lap_p.x;
        rph.y=rph.y + lap_p.y;
        
        rmh.x=rmh.x + lap_m.x;
        rmh.y=rmh.y + lap_m.y;

        rzh.x=rzh.x + lap_z.x;
        rzh.y=rzh.y + lap_z.y;

        rp[h]=rph;
        rm[h]=rmh;
        rz[h]=rzh;

  }

    
} 

static __global__ void addpressure_kernel(double2* rp, double2* rm, double2* rz,
                                          double2* p, double2* dp, double* mesh)                                    
{
  
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int h=j*NT*NZ+k;
  
  if (j<NR && k<NT*NZ)
  {
      
        double kz;
        double kt;

        int i=k/NZ;
        int l=k%NZ;
        
        kt=i<NT/2 ? (double)i : (double)i-(double)NT;
        kz=(double)l;

        //Fraction
        kz=(PI2/LZ)*kz;
        kt=(PI2/LT)*kt; 

        double2 rph=rp[h];
        double2 rmh=rm[h];
        double2 rzh=rz[h];

        double2 ph = p[h];
        double2 dph=dp[h];

        double rad=mesh[j];

        rph.x=rph.x - dph.x + kt*ph.x/rad;
        rph.y=rph.y - dph.y + kt*ph.y/rad;
        
        rmh.x=rmh.x - dph.x - kt*ph.x/rad;
        rmh.y=rmh.y - dph.y - kt*ph.y/rad;

        rzh.x=rzh.x - ( -kz*ph.y);
        rzh.y=rzh.y - (  kz*ph.x);

        rp[h]=rph;
        rm[h]=rmh;
        rz[h]=rzh;

  }

    
} 

static __global__ void calcdiv_kernel(double2* div, double2* ur, double2* ut,
                                      double2* uz, double2* dr, double* mesh)                                    
{
  
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int h=j*NT*NZ+k;
  
  if (j<NR && k<NT*NZ)
  {
      
        double kz;
        double kt;

        int i=k/NZ;
        int l=k%NZ;
        
        kt=i<NT/2 ? (double)i : (double)i-(double)NT;
        kz=(double)l;

        //Fraction
        kz=(PI2/LZ)*kz;
        kt=(PI2/LT)*kt; 

        double2 urh=ur[h];
        double2 uth=ut[h];
        double2 uzh=uz[h];
        double2 drh=dr[h];

        double rad=mesh[j];
        double2 divh;
        
        divh.x=drh.x + urh.x/rad - kz*uzh.y - kt*uth.y/rad;
        divh.y=drh.y + urh.y/rad + kz*uzh.x + kt*uth.x/rad;

        div[h]=divh;

  }

    
} 

static __global__ void forcing_kernel(double2* rz, double* mesh)                                    
{
  
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int h=j*NT*NZ+k;
  
  if (j<NR && k<NT*NZ)
  {
      
        double kz;
        double kt;

        int i=k/NZ;
        int l=k%NZ;
        
        kt=i<NT/2 ? (double)i : (double)i-(double)NT;
        kz=(double)l;

        //Fraction
        kz=(PI2/LZ)*kz;
        kt=(PI2/LT)*kt; 

        double2 rzh=rz[h];

        double rad=mesh[j];

        if(abs(kt)==1 & kz==0){
          rzh.x=rzh.x + 10000*kt;
          rzh.y=rzh.y + 10000*kt;
        }

        rz[h]=rzh;

  }

    
} 

static __global__ void calcvor_kernel(double2* Or, double2* Ot, double2* Oz, double2* ur, double2* ut, double2* uz, double2* dt, double2* dz, double* mesh){
  
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int h=j*NT*NZ+k;
  
  if (j<NR && k<NT*NZ)
  {
      
        double kz;
        double kt;

        int i=k/NZ;
        int l=k%NZ;
        
        kt=i<NT/2 ? (double)i : (double)i-(double)NT;
        kz=(double)l;

        //Fraction
        kz=(PI2/LZ)*kz;
        kt=(PI2/LT)*kt; 

        double2 urh=ur[h];
        double2 uth=ut[h];
        double2 uzh=uz[h];
        double2 dth=dt[h];
        double2 dzh=dz[h];

        double rad=mesh[j];
        
        double2 omer, omet, omez;
      
        
        omer.x= - kt*uzh.y/rad -(- kz*uth.y);
        omer.y=   kt*uzh.x/rad -(  kz*uth.x);

        omet.x= - kz*urh.y - dzh.x;
        omet.y=   kz*urh.x - dzh.y;
     
        omez.x= dth.x + uth.x/rad  -( -kt*urh.y/rad);
        omez.y= dth.y + uth.y/rad  -(  kt*urh.x/rad);
        
        Or[h]=omer;
        Ot[h]=omet;
        Oz[h]=omez;

  }

    
} 


static double* mesh_d;
static double* M, invM, T, Bvalues;

static vfield aux1, aux2, aux3;
static vfield hsolu_i, hsolu_o, hsolp_i, hsolp_o, hsold_i, hsold_o;

static int* pivotArray;
static double2* wbuffer;



void setLinear(double* mesh){
 
    CHECK_CUDART(cudaMalloc((void**)&mesh_d,sizeof(double)*NR));
    CHECK_CUDART(cudaMemcpy(mesh_d,mesh,sizeof(double)*NR,cudaMemcpyHostToDevice)); 
    
    size_t size  =NZ*NT*NR*sizeof(double2);
    size_t size_d=NZ*NT*NR*sizeof(double2);

    CHECK_CUDART(cudaMalloc((void**)&aux1.r,size));
    CHECK_CUDART(cudaMalloc((void**)&aux1.t,size));
    CHECK_CUDART(cudaMalloc((void**)&aux1.z,size));
    
    CHECK_CUDART(cudaMalloc((void**)&aux2.r,size_d));
    CHECK_CUDART(cudaMalloc((void**)&aux2.t,size_d));
    CHECK_CUDART(cudaMalloc((void**)&aux2.z,size_d));
    
    CHECK_CUDART(cudaMalloc((void**)&aux3.r,size_d));
    CHECK_CUDART(cudaMalloc((void**)&aux3.t,size_d));
    CHECK_CUDART(cudaMalloc((void**)&aux3.z,size_d));    
    
    CHECK_CUDART(cudaMalloc((void**)&hsolu_i.r,size));
    CHECK_CUDART(cudaMalloc((void**)&hsolu_i.t,size));
    CHECK_CUDART(cudaMalloc((void**)&hsolu_i.z,size));
    
    CHECK_CUDART(cudaMalloc((void**)&hsolu_o.r,size));
    CHECK_CUDART(cudaMalloc((void**)&hsolu_o.t,size));
    CHECK_CUDART(cudaMalloc((void**)&hsolu_o.z,size));
    
    CHECK_CUDART(cudaMalloc((void**)&hsolp_i.r,size));
    CHECK_CUDART(cudaMalloc((void**)&hsolp_i.t,size));
    CHECK_CUDART(cudaMalloc((void**)&hsolp_i.z,size));
    
    CHECK_CUDART(cudaMalloc((void**)&hsold_o.r,size));
    CHECK_CUDART(cudaMalloc((void**)&hsold_o.t,size));
    CHECK_CUDART(cudaMalloc((void**)&hsold_i.r,size));
    CHECK_CUDART(cudaMalloc((void**)&hsold_i.t,size));
    
    CHECK_CUDART(cudaMalloc((void**)&wbuffer,size));
    
    //Allocate matrices for inversion
    CHECK_CUDART(cudaMalloc((void**)&M,8*8*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&T,8*8*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&invM,8*8*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&pivotArray,8*NT*NZ*sizeof(int)));
    CHECK_CUDART(cudaMalloc((void**)&Bvalues,8*NT*NZ*sizeof(double)));

    
}



void calcLaplacian(vfield rhs, vfield u, double d_implicit, double dt){
   
   //Calc second and first derivative
   
   copyVfield(u,aux1);
   
   deriv_R_9E(aux1.r);
   deriv_R_9E(aux1.t);
   deriv_R_9E(aux1.z);
  
   
   copyVfield(u,aux2);

   deriv_RR_9E(aux2.r);
   deriv_RR_9E(aux2.t);
   deriv_RR_9E(aux2.z);
   
    
   dim3 grid, block;
   block.x = 128;
   grid.x = (NT*NZ + block.x - 1)/block.x;
   grid.y = NR;
   
   laplacian_kernel<<<grid,block>>>
   (rhs.r,rhs.t,rhs.z,u.r ,u.t ,u.z,
   aux1.r,aux1.t,aux1.z,aux2.r,aux2.t,aux2.z,
   mesh_d, d_implicit, dt);
    
    
}

void pressureProject(vfield rhs){
    
    //Calc pressure
    copyVfield(rhs, aux3); 
    decoupleBackward(aux3);
  
    calcPressure(aux1.r, aux3 , aux1.t, aux1.z, (double*)aux2.r, (double*)aux2.t, 
                 (double*)aux2.z); 
    
    //Calc dp/dr
    copyBuffer(aux1.r, aux2.r);
    deriv_R_7E(aux2.r);
    
    //Add pressure to the RHS
    dim3 grid, block;
    block.x = 128;
    grid.x = (NT*NZ + block.x - 1)/block.x;
    grid.y = NR;
    
    addpressure_kernel<<<grid,block>>>(rhs.r, rhs.t, rhs.z, aux1.r, aux2.r, mesh_d);
        
}

void implictStep(vfield u, vfield rhs, double d_implicit, double dt){
    
    implicit_step(rhs.r,aux1.r, (double*)aux1.t, (double*)aux1.z, (double*)aux2.r, d_implicit, dt, 0);
    implicit_step(rhs.t,aux1.r, (double*)aux1.t, (double*)aux1.z, (double*)aux2.r, d_implicit, dt, 1);
    implicit_step(rhs.z,aux1.r, (double*)aux1.t, (double*)aux1.z, (double*)aux2.r, d_implicit, dt, 2);

    copyVfield(rhs,u);
}

void calcDivergence(double2* div, vfield u){
  
    copyBuffer(u.r, aux1.r);
    deriv_R_7E(aux1.r);
    
    
    //Add pressure to the RHS
    dim3 grid, block;
    block.x = 128;
    grid.x = (NT*NZ + block.x - 1)/block.x;
    grid.y = NR;
    
    calcdiv_kernel<<<grid,block>>>(div,u.r,u.t,u.z,aux1.r, mesh_d);
  
}

void calcVorticity(vfield vor, vfield u){
    
    copyVfield(u, aux1);
    deriv_R_9E(aux1.t);
    deriv_R_9E(aux1.z);

    
    //Add pressure to the RHS
    dim3 grid, block;
    block.x = 128;
    grid.x = (NT*NZ + block.x - 1)/block.x;
    grid.y = NR;
    calcvor_kernel<<<grid,block>>>(vor.r,vor.t,vor.z,u.r,u.t,u.z,aux1.t,aux1.z,mesh_d);

}

void forcing(double2* uz)
{
  
    //Add forcing to the RHS
    dim3 grid, block;
    block.x = 128;
    grid.x = (NT*NZ + block.x - 1)/block.x;
    grid.y = NR;
    forcing_kernel<<<grid,block>>>(uz,mesh_d);
  
  
  
}
