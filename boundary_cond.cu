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

static __device__ void evaluateBC(double2* vr, double2* vt, double2* vz, double* BRe, 
                                  double* BIm, const double* rad, double kt, double kz, 
                                  const double* dr1, const double* dr2)
{
  
    double drRe,drIm,urRe,urIm;
    double utRe,utIm,uzRe,uzIm;
    
    int index;

    //Calc derivatives
    for(int i=0;i<2;i++){

      //Inner cilinder
      if(i==0){
        
        drRe=0.0; drIm=0.0;
        index=0;
        
        //Calculate deriv
        for(int j=0;j<4;j++){          
            
            int h=j;
            
            drRe += 0.5*dr1[j]*(vr[h*NT*NZ].x + vt[h*NT*NZ].x);
            drIm += 0.5*dr1[j]*(vr[h*NT*NZ].y + vt[h*NT*NZ].y);
        }  
          
          
      }
      
      //Outer cilinder
      if(i==1){
        
        drRe=0.0; drIm=0.0;
        index=NR-1;

        //Calculate deriv
        for(int j=0;j<4;j++){          
            
            int h=NR-1+(-3+j);
            
            drRe += 0.5*dr2[j]*(vr[h*NT*NZ].x + vt[h*NT*NZ].x);
            drIm += 0.5*dr2[j]*(vr[h*NT*NZ].y + vt[h*NT*NZ].y);
        
        }  
          
          
      }


    double d = 1.0/rad[index];

    urRe =  0.5*(vr[index*NT*NZ].x + vt[index*NT*NZ].x);
    urIm =  0.5*(vr[index*NT*NZ].y + vt[index*NT*NZ].y);
    utRe =  0.5*(vr[index*NT*NZ].y - vt[index*NT*NZ].y);
    utIm = -0.5*(vr[index*NT*NZ].x - vt[index*NT*NZ].x);
    uzRe =  vz[index*NT*NZ].x;
    uzIm =  vz[index*NT*NZ].y;
    
    //evaluate no-slip b.c.   
    BRe[i]   = vr[index*NT*NZ].x;
    BIm[i]   = vr[index*NT*NZ].y;
    BRe[2+i] = vt[index*NT*NZ].x;
    BIm[2+i] = vt[index*NT*NZ].y;
    BRe[4+i] = -uzIm;
    BIm[4+i] =  uzRe;

    BRe[6+i] = urRe*d + drRe - d*kt*utIm - kz*uzIm;
    BIm[6+i] = urIm*d + drIm + d*kt*utRe + kz*uzRe;
    
  }

}

static __global__ void genmatrix_kernel(double* M,
                                        double2* ur_id, double2* ur_od, 
                                        double2* ut_id, double2* ut_od, 
                                        double2* uz_id, double2* uz_od,
                                        double2* dpr_id, double2* dpr_od,
                                        double2* dpt_id, double2* dpt_od,
                                        double2* dpz_id, double2* dpz_od,
                                        double* mesh_d, double* dr1,
                                        double* dr2, double2* zero)                                    
{

  int k = threadIdx.x + blockDim.x*blockIdx.x;

  if(k>=NT*NZ) return;

  M  += k;
  
  ur_id+= k;
  ur_od+= k;
  
  ut_id+= k;
  ut_od+= k;
  
  uz_id+= k;
  uz_od+= k;
  
  dpr_id+= k;
  dpr_od+= k;
  
  dpt_id+= k;
  dpt_od+= k;
  
  dpz_id+= k;
  dpz_od+= k;
  
  int k0=k%NZ;
  int t0=k/NZ;

  // X indices		
  double  kt=t0<NT/2 ? (double)t0 : (double)t0-(double)NT;
  double  kz=k0;
  
  double Bre[MAT], Bim[MAT];
  int n;
  
  //U_PLUS -- INNER
  evaluateBC(ur_id, zero, zero, Bre, Bim, mesh_d, kt, kz, dr1, dr2);

  n=0;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
  
  
  //U_PLUS -- OUTER
  evaluateBC(ur_od, zero, zero, Bre, Bim, mesh_d, kt, kz, dr1, dr2);
  
  n=4;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
  
  //U_MINUS -- INNER
  evaluateBC(zero, ut_id, zero, Bre, Bim, mesh_d, kt, kz, dr1, dr2);
    
  n=1;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
  
  //U_MINUS -- OUTER
  evaluateBC(zero, ut_od, zero, Bre, Bim, mesh_d, kt, kz, dr1, dr2);
  
  n=5;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
  
  //U_Z -- INNER
  evaluateBC(zero, zero, uz_id, Bre, Bim, mesh_d, kt, kz, dr1, dr2);
  
  n=0+2;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
  
  ///U_Z -- OUTTER
  evaluateBC(zero, zero, uz_od, Bre, Bim, mesh_d, kt, kz, dr1, dr2);

  
  n=4+2;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
    
  //PRESSURE  
  
  //Inner cilinder
  evaluateBC(dpr_id, dpt_id, dpz_id, Bre, Bim, mesh_d, kt, kz, dr1, dr2);
  
  n=0+3;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
  
  //Outter cilinder
  evaluateBC(dpr_od, dpt_od, dpz_od, Bre, Bim, mesh_d, kt, kz, dr1, dr2);
  
  n=4+3;
  for(int i=0;i<MAT;i++){
      M[(n*MAT+i)*NT*NZ]=Bre[i];
  }
  
} 

static __global__ void boundary_kernel(double2* ur, double2* ut, double2* uz,
                                       const double* M, 
                                       const double2* hsolur_i, const double2* hsolur_o, 
                                       const double2* hsolut_i, const double2* hsolut_o, 
                                       const double2* hsoluz_i, const double2* hsoluz_o,
                                       const double2* hsolpr_i, const double2* hsolpr_o,
                                       const double2* hsolpt_i, const double2* hsolpt_o,
                                       const double2* hsolpz_i, const double2* hsolpz_o,
                                       const double* mesh_d, const double* dr1, const double* dr2)                                    
{

  int k = threadIdx.x + blockDim.x*blockIdx.x;

  if(k>=NT*NZ) return;

  M  += k;
  
  ur += k;
  ut += k;  
  uz += k;
  
  hsolur_i+= k;
  hsolur_o+= k;
  
  hsolut_i+= k;
  hsolut_o+= k; 
  
  hsoluz_i+= k;
  hsoluz_o+= k;
  
  hsolpr_i+= k;
  hsolpr_o+= k;
  
  hsolpt_i+= k; 
  hsolpt_o+= k;
  
  hsolpz_i+= k; 
  hsolpz_o+= k;
  
  int k0=k%NZ;
  int t0=k/NZ;
  
  double  kt=t0<NT/2 ? (double)t0 : (double)t0-(double)NT;
  double  kz=k0;
  
  double Bre[MAT], Bim[MAT], Ci[MAT], Cr[MAT], Mh[MAT*MAT];
  
  evaluateBC(ur, ut, uz, Bre, Bim, mesh_d, kt, kz, dr1, dr2);

  //Read the matrix inv device
  for(int i=0;i<MAT*MAT;i++){
      Mh[i]=M[i*NT*NZ];
  }
 
  
  //Calculate coefficients
  for(int i=0;i<MAT;i++){
      Cr[i]=0.0;
      Ci[i]=0.0;
  }
  
  for(int j=0;j<MAT;j++){  
    for(int i=0;i<MAT;i++){
      Cr[i]+=-Mh[j*MAT+i]*Bre[j];
      Ci[i]+=-Mh[j*MAT+i]*Bim[j];
    }
  }
  
  double2 ur_i, ut_i, uz_i;
  double2 pr_i, pt_i, pz_i;
  
  double2 ur_o, ut_o, uz_o;
  double2 pr_o, pt_o, pz_o;
  
  double2 ur_d, ut_d, uz_d;
  double  sol[12];
  
  for(int i=0;i<NR;i++){
    
    ur_d = ur[i*NT*NZ];
    ut_d = ut[i*NT*NZ];
    uz_d = uz[i*NT*NZ];
  
    sol[ 0]=hsolur_i[i*NT*NZ].x;
    sol[ 1]=hsolut_i[i*NT*NZ].x;
    sol[ 2]=hsoluz_i[i*NT*NZ].y;
    sol[ 3]=hsolpr_i[i*NT*NZ].x;
    sol[ 4]=hsolpt_i[i*NT*NZ].x;
    sol[ 5]=hsolpz_i[i*NT*NZ].y;

    sol[ 6]=hsolur_o[i*NT*NZ].x;
    sol[ 7]=hsolut_o[i*NT*NZ].x;
    sol[ 8]=hsoluz_o[i*NT*NZ].y;
    sol[ 9]=hsolpr_o[i*NT*NZ].x;
    sol[10]=hsolpt_o[i*NT*NZ].x;
    sol[11]=hsolpz_o[i*NT*NZ].y;
    
    for(int j=0;j<2;j++){

      int m = j*4;
      int p = j*6;
      
      ur_d.x = ur_d.x + Cr[m+0]*sol[p+0]
             + Cr[m+3]*sol[p+3];
      ur_d.y = ur_d.y + Ci[m+0]*sol[p+0]
             + Ci[m+3]*sol[p+3];             
             
      ut_d.x = ut_d.x + Cr[m+1]*sol[p+1]
             + Cr[m+3]*sol[p+4];
      ut_d.y = ut_d.y + Ci[m+1]*sol[p+1]
             + Ci[m+3]*sol[p+4];
             
      uz_d.x = uz_d.x - Ci[m+2]*sol[p+2]
             - Ci[m+3]*sol[p+5];
      uz_d.y = uz_d.y + Cr[m+2]*sol[p+2]
             + Cr[m+3]*sol[p+5];
             
    }
    
    double2 zero={0.0, 0.0};
    
    if((abs(k0)+abs(t0))==0){
    
//       ur[i*NT*NZ] = zero;
    

    }else{
    
      ur[i*NT*NZ] = ur_d;
      ut[i*NT*NZ] = ut_d;
      uz[i*NT*NZ] = uz_d;
    
    }  
    
    }
  
} 


static __global__ void deriv_t_kernel(double2* u,double2* p,double* mesh)
{
  
  int l = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  if (j<NR && l<NT*NZ)
  {
  
  int h=j*(NT*NZ)+l;

  int k=l%NZ;
  int i=l/NZ;
      
  double kz;
  double kt;

  kt=i<NT/2 ? (double)i : (double)i-(double)NT;
  kz=(double)k;

  //Fraction
  kz=(PI2/LZ)*kz;
  kt=(PI2/LT)*kt; 
  
  double rad=mesh[j];
  
  double2 uh=p[h];
  double2 duh;
  
  duh.x=  -kt*uh.y/rad;
  duh.y=   kt*uh.x/rad;
    
  u[h]=duh;

  }
}

static __global__ void deriv_z_kernel(double2* u,double2* p)
{
  
  int l = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  if (j<NR&& l<NT*NZ)
  {
  
  int h=j*(NT*NZ)+l;

  int k=l%NZ;
  int i=l/NZ;

  double kz;
  double kt;

  kt=i<NT/2 ? (double)i : (double)i-(double)NT;
  kz=(double)k;

  //Fraction
  kz=(PI2/LZ)*kz;
  kt=(PI2/LT)*kt; 
  
  double2 uh=p[h];
  double2 duh;

  duh.x=  -kz*uh.y;
  duh.y=   kz*uh.x;

  u[h]=duh;
    
  
  }
}


static double* mesh_d;
static double *M, *invM, *T, *Bvalues, *B;
static double2* coeff;

static vfield aux1, aux2;
static vfield hsolu_i, hsolu_o, hsolp_i, hsolp_o, hsold_i, hsold_o;

static int* pivotArray;
static double2* wbuffer;
static double2* zero;

static double* dr1, *dr2;

void setBoundary(double* mesh){
 
    CHECK_CUDART(cudaMalloc((void**)&mesh_d,sizeof(double)*NR));
    CHECK_CUDART(cudaMemcpy(mesh_d,mesh,sizeof(double)*NR,cudaMemcpyHostToDevice)); 
    
    size_t size=NZ*NT*NR*sizeof(double2);
    
    CHECK_CUDART(cudaMalloc((void**)&aux1.r,size));
    CHECK_CUDART(cudaMalloc((void**)&aux1.t,size));
    CHECK_CUDART(cudaMalloc((void**)&aux1.z,size));
    
    CHECK_CUDART(cudaMalloc((void**)&aux2.r,size));
    CHECK_CUDART(cudaMalloc((void**)&aux2.t,size));
    CHECK_CUDART(cudaMalloc((void**)&aux2.z,size));
    
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
    CHECK_CUDART(cudaMalloc((void**)&M,MAT*MAT*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&T,MAT*MAT*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&invM,MAT*MAT*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&pivotArray,MAT*NT*NZ*sizeof(int)));
    CHECK_CUDART(cudaMalloc((void**)&Bvalues,MAT*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&B,MAT*NT*NZ*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&coeff,MAT*NT*NZ*sizeof(double2)));

    //Coefficients to calculate the first derivative at the walls
    CHECK_CUDART(cudaMalloc((void**)&dr1,5*sizeof(double)));
    CHECK_CUDART(cudaMalloc((void**)&dr2,5*sizeof(double)));

    //Set zero
    CHECK_CUDART(cudaMalloc((void**)&zero,NR*NT*NZ*sizeof(double2)));
    zerosGPU(zero,NR*NT*NZ);

    
    //Calc the weigths
    
    double coef1[4], coef2[4];
    size_t fsize;
    
    fd_weigths(mesh[0], mesh,4,coef1,coef2);
    CHECK_CUDART(cudaMemcpy(dr1,coef1,sizeof(double)*4,cudaMemcpyHostToDevice)); 
    
    FILE * fp;
    
    fp=fopen("dr1.bin","wb");
    fsize =fwrite( (unsigned char*)coef1,sizeof(double),4,fp);
    fclose(fp);
    
    fd_weigths(mesh[NR-1],mesh+(NR-4),4,coef1,coef2);
    CHECK_CUDART(cudaMemcpy(dr2,coef1,sizeof(double)*4,cudaMemcpyHostToDevice)); 
    
    fp=fopen("dr2.bin","wb");
    fsize =fwrite( (unsigned char*)coef1,sizeof(double),4,fp);
    fclose(fp);
    
}

void calcDivp(vfield u){
    
    dim3 grid, block;
    block.x = 128;
    grid.x = (NT*NZ + block.x - 1)/block.x;
    grid.y = NR;
    
    deriv_t_kernel<<<grid,block>>>(u.t,u.r,mesh_d);
    deriv_z_kernel<<<grid,block>>>(u.z,u.r);
    
    return;
}

void calcHsol(double* M, vfield hsolu_i, vfield hsolu_o, vfield hsolp_i, vfield hsolp_o,
              double d_implicit, double dt)
{
    
    double2 bri;
    double2 bro;

    //Uplus
    
    //iner
    bri.x=1.0; bri.y=0.0;
    bro.x=0.0; bro.y=0.0;
    
    implicit_step_hom(hsolu_i.r,wbuffer,(double*)aux1.r, (double*)aux1.t,
                       (double*)aux1.z, d_implicit, dt,0,bri,bro);
    
    //Outter
    bri.x=0.0; bri.y=0.0;
    bro.x=1.0; bro.y=0.0;
    
    implicit_step_hom(hsolu_o.r,wbuffer,(double*)aux1.r, (double*)aux1.t, 
                       (double*)aux1.z, d_implicit, dt,0,bri,bro);

    
    //Uminus
    
    //iner
    bri.x=1.0; bri.y=0.0;
    bro.x=0.0; bro.y=0.0;
    
    implicit_step_hom(hsolu_i.t,wbuffer,(double*)aux1.r, (double*)aux1.t,
                       (double*)aux1.z, d_implicit, dt,1,bri,bro);
    
    
    //Outter
    bri.x=0.0; bri.y=0.0;
    bro.x=1.0; bro.y=0.0;
    
    implicit_step_hom(hsolu_o.t,wbuffer,(double*)aux1.r, (double*)aux1.t, 
                       (double*)aux1.z, d_implicit, dt,1,bri,bro);
    
    
    //u_z
    
    //iner
    bri.x=0.0; bri.y=1.0;
    bro.x=0.0; bro.y=0.0;
    
    implicit_step_hom(hsolu_i.z,wbuffer,(double*)aux1.r, (double*)aux1.t, 
                       (double*)aux1.z, d_implicit, dt,2,bri,bro);

    //Outter
    bri.x=0.0; bri.y=0.0;
    bro.x=0.0; bro.y=1.0;
    
    implicit_step_hom(hsolu_o.z,wbuffer,(double*)aux1.r, (double*)aux1.t, 
                       (double*)aux1.z, d_implicit, dt,2,bri,bro);
    
    //Pressure
    
    //iner
    bri.x=-1.0; bri.y=0.0;
    bro.x= 0.0; bro.y=0.0;
    calcPressure_hom(hsolp_i.r,wbuffer,(double*)aux1.r, (double*)aux1.t, 
                      (double*)aux1.z,bri,bro);
    
    //Calc the divergence of p hsol
    calcDivp(hsolp_i);
    deriv_R_7E(hsolp_i.r);
    decoupleForward(hsolp_i);    

    //outter
    bri.x= 0.0; bri.y=0.0;
    bro.x=-1.0; bro.y=0.0;
    calcPressure_hom(hsolp_o.r,wbuffer,(double*)aux1.r, (double*)aux1.t, 
                      (double*)aux1.z,bri,bro);
    
    calcDivp(hsolp_o);
    deriv_R_7E(hsolp_o.r);
    decoupleForward(hsolp_o);    
        
    dim3 grid, block;
    block.x = 128;
    grid.x = (NT*NZ + block.x - 1)/block.x;

    genmatrix_kernel<<<grid,block>>>(M,hsolu_i.r,hsolu_o.r,hsolu_i.t,hsolu_o.t,
                                       hsolu_i.z,hsolu_o.z,hsolp_i.r,hsolp_o.r,
                                       hsolp_i.t,hsolp_o.t,hsolp_i.z,hsolp_o.z,
                                       mesh_d, dr1, dr2,zero);                            
        
    //Invert infinity matrix
    invert_infmatrix(M);
  
}


void imposeBoundarycond(vfield u, vfield hsolu_i, vfield hsolu_o, vfield hsolp_i,
                        vfield hsolp_o, double* M)
{

    dim3 grid, block;
    block.x = 128;
    grid.x = (NT*NZ + block.x - 1)/block.x;
    //Generates the RHS to obtain the coeff
    boundary_kernel<<<grid,block>>>(u.r,u.t,u.z,(const double*)M,
                                    hsolu_i.r,hsolu_o.r,hsolu_i.t,hsolu_o.t,
                                    hsolu_i.z,hsolu_o.z,hsolp_i.r,hsolp_o.r,
                                    hsolp_i.t,hsolp_o.t,hsolp_i.z,hsolp_o.z,
                                    mesh_d, dr1, dr2); 
}
    
    

