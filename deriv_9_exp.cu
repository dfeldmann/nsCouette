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

//Buffers for the coefficients
const double* aa_f,*bb_f,*cc_f,*dd_f,*ee_f,*ff_f,*gg_f,*hh_f,*jj_f;
const double* aa_s,*bb_s,*cc_s,*dd_s,*ee_s,*ff_s,*gg_s,*hh_s,*jj_s;

__global__ void solve_ep_t_9_inplace(double* __restrict x,
              const double* aa,const double*bb,const double*cc,const double*dd,const double*ee,
              const double*ff,const double*gg,const double* hh,const double*jj)
{

  int k = threadIdx.x + blockDim.x*blockIdx.x;
  
  int sizep=2*NT*NZ;
  
  if(k>=sizep) return;

  //Sets the pointers, each thread at a position

  x += k;
  
  //double2 d[NR];
  double xA, xB, xC, xD, xE, xF, xG, xH, xJ;
  double xT;

  //Reads first nine elements
  xA = x[0*sizep];
  xB = x[1*sizep];
  xC = x[2*sizep];
  xD = x[3*sizep];
  xE = x[4*sizep];
  xF = x[5*sizep];
  xG = x[6*sizep];
  xH = x[7*sizep];
  xJ = x[8*sizep];

  //j=0
   xT=ee[0]*xA + ff[0]*xB + gg[0]*xC + hh[0]*xD + jj[0]*xE ;//RHS
   x[0*sizep]=xT;

  //j=1
   xT=dd[1]*xA + ee[1]*xB + ff[1]*xC + gg[1]*xD + hh[1]*xE + jj[1]*xF ;//RHS
   x[1*sizep]=xT;

  //j=2
   xT=cc[2]*xA + dd[2]*xB + ee[2]*xC + ff[2]*xD + gg[2]*xE + hh[2]*xF + jj[2]*xG ;//RHS
   x[2*sizep]=xT;
  
  //j=3
   xT=bb[3]*xA + cc[3]*xB + dd[3]*xC + ee[3]*xD + ff[3]*xE + gg[3]*xF + hh[3]*xG + jj[3]*xH;//RHS
   x[3*sizep]=xT;

   for(int i=4;i<NR;i++){
     
    xT=aa[i]*xA + bb[i]*xB + cc[i]*xC + dd[i]*xD + ee[i]*xE + ff[i]*xF + gg[i]*xG + hh[i]*xH + jj[i]*xJ;//RHS
    x[i*sizep]=xT;
   
    xA = xB;
    xB = xC;
    xC = xD;
    xD = xE;
    xE = xF;
    xF = xG;
    xG = xH;
    xH = xJ;
    
    if(i<NR-5){ //Read next elements
      xJ = x[(i+5)*sizep];
    }
   
   }
}


void fd_weigths(double x0,double* x,int n,double* L1,double*L2){
    
    //Highly non-compact, non-optimised function to calculate the coefficients
    //of the FD with Lagrange pol. Just to make sure I understand.
    //x0, point to evaluate
    //x[n] , points of the function
    //n , length of the stencil
    //L1[n], coefficients of the first deriv
    //L2[n], coefficients of the secod deriv    
    
    //First polinomia
    for(int i=0;i<n;i++){
      L1[i]=0.0;
      for(int j=0;j<n;j++){
        if(i!=j){
          double prod=1.0;
          for(int m=0;m<n;m++){
            if(m!=i && m!=j){
              prod=prod*(x0-x[m])/(x[i]-x[m]);
            }
          }
          L1[i]=L1[i]+prod/(x[i]-x[j]);
        }
      }
    }
  
    //Second polinomia
    for(int i=0;i<n;i++){
      L2[i]=0.0;
      for(int j=0;j<n;j++){
        if(i!=j){
          double sum=0.0;
          for(int m=0;m<n;m++){
            if(m!=i && m!=j){
              double prod=1.0;
              for(int l=0;l<n;l++){
                if(l!=i && l!=j && l!=m){
                  prod=prod*(x0-x[l])/(x[i]-x[l]);
                }  
              }  
            sum=sum+prod/(x[i]-x[m]);
            }
          }
          L2[i]=L2[i]+sum/(x[i]-x[j]);
        }
      }
    }
    return;
} 



void checkCoef(void){
  
  int n=9;

  double* L1=(double*)malloc(sizeof(double)*n);
  double* L2=(double*)malloc(sizeof(double)*n);
  double* x =(double*)malloc(sizeof(double)*n);

  double x0;
 
  for(int i=0;i<n;i++){
  x[i]=(double)(i-4);
  }
  
  x0=0;
  printf("\nChecking coeff");
  fd_weigths(x0,x,n,L1,L2);

//   x0=-4;
//   n=4;
  
  for(int i=0;i<n;i++){
  x[i]=(double)(i);
  }
  x0=0;
  fd_weigths(x0,x,5,L1,L2);
  
  
   for(int i=0;i<5;i++){
     printf("\n%d,%f",i,L1[i]);
   }

  for(int i=0;i<5;i++){
     printf("\n%d,%f",i,L2[i]);
   }
}


static void setCoefficientsDerivatives(double* mesh_r){

  double* coef1=(double*)malloc(9*sizeof(double));
  double* coef2=(double*)malloc(9*sizeof(double));

  //FIRST
  cudaMalloc((void**)&aa_f,NR*sizeof(double));
  cudaMalloc((void**)&bb_f,NR*sizeof(double));
  cudaMalloc((void**)&cc_f,NR*sizeof(double));
  cudaMalloc((void**)&dd_f,NR*sizeof(double)); 
  cudaMalloc((void**)&ee_f,NR*sizeof(double));
  cudaMalloc((void**)&ff_f,NR*sizeof(double));
  cudaMalloc((void**)&gg_f,NR*sizeof(double));
  cudaMalloc((void**)&hh_f,NR*sizeof(double));
  cudaMalloc((void**)&jj_f,NR*sizeof(double));

  //SECOND
  cudaMalloc((void**)&aa_s,NR*sizeof(double));
  cudaMalloc((void**)&bb_s,NR*sizeof(double));
  cudaMalloc((void**)&cc_s,NR*sizeof(double));
  cudaMalloc((void**)&dd_s,NR*sizeof(double)); 
  cudaMalloc((void**)&ee_s,NR*sizeof(double));
  cudaMalloc((void**)&ff_s,NR*sizeof(double));
  cudaMalloc((void**)&gg_s,NR*sizeof(double));
  cudaMalloc((void**)&hh_s,NR*sizeof(double));
  cudaMalloc((void**)&jj_s,NR*sizeof(double));

  double* aa_h=(double*)malloc(NR*sizeof(double));
  double* bb_h=(double*)malloc(NR*sizeof(double));
  double* cc_h=(double*)malloc(NR*sizeof(double));  
  double* dd_h=(double*)malloc(NR*sizeof(double));
  double* ee_h=(double*)malloc(NR*sizeof(double));
  double* ff_h=(double*)malloc(NR*sizeof(double));
  double* gg_h=(double*)malloc(NR*sizeof(double));
  double* hh_h=(double*)malloc(NR*sizeof(double));
  double* jj_h=(double*)malloc(NR*sizeof(double));

  ///////i=0///////
  
  fd_weigths(mesh_r[0],mesh_r,5,coef1,coef2);
  
  //RHS
  aa_h[0]=0.0;
  bb_h[0]=0.0;
  cc_h[0]=0.0;
  dd_h[0]=0.0;
  ee_h[0]=coef1[0];
  ff_h[0]=coef1[1];
  gg_h[0]=coef1[2];
  hh_h[0]=coef1[3];
  jj_h[0]=coef1[4];
  
  /////////////i=1///////////////////
  
  fd_weigths(mesh_r[1],mesh_r,6,coef1,coef2);
  
  //RHS
  aa_h[1]=0.0;
  bb_h[1]=0.0;
  cc_h[1]=0.0;
  dd_h[1]=coef1[0];
  ee_h[1]=coef1[1];
  ff_h[1]=coef1[2];
  gg_h[1]=coef1[3];
  hh_h[1]=coef1[4];
  jj_h[1]=coef1[5];


  /////////////i=2///////////////////

  fd_weigths(mesh_r[2],mesh_r,7,coef1,coef2);
  
  //RHS
  aa_h[2]=0.0;
  bb_h[2]=0.0;
  cc_h[2]=coef1[0];
  dd_h[2]=coef1[1];
  ee_h[2]=coef1[2];
  ff_h[2]=coef1[3];
  gg_h[2]=coef1[4];
  hh_h[2]=coef1[5];
  jj_h[2]=coef1[6];

  /////////i=3/////////////
  
  fd_weigths(mesh_r[3],mesh_r,8,coef1,coef2);
  
  //RHS
  aa_h[3]=0.0;
  bb_h[3]=coef1[0];
  cc_h[3]=coef1[1];
  dd_h[3]=coef1[2];
  ee_h[3]=coef1[3];
  ff_h[3]=coef1[4];
  gg_h[3]=coef1[5];
  hh_h[3]=coef1[6];
  jj_h[3]=coef1[7];
  
  ////////interior grid points////////

  for(int i=4;i<NR-4;i++){

    fd_weigths(mesh_r[i],mesh_r+(i-4),9,coef1,coef2);

    //RHS
    aa_h[i]=coef1[0];
    bb_h[i]=coef1[1];
    cc_h[i]=coef1[2];
    dd_h[i]=coef1[3];
    ee_h[i]=coef1[4];
    ff_h[i]=coef1[5];
    gg_h[i]=coef1[6];
    hh_h[i]=coef1[7];
    jj_h[i]=coef1[8];

  } 


  /////////////i=NR-4///////////////////

  fd_weigths(mesh_r[NR-4],mesh_r+(NR-8),8,coef1,coef2);

  //RHS
  aa_h[NR-4]=coef1[0];
  bb_h[NR-4]=coef1[1];
  cc_h[NR-4]=coef1[2];
  dd_h[NR-4]=coef1[3];
  ee_h[NR-4]=coef1[4];
  ff_h[NR-4]=coef1[5];
  gg_h[NR-4]=coef1[6];
  hh_h[NR-4]=coef1[7];
  jj_h[NR-4]=0;

  /////////////i=NR-3///////////////////
  
  fd_weigths(mesh_r[NR-3],mesh_r+(NR-7),7,coef1,coef2);

  //RHS
  aa_h[NR-3]=coef1[0];
  bb_h[NR-3]=coef1[1];
  cc_h[NR-3]=coef1[2];
  dd_h[NR-3]=coef1[3];
  ee_h[NR-3]=coef1[4];
  ff_h[NR-3]=coef1[5];
  gg_h[NR-3]=coef1[6];
  hh_h[NR-3]=0;
  jj_h[NR-3]=0;  

  /////////////i=NR-2///////////////////
  
  fd_weigths(mesh_r[NR-2],mesh_r+(NR-6),6,coef1,coef2);

  //RHS
  aa_h[NR-2]=coef1[0];
  bb_h[NR-2]=coef1[1];
  cc_h[NR-2]=coef1[2];
  dd_h[NR-2]=coef1[3];
  ee_h[NR-2]=coef1[4];
  ff_h[NR-2]=coef1[5];
  gg_h[NR-2]=0;
  hh_h[NR-2]=0;
  jj_h[NR-2]=0;
  
  /////////////i=NR-1///////////////////
  
  fd_weigths(mesh_r[NR-1],mesh_r+(NR-5),5,coef1,coef2);

  //RHS
  aa_h[NR-1]=coef1[0];
  bb_h[NR-1]=coef1[1];
  cc_h[NR-1]=coef1[2];
  dd_h[NR-1]=coef1[3];
  ee_h[NR-1]=coef1[4];
  ff_h[NR-1]=0;
  gg_h[NR-1]=0;
  hh_h[NR-1]=0;
  jj_h[NR-1]=0;
  


  //Check coeficients
//   FILE* fp=fopen("coefficients_first_derivative.dat","w");
//   for(int j=0;j<NR;j++){
//     fprintf(fp,"\n %03d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e",j,
//     mesh_r[j],aa_h[j],bb_h[j],cc_h[j],dd_h[j],ee_h[j],ff_h[j],gg_h[j],hh_h[j],jj_h[j]);
//   } 
//   fclose(fp);
  


  //Check coeficients
  size_t fsize; 
  FILE* fp=fopen("coefficients_first_derivative_9.bin","wb");
  fsize =fwrite( (unsigned char*)mesh_r,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)aa_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)bb_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)cc_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)dd_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)ee_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)ff_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)gg_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)hh_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)jj_h,sizeof(double),NR,fp);
  fclose(fp);
  

  //FIRST

  cudaMemcpy((void*)aa_f,aa_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)bb_f,bb_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)cc_f,cc_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)dd_f,dd_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)ee_f,ee_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)ff_f,ff_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)gg_f,gg_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)hh_f,hh_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)jj_f,jj_h,NR*sizeof(double),cudaMemcpyHostToDevice);



  ////Second derivative////////////

  ///////i=0///////
  
  fd_weigths(mesh_r[0],mesh_r,5,coef1,coef2);
  
  //RHS
  aa_h[0]=0.0;
  bb_h[0]=0.0;
  cc_h[0]=0.0;
  dd_h[0]=0.0;
  ee_h[0]=coef2[0];
  ff_h[0]=coef2[1];
  gg_h[0]=coef2[2];
  hh_h[0]=coef2[3];
  jj_h[0]=coef2[4];
  
  /////////////i=1///////////////////
  
  fd_weigths(mesh_r[1],mesh_r,6,coef1,coef2);
  
  //RHS
  aa_h[1]=0.0;
  bb_h[1]=0.0;
  cc_h[1]=0.0;
  dd_h[1]=coef2[0];
  ee_h[1]=coef2[1];
  ff_h[1]=coef2[2];
  gg_h[1]=coef2[3];
  hh_h[1]=coef2[4];
  jj_h[1]=coef2[5];


  /////////////i=2///////////////////

  fd_weigths(mesh_r[2],mesh_r,7,coef1,coef2);
  
  //RHS
  aa_h[2]=0.0;
  bb_h[2]=0.0;
  cc_h[2]=coef2[0];
  dd_h[2]=coef2[1];
  ee_h[2]=coef2[2];
  ff_h[2]=coef2[3];
  gg_h[2]=coef2[4];
  hh_h[2]=coef2[5];
  jj_h[2]=coef2[6];

  /////////i=3/////////////
  
  fd_weigths(mesh_r[3],mesh_r,8,coef1,coef2);
  
  //RHS
  aa_h[3]=0.0;
  bb_h[3]=coef2[0];
  cc_h[3]=coef2[1];
  dd_h[3]=coef2[2];
  ee_h[3]=coef2[3];
  ff_h[3]=coef2[4];
  gg_h[3]=coef2[5];
  hh_h[3]=coef2[6];
  jj_h[3]=coef2[7];
  
  ////////interior grid points////////

  for(int i=4;i<NR-4;i++){

    fd_weigths(mesh_r[i],mesh_r+(i-4),9,coef1,coef2);

    //RHS
    aa_h[i]=coef2[0];
    bb_h[i]=coef2[1];
    cc_h[i]=coef2[2];
    dd_h[i]=coef2[3];
    ee_h[i]=coef2[4];
    ff_h[i]=coef2[5];
    gg_h[i]=coef2[6];
    hh_h[i]=coef2[7];
    jj_h[i]=coef2[8];

  } 


  /////////////i=NR-4///////////////////

  fd_weigths(mesh_r[NR-4],mesh_r+(NR-8),8,coef1,coef2);

  //RHS
  aa_h[NR-4]=coef2[0];
  bb_h[NR-4]=coef2[1];
  cc_h[NR-4]=coef2[2];
  dd_h[NR-4]=coef2[3];
  ee_h[NR-4]=coef2[4];
  ff_h[NR-4]=coef2[5];
  gg_h[NR-4]=coef2[6];
  hh_h[NR-4]=coef2[7];
  jj_h[NR-4]=0;

  /////////////i=NR-3///////////////////
  
  fd_weigths(mesh_r[NR-3],mesh_r+(NR-7),7,coef1,coef2);

  //RHS
  aa_h[NR-3]=coef2[0];
  bb_h[NR-3]=coef2[1];
  cc_h[NR-3]=coef2[2];
  dd_h[NR-3]=coef2[3];
  ee_h[NR-3]=coef2[4];
  ff_h[NR-3]=coef2[5];
  gg_h[NR-3]=coef2[6];
  hh_h[NR-3]=0;
  jj_h[NR-3]=0;  

  /////////////i=NR-2///////////////////
  
  fd_weigths(mesh_r[NR-2],mesh_r+(NR-6),6,coef1,coef2);

  //RHS
  aa_h[NR-2]=coef2[0];
  bb_h[NR-2]=coef2[1];
  cc_h[NR-2]=coef2[2];
  dd_h[NR-2]=coef2[3];
  ee_h[NR-2]=coef2[4];
  ff_h[NR-2]=coef2[5];
  gg_h[NR-2]=0;
  hh_h[NR-2]=0;
  jj_h[NR-2]=0;
  
  /////////////i=NR-1///////////////////
  
  fd_weigths(mesh_r[NR-1],mesh_r+(NR-5),5,coef1,coef2);

  //RHS
  aa_h[NR-1]=coef2[0];
  bb_h[NR-1]=coef2[1];
  cc_h[NR-1]=coef2[2];
  dd_h[NR-1]=coef2[3];
  ee_h[NR-1]=coef2[4];
  ff_h[NR-1]=0;
  gg_h[NR-1]=0;
  hh_h[NR-1]=0;
  jj_h[NR-1]=0;
  
//   //Check coeficients
//   FILE* fp1=fopen("coefficients_second_derivative.dat","w");
//   for(int j=0;j<NR;j++){
//     fprintf(fp1,"\n % %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e",j,
//     mesh_r[j],aa_h[j],bb_h[j],cc_h[j],dd_h[j],ee_h[j],ff_h[j],gg_h[j],hh_h[j],jj_h[j]);
//   } 
//   fclose(fp1);
  
  fp=fopen("coefficients_second_derivative_9.bin","wb");
  fsize =fwrite( (unsigned char*)mesh_r,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)aa_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)bb_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)cc_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)dd_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)ee_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)ff_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)gg_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)hh_h,sizeof(double),NR,fp);
  fsize =fwrite( (unsigned char*)jj_h,sizeof(double),NR,fp);
  fclose(fp);
  
  //SECOND

  cudaMemcpy((void*)aa_s,aa_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)bb_s,bb_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)cc_s,cc_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)dd_s,dd_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)ee_s,ee_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)ff_s,ff_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)gg_s,gg_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)hh_s,hh_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy((void*)jj_s,jj_h,NR*sizeof(double),cudaMemcpyHostToDevice);
  
  
  free(aa_h);
  free(cc_h);
  free(bb_h);
  free(dd_h);
  free(ee_h); 
  free(ff_h);
  free(gg_h);
  free(hh_h);
  free(jj_h);

  free(coef1);
  free(coef2);
}


void setDeriv_9_exp(double* mesh,size_p sizes)
{

    setCoefficientsDerivatives(mesh);
  
}


void deriv_R_9E(double2* u)
{

    dim3 grid,block;
    block.x=128;
    grid.x=(2*NZ*NT + block.x - 1)/block.x;
    solve_ep_t_9_inplace<<<grid,block>>>((double*)u,aa_f,bb_f,cc_f,dd_f,ee_f,ff_f,gg_f,hh_f,jj_f);


}

void deriv_RR_9E(double2* u)
{


    dim3 grid,block;
    block.x=128;
    grid.x=(2*NZ*NT + block.x - 1)/block.x;
    solve_ep_t_9_inplace<<<grid,block>>>((double*)u,aa_s,bb_s,cc_s,dd_s,ee_s,ff_s,gg_s,hh_s,jj_s);


}
