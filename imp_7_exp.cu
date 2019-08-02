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

const double* a1,*b1,*c1,*d1,*e1,*f1,*g1;
const double* a2,*b2,*c2,*d2,*e2,*f2,*g2;

double* mesh_d;

__global__ void solve_ep_t_7_implicit(double2* __restrict x,double2* __restrict y, double* __restrict delta, double* __restrict epsil,
                                      double* __restrict dseta,const double* a1,const double*b1,const double*c1,const double*d1,const double*e1,
                                      const double*f1,const double* g1,const double* a2,const double*b2,const double*c2,const double*d2,
                                      const double*e2,const double*f2,const double* g2, const double* rad,double d_implicit,double dt,int flag,
                                      double Bci, double Bco)
{

  // k=0...NX*NZ

  int k = threadIdx.x + blockDim.x*blockIdx.x;

  if(k>=NT*NZ) return;

  //Sets the pointers, each thread at a position

  double idt=1.0/dt;
  
  x += k;
  y += k;  

  delta += k; 
  epsil += k;
  dseta += k;

  int k0=k%NZ;
  int t0=k/NZ;

  // X indices		
  double  kt=t0<NT/2 ? (double)t0 : (double)t0-(double)NT;
  double  kz=k0;


  //Fraction
  kt=(PI2/LT)*kt;
  kz=(PI2/LZ)*kz;	

  double kk,ri2,diag;

  //double2 d[NY];
  double2 xA, xB, xC, xD, xE, xF, xG;
  double2 xT;

  //Reads first five elements

  //xA.x = (double)x[0*NT*NZ].x;
  //xA.y = (double)x[0*NT*NZ].y;
  xB.x = (double)x[1*NT*NZ].x;
  xB.y = (double)x[1*NT*NZ].y;
  xC.x = (double)x[2*NT*NZ].x;
  xC.y = (double)x[2*NT*NZ].y;
  xD.x = (double)x[3*NT*NZ].x;
  xD.y = (double)x[3*NT*NZ].y;
  xE.x = (double)x[4*NT*NZ].x;
  xE.y = (double)x[4*NT*NZ].y;
  xF.x = (double)x[5*NT*NZ].x;
  xF.y = (double)x[5*NT*NZ].y;
  xG.x = (double)x[6*NT*NZ].x;
  xG.y = (double)x[6*NT*NZ].y;



  //Local variables
  double delta_3,delta_1,delta_2,delta_4;
  double epsil_3,epsil_1,epsil_2,epsil_4;
  double dseta_3,dseta_1,dseta_2,dseta_4;
  
  double alpha_0;
  double gamma_0;
  double betha_0;
  
  double g_1,g_2,g_3;
  double2 y_1,y_2,y_3,y_4;
  double2 x_0,x_1,x_2,x_3;

  double a0,b0,c0,d0,e0,f0,g0;  
       
  //j=0

  d0=1.0;
  e0=0.0;
  f0=0.0;
  g0=0.0;

  delta_1=d0; //D[0]
  epsil_1=e0;
  dseta_1=f0;
  

  g_1=g0;
  
  
  delta[0*NT*NZ]=delta_1;
  epsil[0*NT*NZ]=epsil_1;
  dseta[0*NT*NZ]=dseta_1;

  //Boundary conditions
  
  xT.x=0.0;
  xT.y=0.0; 

  //Boundary conditions for mode 0
  if((abs(k0)+abs(t0))==0){
    xT.y=Bci;
  }
    
  y_1.x=xT.x;
  y_1.y=xT.y;

  y[0*NT*NZ]=y_1;


  //j=1

  ri2=rad[1]*rad[1]; 
  kk=kt*kt/ri2 + kz*kz;
  
  if(flag==0) diag=-1.0/ri2 - 2.0*kt/ri2;
  if(flag==1) diag=-1.0/ri2 + 2.0*kt/ri2;
  if(flag==2) diag=0.0;   
      
  c0=-d_implicit*(c2[1]+1.0/rad[1]*c1[1]);
  d0= idt - d_implicit*(d2[1]+1.0/rad[1]*d1[1]-kk+diag);
  e0=-d_implicit*(e2[1]+1.0/rad[1]*e1[1]);
  f0=-d_implicit*(f2[1]+1.0/rad[1]*f1[1]);
  g0=-d_implicit*(g2[1]+1.0/rad[1]*g1[1]);

  gamma_0=(c0)/delta_1;

  delta_2=d0-epsil_1*gamma_0;//D[1]
  epsil_2=e0-dseta_1*gamma_0;
  dseta_2=f0-g_1*gamma_0;
  g_2=g0;

//   if(k0==1 & t0==1) printf("\ndelta = %d %e",1,delta_2);

  
  delta[1*NT*NZ]=delta_2;
  epsil[1*NT*NZ]=epsil_2;
  dseta[1*NT*NZ]=dseta_2;

// /  if(k0==1 & t0==1) printf("\nepsil_2 %e",delta_2);
  
  xT.x=xB.x;
  xT.y=xB.y;

  y_2.x=xT.x-gamma_0*y_1.x;
  y_2.y=xT.y-gamma_0*y_1.y;

  y[1*NT*NZ]=y_2;

  //j=2
  
  ri2=rad[2]*rad[2]; 
  kk=kt*kt/ri2 + kz*kz;
  
  if(flag==0) diag=-1.0/ri2 - 2.0*kt/ri2;
  if(flag==1) diag=-1.0/ri2 + 2.0*kt/ri2;
  if(flag==2) diag=0.0; 
  
  b0=-d_implicit*(b2[2]+1.0/rad[2]*b1[2]);
  c0=-d_implicit*(c2[2]+1.0/rad[2]*c1[2]);
  d0= idt-d_implicit*(d2[2]+1.0/rad[2]*d1[2]-kk+diag);
  e0=-d_implicit*(e2[2]+1.0/rad[2]*e1[2]);
  f0=-d_implicit*(f2[2]+1.0/rad[2]*f1[2]);
  g0=-d_implicit*(g2[2]+1.0/rad[2]*g1[2]);

	betha_0=(b0)/delta_1;
	gamma_0=(c0-epsil_1*betha_0)/delta_2;

	delta_3=d0-epsil_2*gamma_0-dseta_1*betha_0;//D[2]
	epsil_3=e0-dseta_2*gamma_0-g_1*betha_0;
	dseta_3=f0-g_2*gamma_0;
	g_3=g0;   

//     if(k0==1 & t0==1) printf("\ndelta = %d %e",2,delta_3);


	delta[2*NT*NZ]=delta_3;
	epsil[2*NT*NZ]=epsil_3;
	dseta[2*NT*NZ]=dseta_3;

	xT.x= xC.x;
	xT.y= xC.y;

	y_3.x=xT.x-betha_0*y_1.x-gamma_0*y_2.x;
	y_3.y=xT.y-betha_0*y_1.y-gamma_0*y_2.y;
	y[2*NT*NZ]=y_3;

	//Loop j=3 to j=N-4

    
  for(int j=3; j<NR; j++){
      
    ri2=rad[j]*rad[j]; 
    kk=kt*kt/ri2 + kz*kz;
    
    if(flag==0) diag=-1.0/ri2 - 2.0*kt/ri2;
    if(flag==1) diag=-1.0/ri2 + 2.0*kt/ri2;
    if(flag==2) diag=0.0; 
    
    a0=-d_implicit*(a2[j]+1.0/rad[j]*a1[j]);
    b0=-d_implicit*(b2[j]+1.0/rad[j]*b1[j]);
    c0=-d_implicit*(c2[j]+1.0/rad[j]*c1[j]);
    d0= idt-d_implicit*(d2[j]+1.0/rad[j]*d1[j]-kk+diag);
    e0=-d_implicit*(e2[j]+1.0/rad[j]*e1[j]);
    f0=-d_implicit*(f2[j]+1.0/rad[j]*f1[j]);
    g0=-d_implicit*(g2[j]+1.0/rad[j]*g1[j]);

//     if(k0==1 & t0==1) printf("\nd0 = %d %e",j,d0);

    
    //Boundary conditions

    if(j==NR-1){

      a0=0.0;
      b0=0.0;
      c0=0.0;
      d0=1.0;
      e0=0.0;
      f0=0.0;
      g0=0.0;

    }

    alpha_0= a0/delta_1;   
    betha_0=(b0-epsil_1*alpha_0)/delta_2;
    gamma_0=(c0-epsil_2*betha_0-dseta_1*alpha_0)/delta_3;
        

    delta_4 =d0-epsil_3*gamma_0-dseta_2*betha_0-g_1*alpha_0;
    epsil_4 =e0-dseta_3*gamma_0-g_2*betha_0;
    dseta_4 =f0-g_3*gamma_0;

//     if(k0==1 & t0==1) printf("\ndelta = %d %e",j,delta_4);

    
    //write
    delta[j*NT*NZ]=delta_4;
    delta_1=delta_2;
    delta_2=delta_3;
    delta_3=delta_4;

    //write
    epsil[j*NT*NZ]=epsil_4;
    epsil_1 =epsil_2;
    epsil_2 =epsil_3;
    epsil_3 =epsil_4;


    //write
    dseta[j*NT*NZ]=dseta_4;
    dseta_1 =dseta_2;
    dseta_2 =dseta_3;
    dseta_3 =dseta_4;

    //Read g
    g_1=g_2;
    g_2=g_3;
    g_3=g0;

    xT.x=xD.x;
    xT.y=xD.y;

    //Boundary conditions

    if(j==NR-1){

      xT.x=0.0;
      xT.y=0.0;

      //Boundary conditions for mode 0
      if((abs(k0)+abs(t0))==0){
        xT.y=Bco;
      } 
      
    }

    y_4.x=xT.x-alpha_0*y_1.x-betha_0*y_2.x-gamma_0*y_3.x;
    y_4.y=xT.y-alpha_0*y_1.y-betha_0*y_2.y-gamma_0*y_3.y;

    y[j*NT*NZ]=y_4;
    y_1=y_2;
    y_2=y_3;
    y_3=y_4;

    xA = xB;
    xB = xC;
    xC = xD;
    xD = xE;
    xE = xF;
    xF = xG;

    if(j<NR-4){ //Read next elements
      double2 next_x = x[(j+4)*NT*NZ]; //Read H
      xG.x = (double)next_x.x;
      xG.y = (double)next_x.y;   
    }
  }
 

	//BACKWARD SUSTITUTION

	//N
    x_3.x=y_3.x/delta_3;
    x_3.y=y_3.y/delta_3;

    x[(NR-1)*NT*NZ].x=(double)x_3.x;
    x[(NR-1)*NT*NZ].y=(double)x_3.y;

    
    //N-1
    x_2.x=(y_2.x-epsil_2*x_3.x)/delta_2;
    x_2.y=(y_2.y-epsil_2*x_3.y)/delta_2;

    x[(NR-2)*NT*NZ].x=(double)x_2.x;    
    x[(NR-2)*NT*NZ].y=(double)x_2.y;    

    //N-3
    x_1.x=(y_1.x-epsil_1*x_2.x-dseta_1*x_3.x)/delta_1;
    x_1.y=(y_1.y-epsil_1*x_2.y-dseta_1*x_3.y)/delta_1;

    x[(NR-3)*NT*NZ].x=(double)x_1.x;    
    x[(NR-3)*NT*NZ].y=(double)x_1.y;    

 	for(int j=NR-4;j>=0;j--){

        g0=-d_implicit*(g2[j]+1.0/rad[j]*g1[j]);      

        if(j==0) g0=0.0;

        x_0.x =(y[j*NT*NZ].x-epsil[j*NT*NZ]*x_1.x
        -dseta[j*NT*NZ]*x_2.x-g0*x_3.x)/delta[j*NT*NZ];

        x_0.y =(y[j*NT*NZ].y-epsil[j*NT*NZ]*x_1.y
        -dseta[j*NT*NZ]*x_2.y-g0*x_3.y)/delta[j*NT*NZ];
        //Write
        x[j*NT*NZ].x=(double)x_0.x;
        x[j*NT*NZ].y=(double)x_0.y;

        x_3 =x_2;
        x_2 =x_1;
        x_1 =x_0;
    
    }  


}

__global__ void solve_ep_t_7_implicit_hom(double2* __restrict x,double2* __restrict y, double* __restrict delta, double* __restrict epsil, double* __restrict dseta,const double* a1,const double*b1,const double*c1,const double*d1,const double*e1,const double*f1,const double* g1,const double* a2,const double*b2,const double*c2,const double*d2,const double*e2,const double*f2,const double* g2, const double* rad,double d_implicit,double dt,int flag, double2 Bci, double2 Bco)
{

  // k=0...NX*NZ

  int k = threadIdx.x + blockDim.x*blockIdx.x;

  if(k>=NT*NZ) return;

  //Sets the pointers, each thread at a position

  double idt=1.0/dt;
  
  x += k;
  y += k;  

  delta += k; 
  epsil += k;
  dseta += k;

  int k0=k%NZ;
  int t0=k/NZ;

  // X indices		
  double  kt=t0<NT/2 ? (double)t0 : (double)t0-(double)NT;
  double  kz=k0;


  //Fraction
  kt=(PI2/LT)*kt;
  kz=(PI2/LZ)*kz;	

  double kk,ri2,diag;

  //double2 d[NY];
  double2 xT;

  //Local variables
  double delta_3,delta_1,delta_2,delta_4;
  double epsil_3,epsil_1,epsil_2,epsil_4;
  double dseta_3,dseta_1,dseta_2,dseta_4;
  
  double alpha_0;
  double gamma_0;
  double betha_0;
  
  double g_1,g_2,g_3;
  double2 y_1,y_2,y_3,y_4;
  double2 x_0,x_1,x_2,x_3;

  double a0,b0,c0,d0,e0,f0,g0;  
       
  //j=0

  d0=1.0;
  e0=0.0;
  f0=0.0;
  g0=0.0;

  delta_1=d0; //D[0]
  epsil_1=e0;
  dseta_1=f0;
  

  g_1=g0;
  
  
  delta[0*NT*NZ]=delta_1;
  epsil[0*NT*NZ]=epsil_1;
  dseta[0*NT*NZ]=dseta_1;

  //Boundary conditions
  
  xT.x=Bci.x;
  xT.y=Bci.y; 

  y_1.x=xT.x;
  y_1.y=xT.y;

  y[0*NT*NZ]=y_1;


  //j=1

  ri2=rad[1]*rad[1]; 
  kk=kt*kt/ri2 + kz*kz;
  
  if(flag==0) diag=-1.0/ri2 - 2.0*kt/ri2;
  if(flag==1) diag=-1.0/ri2 + 2.0*kt/ri2;
  if(flag==2) diag=0.0;   
      
  c0=-d_implicit*(c2[1]+1.0/rad[1]*c1[1]);
  d0= idt - d_implicit*(d2[1]+1.0/rad[1]*d1[1]-kk+diag);
  e0=-d_implicit*(e2[1]+1.0/rad[1]*e1[1]);
  f0=-d_implicit*(f2[1]+1.0/rad[1]*f1[1]);
  g0=-d_implicit*(g2[1]+1.0/rad[1]*g1[1]);

  gamma_0=(c0)/delta_1;

  delta_2=d0-epsil_1*gamma_0;//D[1]
  epsil_2=e0-dseta_1*gamma_0;
  dseta_2=f0-g_1*gamma_0;
  g_2=g0;

//   if(k0==1 & t0==1) printf("\ndelta = %d %e",1,delta_2);

  
  delta[1*NT*NZ]=delta_2;
  epsil[1*NT*NZ]=epsil_2;
  dseta[1*NT*NZ]=dseta_2;

// /  if(k0==1 & t0==1) printf("\nepsil_2 %e",delta_2);
  
  xT.x=0.0;
  xT.y=0.0;

  y_2.x=xT.x-gamma_0*y_1.x;
  y_2.y=xT.y-gamma_0*y_1.y;

  y[1*NT*NZ]=y_2;

  //j=2
  
  ri2=rad[2]*rad[2]; 
  kk=kt*kt/ri2 + kz*kz;
  
  if(flag==0) diag=-1.0/ri2 - 2.0*kt/ri2;
  if(flag==1) diag=-1.0/ri2 + 2.0*kt/ri2;
  if(flag==2) diag=0.0; 
  
  b0=-d_implicit*(b2[2]+1.0/rad[2]*b1[2]);
  c0=-d_implicit*(c2[2]+1.0/rad[2]*c1[2]);
  d0= idt-d_implicit*(d2[2]+1.0/rad[2]*d1[2]-kk+diag);
  e0=-d_implicit*(e2[2]+1.0/rad[2]*e1[2]);
  f0=-d_implicit*(f2[2]+1.0/rad[2]*f1[2]);
  g0=-d_implicit*(g2[2]+1.0/rad[2]*g1[2]);

    betha_0=(b0)/delta_1;
    gamma_0=(c0-epsil_1*betha_0)/delta_2;

    delta_3=d0-epsil_2*gamma_0-dseta_1*betha_0;//D[2]
    epsil_3=e0-dseta_2*gamma_0-g_1*betha_0;
    dseta_3=f0-g_2*gamma_0;
    g_3=g0;   

    //     if(k0==1 & t0==1) printf("\ndelta = %d %e",2,delta_3);


    delta[2*NT*NZ]=delta_3;
    epsil[2*NT*NZ]=epsil_3;
    dseta[2*NT*NZ]=dseta_3;

    xT.x= 0.0;
    xT.y= 0.0;

    y_3.x=xT.x-betha_0*y_1.x-gamma_0*y_2.x;
    y_3.y=xT.y-betha_0*y_1.y-gamma_0*y_2.y;
    y[2*NT*NZ]=y_3;

    //Loop j=3 to j=N-4

    
  for(int j=3; j<NR; j++){
      
    ri2=rad[j]*rad[j]; 
    kk=kt*kt/ri2 + kz*kz;
    
    if(flag==0) diag=-1.0/ri2 - 2.0*kt/ri2;
    if(flag==1) diag=-1.0/ri2 + 2.0*kt/ri2;
    if(flag==2) diag=0.0; 
    
    a0=-d_implicit*(a2[j]+1.0/rad[j]*a1[j]);
    b0=-d_implicit*(b2[j]+1.0/rad[j]*b1[j]);
    c0=-d_implicit*(c2[j]+1.0/rad[j]*c1[j]);
    d0= idt-d_implicit*(d2[j]+1.0/rad[j]*d1[j]-kk+diag);
    e0=-d_implicit*(e2[j]+1.0/rad[j]*e1[j]);
    f0=-d_implicit*(f2[j]+1.0/rad[j]*f1[j]);
    g0=-d_implicit*(g2[j]+1.0/rad[j]*g1[j]);

//     if(k0==1 & t0==1) printf("\nd0 = %d %e",j,d0);

    
    //Boundary conditions

    if(j==NR-1){

      a0=0.0;
      b0=0.0;
      c0=0.0;
      d0=1.0;
      e0=0.0;
      f0=0.0;
      g0=0.0;

    }

    alpha_0= a0/delta_1;   
    betha_0=(b0-epsil_1*alpha_0)/delta_2;
    gamma_0=(c0-epsil_2*betha_0-dseta_1*alpha_0)/delta_3;
        

    delta_4 =d0-epsil_3*gamma_0-dseta_2*betha_0-g_1*alpha_0;
    epsil_4 =e0-dseta_3*gamma_0-g_2*betha_0;
    dseta_4 =f0-g_3*gamma_0;

//     if(k0==1 & t0==1) printf("\ndelta = %d %e",j,delta_4);

    
    //write
    delta[j*NT*NZ]=delta_4;
    delta_1=delta_2;
    delta_2=delta_3;
    delta_3=delta_4;

    //write
    epsil[j*NT*NZ]=epsil_4;
    epsil_1 =epsil_2;
    epsil_2 =epsil_3;
    epsil_3 =epsil_4;


    //write
    dseta[j*NT*NZ]=dseta_4;
    dseta_1 =dseta_2;
    dseta_2 =dseta_3;
    dseta_3 =dseta_4;

    //Read g
    g_1=g_2;
    g_2=g_3;
    g_3=g0;

    xT.x=0.0;
    xT.y=0.0;

    //Boundary conditions

    if(j==NR-1){

      xT.x=Bco.x;
      xT.y=Bco.y;

    }

    y_4.x=xT.x-alpha_0*y_1.x-betha_0*y_2.x-gamma_0*y_3.x;
    y_4.y=xT.y-alpha_0*y_1.y-betha_0*y_2.y-gamma_0*y_3.y;

    y[j*NT*NZ]=y_4;
    y_1=y_2;
    y_2=y_3;
    y_3=y_4;
    
  }
 

    //BACKWARD SUSTITUTION

    //N
    x_3.x=y_3.x/delta_3;
    x_3.y=y_3.y/delta_3;

    x[(NR-1)*NT*NZ].x=(double)x_3.x;
    x[(NR-1)*NT*NZ].y=(double)x_3.y;

//     if(k==30)printf("\n%d %e",NR,x_0.x);

    //N-1
    x_2.x=(y_2.x-epsil_2*x_3.x)/delta_2;
    x_2.y=(y_2.y-epsil_2*x_3.y)/delta_2;

    x[(NR-2)*NT*NZ].x=(double)x_2.x;    
    x[(NR-2)*NT*NZ].y=(double)x_2.y;    

    //N-3
    x_1.x=(y_1.x-epsil_1*x_2.x-dseta_1*x_3.x)/delta_1;
    x_1.y=(y_1.y-epsil_1*x_2.y-dseta_1*x_3.y)/delta_1;

    x[(NR-3)*NT*NZ].x=(double)x_1.x;    
    x[(NR-3)*NT*NZ].y=(double)x_1.y;    

 	for(int j=NR-4;j>=0;j--){

        g0=-d_implicit*(g2[j]+1.0/rad[j]*g1[j]);      

        if(j==0) g0=0.0;

        x_0.x =(y[j*NT*NZ].x-epsil[j*NT*NZ]*x_1.x
        -dseta[j*NT*NZ]*x_2.x-g0*x_3.x)/delta[j*NT*NZ];

        x_0.y =(y[j*NT*NZ].y-epsil[j*NT*NZ]*x_1.y
        -dseta[j*NT*NZ]*x_2.y-g0*x_3.y)/delta[j*NT*NZ];
        //Write
        x[j*NT*NZ].x=(double)x_0.x;
        x[j*NT*NZ].y=(double)x_0.y;

        x_3 =x_2;
        x_2 =x_1;
        x_1 =x_0;
    
//         if(k==30 && j==0)printf("\n%d %e",j,x_0.x);
    }  


}


__global__ void solve_ep_t_7_pressure(double2* __restrict x,double2* __restrict y, double2* __restrict ur, double2* __restrict ut, double2* __restrict uz
                                      ,double2* __restrict du, double* __restrict delta, double* __restrict epsil, double* __restrict dseta
                                      ,const double* a1,const double*b1,const double*c1,const double*d1,const double*e1,const double*f1,const double* g1
                                      ,const double* a2,const double*b2,const double*c2,const double*d2,const double*e2,const double*f2,const double* g2
                                      ,const double* rad)
{

  // k=0...NX*NZ

  int k = threadIdx.x + blockDim.x*blockIdx.x;

  if(k>=NT*NZ) return;

  //Sets the pointers, each thread at a position

  x += k;
  y += k;  

  delta += k; 
  epsil += k;
  dseta += k;

  ur += k;
  ut += k;  
  uz += k;
  du += k;  
  
  int k0=k%NZ;
  int t0=k/NZ;

  // X indices		
  double  kt=t0<NT/2 ? (double)t0 : (double)t0-(double)NT;
  double  kz=k0;


  //Fraction
  kt=(PI2/LT)*kt;
  kz=(PI2/LZ)*kz;	

  double	kk,ri2;

  //Local variables
  double delta_3,delta_1,delta_2,delta_4;
  double epsil_3,epsil_1,epsil_2,epsil_4;
  double dseta_3,dseta_1,dseta_2,dseta_4;
  
  double alpha_0;
  double gamma_0;
  double betha_0;

  double2 xT;

  double g_1,g_2,g_3;
  double2 y_1,y_2,y_3,y_4;
  double2 x_0,x_1,x_2,x_3;

  double a0,b0,c0,d0,e0,f0,g0;  
  double2 urh, uth, uzh, duh;

  
  //j=0 //Neuman boundary condtions

  d0=d1[0];
  e0=e1[0];
  f0=f1[0];
  g0=g1[0];
  
  delta_1=d0; //D[0]
  epsil_1=e0;
  dseta_1=f0;
  
  g_1=g0;
  
//   if(k0==1 & t0==1) printf("\ndelta = %d %e",1,delta_1);
//   if(k0==1 & t0==1) printf("\nepsil = %d %e",1,epsil_1);
//   if(k0==1 & t0==1) printf("\ndseda = %d %e",1,dseta_1);
//   if(k0==1 & t0==1) printf("\ng = %d %e",1,g0);

  delta[0*NT*NZ]=delta_1;
  epsil[0*NT*NZ]=epsil_1;
  dseta[0*NT*NZ]=dseta_1;

  //Boundary conditions dp/dr=0
  
  xT.x=0.0;
  xT.y=0.0; 

  y_1.x=xT.x;
  y_1.y=xT.y;

  y[0*NT*NZ]=y_1;


  //j=1

  ri2=rad[1]*rad[1]; 
//   if(k0==1 & t0==1) printf("\nc_0 = %d %e",1,ri2);

  kk=kt*kt/ri2 + kz*kz;
  
  c0=(c2[1]+1.0/rad[1]*c1[1]);
  d0=(d2[1]+1.0/rad[1]*d1[1]-kk);
  e0=(e2[1]+1.0/rad[1]*e1[1]);
  f0=(f2[1]+1.0/rad[1]*f1[1]);
  g0=(g2[1]+1.0/rad[1]*g1[1]);

  gamma_0=(c0)/delta_1;
//   if(k0==1 & t0==1) printf("\ngamma = %d %e",1,gamma_0);
//   if(k0==1 & t0==1) printf("\nc_0 = %d %e %e %e %e",1,c0,delta_1,gamma_0,d0);
//   if(k0==1 & t0==1) printf("\nall = %d %e %e %e",1,c2[1],rad[1],c1[1]);

  delta_2=d0-epsil_1*gamma_0;//D[1]
  epsil_2=e0-dseta_1*gamma_0;
  dseta_2=f0-g_1*gamma_0;
  g_2=g0;

  
//   if(k0==1 & t0==1) printf("\ndelta = %d %e",2,delta_2);
//   if(k0==1 & t0==1) printf("\nepsil = %d %e",2,epsil_2);
//   if(k0==1 & t0==1) printf("\nepsil = %d %e",2,dseta_2);

  
  delta[1*NT*NZ]=delta_2;
  epsil[1*NT*NZ]=epsil_2;
  dseta[1*NT*NZ]=dseta_2;

  urh=ur[1*NT*NZ];
  uth=ut[1*NT*NZ];
  uzh=uz[1*NT*NZ];
  duh=du[1*NT*NZ];
    
  //RHS
  xT.x=duh.x+urh.x/rad[1] - kz*uzh.y - kt*uth.y/rad[1];
  xT.y=duh.y+urh.y/rad[1] + kz*uzh.x + kt*uth.x/rad[1];

  
  y_2.x=xT.x-gamma_0*y_1.x;
  y_2.y=xT.y-gamma_0*y_1.y;

  y[1*NT*NZ]=y_2;

  //j=2
  
  ri2=rad[2]*rad[2]; 
  kk=kt*kt/ri2 + kz*kz;
  
  b0=(b2[2]+1.0/rad[2]*b1[2]);
  c0=(c2[2]+1.0/rad[2]*c1[2]);
  d0=(d2[2]+1.0/rad[2]*d1[2]-kk);
  e0=(e2[2]+1.0/rad[2]*e1[2]);
  f0=(f2[2]+1.0/rad[2]*f1[2]);
  g0=(g2[2]+1.0/rad[2]*g1[2]);

  betha_0=(b0)/delta_1;
  gamma_0=(c0-epsil_1*betha_0)/delta_2;

  delta_3=d0-epsil_2*gamma_0-dseta_1*betha_0;//D[2]
  epsil_3=e0-dseta_2*gamma_0-g_1*betha_0;
  dseta_3=f0-g_2*gamma_0;
  g_3=g0;   


  delta[2*NT*NZ]=delta_3;
  epsil[2*NT*NZ]=epsil_3;
  dseta[2*NT*NZ]=dseta_3;

  urh=ur[2*NT*NZ];
  uth=ut[2*NT*NZ];
  uzh=uz[2*NT*NZ];
  duh=du[2*NT*NZ];
    
  //RHS
  xT.x=duh.x+urh.x/rad[2] - kz*uzh.y - kt*uth.y/rad[2];
  xT.y=duh.y+urh.y/rad[2] + kz*uzh.x + kt*uth.x/rad[2];
  

  y_3.x=xT.x-betha_0*y_1.x-gamma_0*y_2.x;
  y_3.y=xT.y-betha_0*y_1.y-gamma_0*y_2.y;
  y[2*NT*NZ]=y_3;

  //Loop 2=3 to j=N-4

  for(int j=3; j<NR; j++){
  
    ri2=rad[j]*rad[j]; 
    kk=kt*kt/ri2 + kz*kz;
      
    a0=(a2[j]+1.0/rad[j]*a1[j]);
    b0=(b2[j]+1.0/rad[j]*b1[j]);
    c0=(c2[j]+1.0/rad[j]*c1[j]);
    d0=(d2[j]+1.0/rad[j]*d1[j]-kk);
    e0=(e2[j]+1.0/rad[j]*e1[j]);
    f0=(f2[j]+1.0/rad[j]*f1[j]);
    g0=(g2[j]+1.0/rad[j]*g1[j]);

    //Neuman boundary conditions dp/dr=0
    
    if(j==NR-1){

      a0=a1[NR-1];
      b0=b1[NR-1];
      c0=c1[NR-1];
      d0=d1[NR-1];
      e0=0.0;
      f0=0.0;
      g0=0.0;

//       if(k==50)printf("\n %e %e %e %e",a0,b0,c0,d0);
      
    }

    alpha_0= a0/delta_1;   
    betha_0=(b0-epsil_1*alpha_0)/delta_2;
    gamma_0=(c0-epsil_2*betha_0-dseta_1*alpha_0)/delta_3;
        

    delta_4 =d0-epsil_3*gamma_0-dseta_2*betha_0-g_1*alpha_0;
    epsil_4 =e0-dseta_3*gamma_0-g_2*betha_0;
    dseta_4 =f0-g_3*gamma_0;

//     if(k0==1 & t0==1) printf("\ndelta = %d %e",j,delta_4);

    
    //write
    delta[j*NT*NZ]=delta_4;
    delta_1=delta_2;
    delta_2=delta_3;
    delta_3=delta_4;

    //write
    epsil[j*NT*NZ]=epsil_4;
    epsil_1 =epsil_2;
    epsil_2 =epsil_3;
    epsil_3 =epsil_4;


    //write
    dseta[j*NT*NZ]=dseta_4;
    dseta_1 =dseta_2;
    dseta_2 =dseta_3;
    dseta_3 =dseta_4;

    //Read g
    g_1=g_2;
    g_2=g_3;
    g_3=g0;

    urh=ur[j*NT*NZ];
    uth=ut[j*NT*NZ];
    uzh=uz[j*NT*NZ];
    duh=du[j*NT*NZ];
    
    //RHS
    xT.x=duh.x + urh.x/rad[j] - kz*uzh.y - kt*uth.y/rad[j];
    xT.y=duh.y + urh.y/rad[j] + kz*uzh.x + kt*uth.x/rad[j];
        
    //Boundary conditions

    if(j==NR-1){

      xT.x=0.0;
      xT.y=0.0;

    }

//     if(k0==1 & t0==1) printf("\nRHS = %d %e %e",j,xT.x,xT.y);
//     if(k0==1 & t0==1) printf("\nRHS = %d %e %e",j,kt,kz);

    
    y_4.x=xT.x-alpha_0*y_1.x-betha_0*y_2.x-gamma_0*y_3.x;
    y_4.y=xT.y-alpha_0*y_1.y-betha_0*y_2.y-gamma_0*y_3.y;

    y[j*NT*NZ]=y_4;
    y_1=y_2;
    y_2=y_3;
    y_3=y_4;
  
  }
 

    //BACKWARD SUSTITUTION

    //N
    x_3.x=y_3.x/delta_3;
    x_3.y=y_3.y/delta_3;

    x[(NR-1)*NT*NZ].x=(double)x_3.x;
    x[(NR-1)*NT*NZ].y=(double)x_3.y;

    
    //N-1
    x_2.x=(y_2.x-epsil_2*x_3.x)/delta_2;
    x_2.y=(y_2.y-epsil_2*x_3.y)/delta_2;

    x[(NR-2)*NT*NZ].x=(double)x_2.x;    
    x[(NR-2)*NT*NZ].y=(double)x_2.y;    

    //N-3
    x_1.x=(y_1.x-epsil_1*x_2.x-dseta_1*x_3.x)/delta_1;
    x_1.y=(y_1.y-epsil_1*x_2.y-dseta_1*x_3.y)/delta_1;

    x[(NR-3)*NT*NZ].x=(double)x_1.x;    
    x[(NR-3)*NT*NZ].y=(double)x_1.y;    

    for(int j=NR-4;j>=0;j--){

    g0=(g2[j]+1.0/rad[j]*g1[j]);
	
    if(j==0) g0=g1[0];

    x_0.x =(y[j*NT*NZ].x-epsil[j*NT*NZ]*x_1.x
		-dseta[j*NT*NZ]*x_2.x-g0*x_3.x)/delta[j*NT*NZ];

    x_0.y =(y[j*NT*NZ].y-epsil[j*NT*NZ]*x_1.y
		-dseta[j*NT*NZ]*x_2.y-g0*x_3.y)/delta[j*NT*NZ];
    //Write
    x[j*NT*NZ].x=(double)x_0.x;
    x[j*NT*NZ].y=(double)x_0.y;

    x_3 =x_2;
    x_2 =x_1;
    x_1 =x_0;
    
    }  


}

__global__ void solve_ep_t_7_pressure_hom(double2* __restrict x,double2* __restrict y,
                                      double* __restrict delta, double* __restrict epsil, double* __restrict dseta
                                      ,const double* a1,const double*b1,const double*c1,const double*d1,const double*e1,const double*f1,const double* g1
                                      ,const double* a2,const double*b2,const double*c2,const double*d2,const double*e2,const double*f2,const double* g2
                                      ,const double* rad,double2 Bci, double2 Bco)
{

  // k=0...NX*NZ

  int k = threadIdx.x + blockDim.x*blockIdx.x;

  if(k>=NT*NZ) return;

  //Sets the pointers, each thread at a position

  x += k;
  y += k;  

  delta += k; 
  epsil += k;
  dseta += k;

  int k0=k%NZ;
  int t0=k/NZ;

  // X indices		
  double  kt=t0<NT/2 ? (double)t0 : (double)t0-(double)NT;
  double  kz=k0;


  //Fraction
  kt=(PI2/LT)*kt;
  kz=(PI2/LZ)*kz;	

  double	kk,ri2;

  
  //Local variables
  double delta_3,delta_1,delta_2,delta_4;
  double epsil_3,epsil_1,epsil_2,epsil_4;
  double dseta_3,dseta_1,dseta_2,dseta_4;
  
  double alpha_0;
  double gamma_0;
  double betha_0;

  double2 xT;

  double g_1,g_2,g_3;
  double2 y_1,y_2,y_3,y_4;
  double2 x_0,x_1,x_2,x_3;

  double a0,b0,c0,d0,e0,f0,g0;  

  //j=0 //Neuman boundary condtions

  d0=d1[0];
  e0=e1[0];
  f0=f1[0];
  g0=g1[0];
  
//   if(k==50)printf("\ncoef=%e,%e,%e,%e",d0,e0,f0,g0);
  
  delta_1=d0; //D[0]
  epsil_1=e0;
  dseta_1=f0;
  

  g_1=g0;
  
//   if(k0==1 & t0==1) printf("\ndelta = %d %e",1,delta_1);
//   if(k0==1 & t0==1) printf("\nepsil = %d %e",1,epsil_1);
//   if(k0==1 & t0==1) printf("\ndseda = %d %e",1,dseta_1);
//   if(k0==1 & t0==1) printf("\ng = %d %e",1,g0);

  delta[0*NT*NZ]=delta_1;
  epsil[0*NT*NZ]=epsil_1;
  dseta[0*NT*NZ]=dseta_1;

  //Boundary conditions dp/dr=0
  
  xT.x=Bci.x;
  xT.y=Bci.y; 

//   if(k==20);printf("\nBci=%e,%e",Bci.x,Bci.y);
  
  y_1.x=xT.x;
  y_1.y=xT.y;

  y[0*NT*NZ]=y_1;


  //j=1

  ri2=rad[1]*rad[1]; 
//   if(k0==1 & t0==1) printf("\nc_0 = %d %e",1,ri2);

  kk=kt*kt/ri2 + kz*kz;
  
  c0=(c2[1]+1.0/rad[1]*c1[1]);
  d0=(d2[1]+1.0/rad[1]*d1[1]-kk);
  e0=(e2[1]+1.0/rad[1]*e1[1]);
  f0=(f2[1]+1.0/rad[1]*f1[1]);
  g0=(g2[1]+1.0/rad[1]*g1[1]);

  gamma_0=(c0)/delta_1;
//   if(k0==1 & t0==1) printf("\ngamma = %d %e",1,gamma_0);
//   if(k0==1 & t0==1) printf("\nc_0 = %d %e %e %e %e",1,c0,delta_1,gamma_0,d0);
//   if(k0==1 & t0==1) printf("\nall = %d %e %e %e",1,c2[1],rad[1],c1[1]);

  delta_2=d0-epsil_1*gamma_0;//D[1]
  epsil_2=e0-dseta_1*gamma_0;
  dseta_2=f0-g_1*gamma_0;
  g_2=g0;

  delta[1*NT*NZ]=delta_2;
  epsil[1*NT*NZ]=epsil_2;
  dseta[1*NT*NZ]=dseta_2;
  
//   if(k0==1 & t0==1) printf("\ndelta = %d %e",2,delta_2);
//   if(k0==1 & t0==1) printf("\nepsil = %d %e",2,epsil_2);
//   if(k0==1 & t0==1) printf("\nepsil = %d %e",2,dseta_2);

  
  xT.x=0.0;
  xT.y=0.0;

  
  y_2.x=xT.x-gamma_0*y_1.x;
  y_2.y=xT.y-gamma_0*y_1.y;

  y[1*NT*NZ]=y_2;

  //j=2
  
  ri2=rad[2]*rad[2]; 
  kk=kt*kt/ri2 + kz*kz;
  
  b0=(b2[2]+1.0/rad[2]*b1[2]);
  c0=(c2[2]+1.0/rad[2]*c1[2]);
  d0=(d2[2]+1.0/rad[2]*d1[2]-kk);
  e0=(e2[2]+1.0/rad[2]*e1[2]);
  f0=(f2[2]+1.0/rad[2]*f1[2]);
  g0=(g2[2]+1.0/rad[2]*g1[2]);

  betha_0=(b0)/delta_1;
  gamma_0=(c0-epsil_1*betha_0)/delta_2;

  delta_3=d0-epsil_2*gamma_0-dseta_1*betha_0;//D[2]
  epsil_3=e0-dseta_2*gamma_0-g_1*betha_0;
  dseta_3=f0-g_2*gamma_0;
  g_3=g0;   


  delta[2*NT*NZ]=delta_3;
  epsil[2*NT*NZ]=epsil_3;
  dseta[2*NT*NZ]=dseta_3;

  xT.x=0.0;
  xT.y=0.0;
  

  y_3.x=xT.x-betha_0*y_1.x-gamma_0*y_2.x;
  y_3.y=xT.y-betha_0*y_1.y-gamma_0*y_2.y;
  y[2*NT*NZ]=y_3;

  //Loop 2=3 to j=N-4

  for(int j=3; j<NR; j++){
  
    ri2=rad[j]*rad[j]; 
    kk=kt*kt/ri2 + kz*kz;
      
    a0=(a2[j]+1.0/rad[j]*a1[j]);
    b0=(b2[j]+1.0/rad[j]*b1[j]);
    c0=(c2[j]+1.0/rad[j]*c1[j]);
    d0=(d2[j]+1.0/rad[j]*d1[j]-kk);
    e0=(e2[j]+1.0/rad[j]*e1[j]);
    f0=(f2[j]+1.0/rad[j]*f1[j]);
    g0=(g2[j]+1.0/rad[j]*g1[j]);

    //Neuman boundary conditions dp/dr=0
    
    if(j==NR-1){

      a0=a1[NR-1];
      b0=b1[NR-1];
      c0=c1[NR-1];
      d0=d1[NR-1];
      e0=0.0;
      f0=0.0;
      g0=0.0;

    }

    alpha_0= a0/delta_1;   
    betha_0=(b0-epsil_1*alpha_0)/delta_2;
    gamma_0=(c0-epsil_2*betha_0-dseta_1*alpha_0)/delta_3;
        

    delta_4 =d0-epsil_3*gamma_0-dseta_2*betha_0-g_1*alpha_0;
    epsil_4 =e0-dseta_3*gamma_0-g_2*betha_0;
    dseta_4 =f0-g_3*gamma_0;

//     if(k0==1 & t0==1) printf("\ndelta = %d %e",j,delta_4);

    
    //write
    delta[j*NT*NZ]=delta_4;
    delta_1=delta_2;
    delta_2=delta_3;
    delta_3=delta_4;

    //write
    epsil[j*NT*NZ]=epsil_4;
    epsil_1 =epsil_2;
    epsil_2 =epsil_3;
    epsil_3 =epsil_4;


    //write
    dseta[j*NT*NZ]=dseta_4;
    dseta_1 =dseta_2;
    dseta_2 =dseta_3;
    dseta_3 =dseta_4;

    //Read g
    g_1=g_2;
    g_2=g_3;
    g_3=g0;

    xT.x=0.0;
    xT.y=0.0;
        
    //Boundary conditions

    if(j==NR-1){

      xT.x=Bco.x;
      xT.y=Bco.y;

    }

//     if(k0==1 & t0==1) printf("\nRHS = %d %e %e",j,xT.x,xT.y);
//     if(k0==1 & t0==1) printf("\nRHS = %d %e %e",j,kt,kz);

    
    y_4.x=xT.x-alpha_0*y_1.x-betha_0*y_2.x-gamma_0*y_3.x;
    y_4.y=xT.y-alpha_0*y_1.y-betha_0*y_2.y-gamma_0*y_3.y;

    y[j*NT*NZ]=y_4;
    y_1=y_2;
    y_2=y_3;
    y_3=y_4;

  }
 

    //BACKWARD SUSTITUTION

    //N
    x_3.x=y_3.x/delta_3;
    x_3.y=y_3.y/delta_3;

    x[(NR-1)*NT*NZ].x=(double)x_3.x;
    x[(NR-1)*NT*NZ].y=(double)x_3.y;

    
    //N-1
    x_2.x=(y_2.x-epsil_2*x_3.x)/delta_2;
    x_2.y=(y_2.y-epsil_2*x_3.y)/delta_2;

    x[(NR-2)*NT*NZ].x=(double)x_2.x;    
    x[(NR-2)*NT*NZ].y=(double)x_2.y;    

    //N-3
    x_1.x=(y_1.x-epsil_1*x_2.x-dseta_1*x_3.x)/delta_1;
    x_1.y=(y_1.y-epsil_1*x_2.y-dseta_1*x_3.y)/delta_1;

    x[(NR-3)*NT*NZ].x=(double)x_1.x;    
    x[(NR-3)*NT*NZ].y=(double)x_1.y;    

    for(int j=NR-4;j>=0;j--){

    g0=(g2[j]+1.0/rad[j]*g1[j]);
	
    if(j==0) g0=g1[0];

    x_0.x =(y[j*NT*NZ].x-epsil[j*NT*NZ]*x_1.x
		-dseta[j*NT*NZ]*x_2.x-g0*x_3.x)/delta[j*NT*NZ];

    x_0.y =(y[j*NT*NZ].y-epsil[j*NT*NZ]*x_1.y
		-dseta[j*NT*NZ]*x_2.y-g0*x_3.y)/delta[j*NT*NZ];
    //Write
    x[j*NT*NZ].x=(double)x_0.x;
    x[j*NT*NZ].y=(double)x_0.y;

    x_3 =x_2;
    x_2 =x_1;
    x_1 =x_0;
    
    }  


}


__global__ void solve_ep_t_7_inplace(double* __restrict x,
              const double* aa,const double*bb,const double*cc,const double*dd,const double*ee,
              const double*ff,const double*gg)
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

//   if(k==50)printf("\ncoef=%e,%e,%e,%e",dd[0],ee[0],ff[0],gg[0]);
//   if(k==50)printf("\ncoef=%e,%e,%e,%e",xA,xB,xC,xD);

  //j=0
   xT=dd[0]*xA + ee[0]*xB + ff[0]*xC + gg[0]*xD;//RHS
//    if(xT<-0.5)printf("\nxT=%e",xT);
   x[0*sizep]=xT;

  //j=1
   xT=cc[1]*xA + dd[1]*xB + ee[1]*xC + ff[1]*xD + gg[1]*xE;//RHS
   x[1*sizep]=xT;

  //j=2
   xT=bb[2]*xA + cc[2]*xB + dd[2]*xC + ee[2]*xD + ff[2]*xE + gg[2]*xF;//RHS
   x[2*sizep]=xT;
  
   for(int i=3;i<NR;i++){
     
    xT=aa[i]*xA + bb[i]*xB + cc[i]*xC + dd[i]*xD + ee[i]*xE + ff[i]*xF + gg[i]*xG;//RHS
    x[i*sizep]=xT;
   
    xA = xB;
    xB = xC;
    xC = xD;
    xD = xE;
    xE = xF;
    xF = xG;
    
    if(i<NR-4){ //Read next elements
      xG = x[(i+4)*sizep];
    }
   
   }
}




static void setCoefficientsDerivatives(double* mesh_r){

  double* coef1=(double*)malloc(7*sizeof(double));
  double* coef2=(double*)malloc(7*sizeof(double));

  //FIRST
  CHECK_CUDART(cudaMalloc((void**)&a1,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&b1,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&c1,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&d1,NR*sizeof(double))); 
  CHECK_CUDART(cudaMalloc((void**)&e1,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&f1,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&g1,NR*sizeof(double)));


  //SECOND
  CHECK_CUDART(cudaMalloc((void**)&a2,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&b2,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&c2,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&d2,NR*sizeof(double))); 
  CHECK_CUDART(cudaMalloc((void**)&e2,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&f2,NR*sizeof(double)));
  CHECK_CUDART(cudaMalloc((void**)&g2,NR*sizeof(double)));

  double* aa_h=(double*)malloc(NR*sizeof(double));
  double* bb_h=(double*)malloc(NR*sizeof(double));
  double* cc_h=(double*)malloc(NR*sizeof(double));  
  double* dd_h=(double*)malloc(NR*sizeof(double));
  double* ee_h=(double*)malloc(NR*sizeof(double));
  double* ff_h=(double*)malloc(NR*sizeof(double));
  double* gg_h=(double*)malloc(NR*sizeof(double));

  ///////i=0///////
  
  fd_weigths(mesh_r[0],mesh_r,4,coef1,coef2);
  
  //RHS
  aa_h[0]=0.0;
  bb_h[0]=0.0;
  cc_h[0]=0.0;
  dd_h[0]=coef1[0];
  ee_h[0]=coef1[1];
  ff_h[0]=coef1[2];
  gg_h[0]=coef1[3];
  
  /////////////i=1///////////////////
  
  fd_weigths(mesh_r[1],mesh_r,5,coef1,coef2);
  
  //RHS
  aa_h[1]=0.0;
  bb_h[1]=0.0;
  cc_h[1]=coef1[0];
  dd_h[1]=coef1[1];
  ee_h[1]=coef1[2];
  ff_h[1]=coef1[3];
  gg_h[1]=coef1[4];


  /////////////i=2///////////////////

  fd_weigths(mesh_r[2],mesh_r,6,coef1,coef2);
  
  //RHS
  aa_h[2]=0.0;
  bb_h[2]=coef1[0];
  cc_h[2]=coef1[1];
  dd_h[2]=coef1[2];
  ee_h[2]=coef1[3];
  ff_h[2]=coef1[4];
  gg_h[2]=coef1[5];

  for(int i=3;i<NR-3;i++){

    fd_weigths(mesh_r[i],mesh_r+(i-3),7,coef1,coef2);

    //RHS
    aa_h[i]=coef1[0];
    bb_h[i]=coef1[1];
    cc_h[i]=coef1[2];
    dd_h[i]=coef1[3];
    ee_h[i]=coef1[4];
    ff_h[i]=coef1[5];
    gg_h[i]=coef1[6];


  } 


  /////////////i=NR-3///////////////////
  
  fd_weigths(mesh_r[NR-3],mesh_r+(NR-6),6,coef1,coef2);

  //RHS
  aa_h[NR-3]=coef1[0];
  bb_h[NR-3]=coef1[1];
  cc_h[NR-3]=coef1[2];
  dd_h[NR-3]=coef1[3];
  ee_h[NR-3]=coef1[4];
  ff_h[NR-3]=coef1[5];
  gg_h[NR-3]=0.0;


  /////////////i=NR-2///////////////////
  
  fd_weigths(mesh_r[NR-2],mesh_r+(NR-5),5,coef1,coef2);

  //RHS
  aa_h[NR-2]=coef1[0];
  bb_h[NR-2]=coef1[1];
  cc_h[NR-2]=coef1[2];
  dd_h[NR-2]=coef1[3];
  ee_h[NR-2]=coef1[4];
  ff_h[NR-2]=0.0;
  gg_h[NR-2]=0.0;

  
  /////////////i=NR-1///////////////////
  
  fd_weigths(mesh_r[NR-1],mesh_r+(NR-4),4,coef1,coef2);

  //RHS
  aa_h[NR-1]=coef1[0];
  bb_h[NR-1]=coef1[1];
  cc_h[NR-1]=coef1[2];
  dd_h[NR-1]=coef1[3];
  ee_h[NR-1]=0;
  ff_h[NR-1]=0;
  gg_h[NR-1]=0;

  


  //Check coeficients
  FILE* fp=fopen("coefficients_first_derivative_7.dat","w");
  for(int j=0;j<NR;j++){
    fprintf(fp,"\n %03d %+e %+e %+e %+e %+e %+e %+e %+e",j,
    mesh_r[j],aa_h[j],bb_h[j],cc_h[j],dd_h[j],ee_h[j],ff_h[j],gg_h[j]);
  } 
  fclose(fp);
  


//   //Check coeficients
//   size_t fsize; 
//   FILE* fp=fopen("coefficients_first_derivative.bin","wb");
//   fsize =fwrite( (unsigned char*)mesh_points,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)aa_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)bb_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)cc_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)dd_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)ee_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)ff_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)gg_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)hh_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)jj_h,sizeof(double),NR,fp);
//   fclose(fp);
  

  //FIRST

  CHECK_CUDART(cudaMemcpy((void*)a1,aa_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)b1,bb_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)c1,cc_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)d1,dd_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)e1,ee_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)f1,ff_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)g1,gg_h,NR*sizeof(double),cudaMemcpyHostToDevice));

  //SECOND DERIVATIVE
  
  ///////i=0///////
  
  fd_weigths(mesh_r[0],mesh_r,4,coef1,coef2);
  
  //RHS
  aa_h[0]=0.0;
  bb_h[0]=0.0;
  cc_h[0]=0.0;
  dd_h[0]=coef2[0];
  ee_h[0]=coef2[1];
  ff_h[0]=coef2[2];
  gg_h[0]=coef2[3];
  
  /////////////i=1///////////////////
  
  fd_weigths(mesh_r[1],mesh_r,5,coef1,coef2);
  
  //RHS
  aa_h[1]=0.0;
  bb_h[1]=0.0;
  cc_h[1]=coef2[0];
  dd_h[1]=coef2[1];
  ee_h[1]=coef2[2];
  ff_h[1]=coef2[3];
  gg_h[1]=coef2[4];


  /////////////i=2///////////////////

  fd_weigths(mesh_r[2],mesh_r,6,coef1,coef2);
  
  //RHS
  aa_h[2]=0.0;
  bb_h[2]=coef2[0];
  cc_h[2]=coef2[1];
  dd_h[2]=coef2[2];
  ee_h[2]=coef2[3];
  ff_h[2]=coef2[4];
  gg_h[2]=coef2[5];

  for(int i=3;i<NR-3;i++){

    fd_weigths(mesh_r[i],mesh_r+(i-3),7,coef1,coef2);

    //RHS
    aa_h[i]=coef2[0];
    bb_h[i]=coef2[1];
    cc_h[i]=coef2[2];
    dd_h[i]=coef2[3];
    ee_h[i]=coef2[4];
    ff_h[i]=coef2[5];
    gg_h[i]=coef2[6];


  } 


  /////////////i=NR-3///////////////////
  
  fd_weigths(mesh_r[NR-3],mesh_r+(NR-6),6,coef1,coef2);

  //RHS
  aa_h[NR-3]=coef2[0];
  bb_h[NR-3]=coef2[1];
  cc_h[NR-3]=coef2[2];
  dd_h[NR-3]=coef2[3];
  ee_h[NR-3]=coef2[4];
  ff_h[NR-3]=coef2[5];
  gg_h[NR-3]=0.0;


  /////////////i=NR-2///////////////////
  
  fd_weigths(mesh_r[NR-2],mesh_r+(NR-5),5,coef1,coef2);

  //RHS
  aa_h[NR-2]=coef2[0];
  bb_h[NR-2]=coef2[1];
  cc_h[NR-2]=coef2[2];
  dd_h[NR-2]=coef2[3];
  ee_h[NR-2]=coef2[4];
  ff_h[NR-2]=0.0;
  gg_h[NR-2]=0.0;

  
  /////////////i=NR-1///////////////////
  
  fd_weigths(mesh_r[NR-1],mesh_r+(NR-4),4,coef1,coef2);

  //RHS
  aa_h[NR-1]=coef2[0];
  bb_h[NR-1]=coef2[1];
  cc_h[NR-1]=coef2[2];
  dd_h[NR-1]=coef2[3];
  ee_h[NR-1]=0;
  ff_h[NR-1]=0;
  gg_h[NR-1]=0;

  


  //Check coeficients
  fp=fopen("coefficients_second_derivative_7.dat","w");
  for(int j=0;j<NR;j++){
    fprintf(fp,"\n %03d %+e %+e %+e %+e %+e %+e %+e %+e",j,
    mesh_r[j],aa_h[j],bb_h[j],cc_h[j],dd_h[j],ee_h[j],ff_h[j],gg_h[j]);
  } 
  fclose(fp);
  


//   //Check coeficients
//   size_t fsize; 
//   FILE* fp=fopen("coefficients_first_derivative.bin","wb");
//   fsize =fwrite( (unsigned char*)mesh_points,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)aa_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)bb_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)cc_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)dd_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)ee_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)ff_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)gg_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)hh_h,sizeof(double),NR,fp);
//   fsize =fwrite( (unsigned char*)jj_h,sizeof(double),NR,fp);
//   fclose(fp);
  

  //FIRST

  CHECK_CUDART(cudaMemcpy((void*)a2,aa_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)b2,bb_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)c2,cc_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)d2,dd_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)e2,ee_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)f2,ff_h,NR*sizeof(double),cudaMemcpyHostToDevice));
  CHECK_CUDART(cudaMemcpy((void*)g2,gg_h,NR*sizeof(double),cudaMemcpyHostToDevice));

  //Copy the mesh
  CHECK_CUDART(cudaMalloc((void**)&mesh_d,NR*sizeof(double)));
  CHECK_CUDART(cudaMemcpy((void*)mesh_d,mesh_r,NR*sizeof(double),cudaMemcpyHostToDevice));

  free(aa_h);
  free(cc_h);
  free(bb_h);
  free(dd_h);
  free(ee_h); 
  free(ff_h);
  free(gg_h);


  free(coef1);
  free(coef2);
}

void setImplic_7_exp(double* mesh,size_p sizes)
{
    setCoefficientsDerivatives(mesh);

}

void deriv_R_7E(double2* u)
{

    dim3 grid,block;
    block.x=128;
    grid.x=(2*NZ*NT + block.x - 1)/block.x;
    solve_ep_t_7_inplace<<<grid,block>>>((double*)u,a1,b1,c1,d1,e1,f1,g1);


}

void implicit_step(double2* u, double2* aux1,double* aux2, double* aux3, double* aux4, double d_implicit,double dt,int flag){
    

  double Bco, Bci;
  
  //Boundary conditions for mode 0
  
  if(flag==0){
    Bci= REYNOLDS_INNER;
    Bco= REYNOLDS_OUTER;
  }
  
  if(flag==1){
    Bci= -REYNOLDS_INNER;
    Bco= -REYNOLDS_OUTER;
  }
  
  if(flag==2){
    Bco=0.0;
    Bci=0.0;
  }
  
  dim3 grid,block;
  block.x=128;
  grid.x=(NT*NZ + block.x - 1)/block.x;
  solve_ep_t_7_implicit<<<grid,block>>>(u,aux1,aux2,aux3,aux4
  /*Coefficients*/,a1,b1,c1,d1,e1,f1,g1,a2,b2,c2,d2,e2,f2,g2,mesh_d,d_implicit,dt,flag,Bci,Bco);

}

void implicit_step_hom(double2* u, double2* aux1,double* aux2, double* aux3, double* aux4, double d_implicit,double dt,int flag,double2 Bci, double2 Bco){
    

  dim3 grid,block;
  block.x=128;
  grid.x=(NT*NZ + block.x - 1)/block.x;
  solve_ep_t_7_implicit_hom<<<grid,block>>>(u,aux1,aux2,aux3,aux4
  /*Coefficients*/,a1,b1,c1,d1,e1,f1,g1,a2,b2,c2,d2,e2,f2,g2,
  mesh_d,d_implicit,dt,flag,Bci,Bco);

}

void calcPressure(double2* p, vfield u, double2* dr,
                   double2* aux1,double* aux2, double* aux3, double* aux4){

  copyBuffer(u.r,dr);
  deriv_R_9E(dr);
  

  CHECK_CUDART( cudaDeviceSetCacheConfig(cudaFuncCachePreferL1) );

  dim3 grid,block;
  block.x=128;
  grid.x=(NT*NZ + block.x - 1)/block.x;
  solve_ep_t_7_pressure<<<grid,block>>>(p,aux1,u.r,u.t,u.z,dr,aux2,aux3,aux4
  /*Coefficients*/,a1,b1,c1,d1,e1,f1,g1,a2,b2,c2,d2,e2,f2,g2,mesh_d);

  CHECK_CUDART( cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));

  //Set zero mode to 0
  setZeromode(p);  
}

void calcPressure_hom(double2* u,double2* aux1,double* aux2, double* aux3, double* aux4
                      ,double2 Bci, double2 Bco){
  
  CHECK_CUDART( cudaDeviceSetCacheConfig(cudaFuncCachePreferL1) );
  
  dim3 grid,block;
  block.x=128;
  grid.x=(NT*NZ + block.x - 1)/block.x;
  solve_ep_t_7_pressure_hom<<<grid,block>>>(u,aux1,aux2,aux3,aux4
  /*Coefficients*/,a1,b1,c1,d1,e1,f1,g1,a2,b2,c2,d2,e2,f2,g2,mesh_d,Bci,Bco);

  CHECK_CUDART( cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));
  
  //Set zero mode to 0
  setZeromode(u);
  
}


