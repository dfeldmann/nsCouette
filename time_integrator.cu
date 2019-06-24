#include "TayC.h"

vfield u_w, rhs_w;
vfield hsolu_i, hsolu_o, hsolp_i, hsolp_o;
vfield u_init;

double* M;
double2* aux;

//Time integrator 
static double dterr=0, corr_dt=0, cfl=0;
static double lasterr=0;
static int counter;

//Cublas handle
static cublasHandle_t cublasHandle;

void setIntegrator(size_p sizes){
       
  
    size_t size_p=NR*NT*NZ*sizeof(double2);

    CHECK_CUDART(cudaMalloc(&hsolu_i.r,size_p));
    CHECK_CUDART(cudaMalloc(&hsolu_i.t,size_p));
    CHECK_CUDART(cudaMalloc(&hsolu_i.z,size_p));
    
    CHECK_CUDART(cudaMalloc(&hsolu_o.r,size_p));
    CHECK_CUDART(cudaMalloc(&hsolu_o.t,size_p));
    CHECK_CUDART(cudaMalloc(&hsolu_o.z,size_p));
    
    CHECK_CUDART(cudaMalloc(&hsolp_i.r,size_p));
    CHECK_CUDART(cudaMalloc(&hsolp_i.t,size_p));
    CHECK_CUDART(cudaMalloc(&hsolp_i.z,size_p));
    
    CHECK_CUDART(cudaMalloc(&hsolp_o.r,size_p));
    CHECK_CUDART(cudaMalloc(&hsolp_o.t,size_p));
    CHECK_CUDART(cudaMalloc(&hsolp_o.z,size_p));
    
    CHECK_CUDART(cudaMalloc(&u_init.r,size_p));
    CHECK_CUDART(cudaMalloc(&u_init.t,size_p));
    CHECK_CUDART(cudaMalloc(&u_init.z,size_p));
    
    CHECK_CUDART(cudaMalloc(&aux,size_p));
    
    CHECK_CUDART(cudaMalloc(&M,NZ*NT*MAT*MAT*sizeof(double)));
    
    cublasCheck(cublasCreate(&cublasHandle),"Cre");

    
}

void pretime(vfield hsolu_i, vfield hsolu_o, 
             vfield hsolp_i, vfield hsolp_o,
             double* M, double dt)
{
  
     // Calc homogenous solutions
     calcHsol(M, hsolu_i, hsolu_o, hsolp_i,
              hsolp_o, C_IMPLICIT, dt);
  
} 

void predictor(vfield u, vfield u_w, vfield rhs, 
               vfield rhs_w, double dt)
{

    //Save buffers
    copyVfield(u, u_w);
    copyVfield(rhs, rhs_w);
    
    //Decouple u.r and u.t in u_plus and u_minus
    decoupleForward(u);
    decoupleForward(u_w);
    decoupleForward(rhs);

    //Add the Laplacian
    calcLaplacian(rhs,u,C_IMPLICIT,dt);
  
    //Project pressure
    pressureProject(rhs);
    
    //Implicit step
    implictStep(u,rhs,C_IMPLICIT,dt);
    
    //Impose boundary conditions
    imposeBoundarycond(u,hsolu_i,hsolu_o,
                       hsolp_i,hsolp_o,M);
    
    //Move back to ur and ut
    decoupleBackward(u);
    
    //Set ur(k=0) to zero
    setZeromode(u.r);

}


void corrector(vfield u, vfield u_w, vfield rhs, 
               vfield rhs_w, double dt)
{
    //Copy field to check convergence
    copyVfield(u, u_init);
    
    //Nnew=c*Nnew + (1-c)*Nold
    updateNonlinear(rhs,rhs_w,C_IMPLICIT);
    
    //Decouple non-linear term
    decoupleForward(rhs);

    //Add the laplacian
    calcLaplacian(rhs,u_w,C_IMPLICIT,dt);

    //Project pressure
    pressureProject(rhs);
    
    //Implicit step
    implictStep(u,rhs,C_IMPLICIT,dt);
    
    //Impose boundary conditions
    imposeBoundarycond(u,hsolu_i,hsolu_o,
                       hsolp_i,hsolp_o,M);
    
    //Move back to ur and ut
    decoupleBackward(u);
        
    //Set ur(k=0) to zero
    setZeromode(u.r);
    
    //Check whether the difference is within the tolerance demanded
    measurecorr(u,u_init);
    
}

void integrate(vfield u, vfield u_w, vfield rhs, vfield rhs_w,
               int nsteps, double dt,double* time)
{

    pretime(hsolu_i, hsolu_o, hsolp_i, hsolp_o, M, dt);
    
    counter=0;

    while(counter<nsteps){  
    
      //Calc non-linear terms
      nonLinear(u,rhs,1,*time);
      
      if(VARIABLE_DT){
        new_step(&dt);
        pretime(hsolu_i, hsolu_o, hsolp_i, hsolp_o, M, dt);
      }

      //Predictor
      predictor(u, u_w, rhs, rhs_w, dt);

      int iter=1;
      
      while(iter!=0){
        
        //Calc non-linear
        nonLinear(u,rhs,0,*time);

        //Predictor
        corrector(u,u_w,rhs,rhs_w,dt);
        
        check_convergence(&iter,dt);
  
      }
      
      counter++;
      *time=*time + dt;
      
    }

    
    return;
    
}

void new_step(double* dt)
{

  static int ind = 0;
  int new_dt=0;
  double d  = max(abs(REYNOLDS_INNER-REYNOLDS_OUTER),1.0);
  double deltat = min(*dt*1.11, MAXDT/d);

  printf("\nCfl(t,r,z)=(%e,%e,%e)",Cfl_th,Cfl_r,Cfl_z);

  cfl=min(Cfl_th,Cfl_r); 
  cfl=min(cfl,Cfl_z);
  
  if(counter==0 & corr_dt==0) deltat = min(deltat, cfl*0.1);
  if(cfl > 0.0)  deltat = min(deltat, cfl*COURANT);
  if(corr_dt > 0.0)  deltat = min(deltat, corr_dt*0.95);

  ind = ind - 1;
  if(deltat < *dt*0.95 || (ind < 0 & deltat >  *dt*1.1)) new_dt=1;
      
  if(new_dt){
      if(deltat > *dt*1.1) ind = 100;
      *dt = deltat;
  }
  
  printf("\nDt=%e",*dt);

}

void check_convergence(int* iter, double dt)
{

  //Measure correction of dt                                                                                                                                                                                                                                             
  double d  = max(abs(REYNOLDS_INNER-REYNOLDS_OUTER),1.0);
  
  if(*iter==1){ 
    corr_dt = dt*sqrt(TOLERANCE_DTERR*d/dterr);
    lasterr = 1e50;
  }
  
  if(dt<1e-9 & counter>30){
     printf("check_convergence: dt --> 0 ");
     exit(1);
  }else if(*iter> 10){
     printf("check_convergence: too many its");
     exit(1);
  }else if(dterr>lasterr){
     printf("check_convergence: increasing error!!!");
     if(dterr > 2.0*TOLERANCE_DTERR*d) exit(1);
     if(dterr < 2.0*TOLERANCE_DTERR*d) corr_dt = dt/(10.0*COURANT);
     *iter = 0;
  }else if(dterr > TOLERANCE_DTERR){
     lasterr = dterr;
     *iter = *iter + 1;
  }else{
     if(VARIABLE_DT){
          printf("\nStep,dt=%d, %e",counter,dt);
     }else{
          printf("\nStep,it=%d, %d",counter,dt);
     }
     *iter = 0;
  }
  
  dterr = 0.0;
  
  return;
}

void measurecorr(vfield u, vfield u_init)
{
  
  int size_l=2*NT*NR*NZ;
  int index;
  
  double dt_max,ddt_max;
  double dz_max,ddz_max;
  double dr_max,ddr_max;
  double dd_max,d_max;
  
  fieldAbs(aux,u.r);
  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)aux,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&dr_max,(double*)aux + index-1, sizeof(double), cudaMemcpyDeviceToHost));
  
  fieldAbs(aux,u.t);
  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)aux,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&dt_max,(double*)aux + index-1, sizeof(double), cudaMemcpyDeviceToHost));
  
  fieldAbs(aux,u.z);
  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)aux,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&dz_max,(double*)aux + index-1, sizeof(double), cudaMemcpyDeviceToHost));

  diffAbs(u_init.r,u.r);
  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)u_init.r,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&ddr_max,(double*)u_init.r + index-1, sizeof(double), cudaMemcpyDeviceToHost));
  
  diffAbs(u_init.t,u.t);
  cublasCheck(cublasIdamax (cublasHandle,size_l, (const double *)u_init.t,1,&index),"Isa");
  CHECK_CUDART(cudaMemcpy(&ddt_max,(double*)u_init.t + index-1, sizeof(double), cudaMemcpyDeviceToHost));
  
  diffAbs(u_init.z,u.z);
  CHECK_CUDART(cudaMemcpy(&ddz_max,(double*)u_init.z + index-1, sizeof(double), cudaMemcpyDeviceToHost));

  dd_max=max(ddr_max,ddt_max);dd_max=max(dd_max,ddz_max);
  
  d_max =max(dr_max,dt_max);d_max=max(d_max,dz_max);

//   printf("\ndd_max,d_max=%e,%e",d_max,dd_max);
  
  if(d_max>0){
      dterr= dd_max/d_max;
  }else{
      dterr= 0.0;
  }

return;
  
}
