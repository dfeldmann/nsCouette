#include"TayC.h"

int main( int argc, const char* argv[]){
  
  printf("\n+++++++++++++++++++++++++++");
  printf("\nStarting GPU Taylor-Couette");
  printf("\n+++++++++++++++++++++++++++");
  
  int dev=7;
  double time;

  printf("\n");
  printf("\nSetting device: %d",dev);
  CHECK_CUDART(cudaSetDevice(dev));
  
  //Timing variables
  clock_t start, end;
  double gpu_time_used;

  //Set up
  size_p sizes;
    
  sizes.Nr=NR;
  sizes.Nt=NT;
  sizes.Nz=NZ;
    
  //Initialize the mesh  
  double* grid_mesh=(double*)malloc(NR*sizeof(double));

  double r_i = 1.0;         //ETA/(1.0-ETA);
  double r_o = 1.0/ETA;     //1.0/(1.0-ETA);
  
  
  //Chebyshev mesh in r
  
  for(int j=0;j<NR;j++){
      //grid_mesh[j]=2.0*(double)j/(double)(NR-1)+1.0;
        grid_mesh[j]= (r_i+r_o)/2.0 - cos(PI2/2.0*j/(NR-1))/2.0;
    
  }
  
  //Write mesh

  //Set modules
  setImplic_7_exp(grid_mesh,sizes);
  setDeriv_9_exp(grid_mesh,sizes);
  setBoundary(grid_mesh);
  setNonlinear(grid_mesh);
  setFft(sizes);
  setCublas();
  setIntegrator(sizes);
  setLinear(grid_mesh);
  setStatistics(grid_mesh);
  
  //Allocate memory buffers
  vfield u, uw, rhs, rhsw;
  
  size_t size_p=NR*NT*NZ*sizeof(double2);

  CHECK_CUDART(cudaMalloc(&u.r,size_p));
  CHECK_CUDART(cudaMalloc(&u.t,size_p));
  CHECK_CUDART(cudaMalloc(&u.z,size_p));
  
  CHECK_CUDART(cudaMalloc(&uw.r,size_p));
  CHECK_CUDART(cudaMalloc(&uw.t,size_p));
  CHECK_CUDART(cudaMalloc(&uw.z,size_p));
  
  CHECK_CUDART(cudaMalloc(&rhs.r,size_p));
  CHECK_CUDART(cudaMalloc(&rhs.t,size_p));
  CHECK_CUDART(cudaMalloc(&rhs.z,size_p));
  
  CHECK_CUDART(cudaMalloc(&rhsw.r,size_p));
  CHECK_CUDART(cudaMalloc(&rhsw.t,size_p));
  CHECK_CUDART(cudaMalloc(&rhsw.z,size_p));

  
  //Start initial field
  initField(u, grid_mesh);
  
  //Or read check point
  //readCheckpoint(u,grid_mesh,&time,"./checkpoint.h5");

  //Time-step
  printf("\nParameter");
  printf("\nr_o,r_i=%e,%e",r_o,r_i);
  printf("\nRe_i,Re_o=%e,%e",REYNOLDS_OUTER,REYNOLDS_INNER);
  
  double U_i=REYNOLDS_OUTER*r_o;
  double U_o=REYNOLDS_INNER*r_i;
  
  double t_i=PI2*r_i/NT/U_i;
  double t_o=PI2*r_o/NT/U_o;
 
  printf("\nt_i,t_o=%e,%e",t_o,t_i);
  
  double dt=0.1*min(abs(t_i),abs(t_o));
  
  if(VARIABLE_DT){
    printf("\nRunning with fixed Courant number=%e",COURANT);
  }else{
    printf("\nRunning with fixed Dt=%e",dt);
  }
  
  //Number of steps
  int nsteps=200;
  
  start=clock(); 

  integrate(u, uw, rhs, rhsw, nsteps, dt,&time);
  
  CHECK_CUDART(cudaDeviceSynchronize());
  end=clock();  

  time=dt*nsteps;

  printf("\nTime_per_step=%e s",((double)(end-start))/CLOCKS_PER_SEC/nsteps);  

  //writeCheckpoint(u,grid_mesh,&time,"./checkpoint.h5");
  //writeFieldVis(u,grid_mesh,&time,"./visualiza.h5");


  return 0;
  
}
