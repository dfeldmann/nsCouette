#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <time.h>
#include <sys/time.h>
#include <hdf5.h>
#include <hdf5_hl.h>


#define CHECK_CUDART(x) do { \
  cudaError_t res = (x); \
  if(res != cudaSuccess) { \
    fprintf(stderr, " CUDART: %s = %d (%s) at (%s:%d)\n",  #x, res, cudaGetErrorString(res),__FILE__,__LINE__); \
    exit(1); \
  } \
} while(0)

#define NZ 65 //Number of grid points in real space is 2*(NZ-1)
#define NT 64 //Number of grid points in real space is NT
#define NR 64 //Number of grid points in real space is NR

#define ETA (0.5) //eta=inner_radius/outer_radius, eta<1

#define NZP (3*(NZ-1)/2+1)
#define NTP (3*NT/2)

#define REYNOLDS_INNER (-700.0)
#define REYNOLDS_OUTER ( 700.0)

#define PI2 6.283185307179586
#define LT PI2
#define LZ (PI2/2)

//Number of steps to perform non-linear convolution: save memory
#define NSTEPS_CONV 1 

//Size of the influence matrix
#define MAT 8

//Parameter of the implicit time-integrator
#define C_IMPLICIT (0.51)

//Variable dt 
#define VARIABLE_DT 1
#define MAXDT 0.01
#define COURANT 0.25
#define TOLERANCE_DTERR 5e-5

typedef struct { double2* r;double2* t;double2* z;} vfield;
typedef struct { int Nr; int Nt; int Nz; double Lt ; double Lz;} size_p;

//Globals
extern double2* AUX_T;
extern double Cfl_r, Cfl_z, Cfl_th;

//Finite differences
void fd_weigths(double x0,double* x,int n,double* L1,double*L2);

//Cublas
void setCublas(void);
void transpose_A(double2* u_2,double2* u_1);
void transpose_B(double2* u_2,double2* u_1);
void transpose_infbackward(double* u_2,double* u_1,int mr);
void transpose_infforward(double* u_2,double* u_1,int mr);
void invert_infmatrix(double* src);

//Derivatives compact 3-3
void setDeriv_3_comp(double* mesh,size_p sizes);
void deriv_R_3C(double2* u,double2* aux);
void deriv_RR_3C(double2* u, double2* aux);

//Derivatives explicit 9
void setDeriv_9_exp(double* mesh,size_p sizes);
void deriv_R_9E(double2* u);
void deriv_RR_9E(double2* u);

//Implicit steps explicit 7
void setImplic_7_exp(double* mesh,size_p sizes);
void implicit_step(double2* u, double2* aux1,double* aux2, double* aux3, double* aux4, 
                   double d_implicit,double dt,int flag);
void calcPressure(double2* p, vfield u, double2* dr,
                   double2* aux1,double* aux2, double* aux3, double* aux4);
void deriv_R_7E(double2* u);

//IO
void writeBufferBinary(double* w,const char* file,size_t elsize, size_t elements);
void writeBufferBinary(double* w,const char* file,size_t elsize, size_t elements);
void writeBufferBinaryGPU(double* w,const char* file,size_t elsize,size_t elements);
void writeCheckpoint(vfield u, double* grid_mesh,double* time,char* fileName);
void readCheckpoint(vfield u, double* grid_mesh,double* time,char* fileName);
void writeFieldVis(vfield u, double* grid_mesh,double* time,char* fileName);

//fft
void setFft(size_p sizes);
void fftDestroy(void);
void fftForward(double2* buffer);
void fftBackward(double2* buffer);
void fftForward_reduced(double2* buffer);
void fftBackward_reduced(double2* buffer);


//Utils
void copyVfield(vfield u1, vfield u2);
void onesGPU(double2* u1, size_t elements);
void normalize(double2* u, double norm, size_t elements);
void decoupleForward(vfield u);
void decoupleBackward(vfield u);
void copyBuffer(double2* u1, double2* u2);
void updateNonlinear(vfield n1, vfield n2, double d_implicit);
void calcnorm(double2* u_out, double2* u_in);
void setZeromode(double2* u);
void initField(vfield u, double* mesh);
void writeBufferBinaryGPUreal(double* w,const char* file,
                              size_t elsize,size_t elements);
void zerosGPU(double2* u1, size_t elements);
void diffAbs(double2* u_diff, double2* u);
void fieldAbs(double2* u_diff, double2* u);


//Non-linear
void setNonlinear(double* mesh);
void nonLinear(vfield u, vfield dur,int flag, double time);


//Check routines
void checkDerivatives(void);
void checkNonlinear(void);
void checkImplicit(void);
void checkPressure(void);
void checkLaplacian(void);
void checkTranspose(void);
void checkBoundaries(void);

//Linear 
void setLinear(double* mesh);
void pressureProject(vfield rhs);
void implictStep(vfield u, vfield rhs, double d_implicit, double dt);
void calcLaplacian(vfield rhs, vfield u, double d_implicit, double dt);
void calcPressure_hom(double2* u,double2* aux1,double* aux2, double* aux3, double* aux4
                      ,double2 Bci, double2 Bco);
void implicit_step_hom(double2* u, double2* aux1,double* aux2, double* aux3, 
                        double* aux4, double d_implicit,double dt,int flag,double2 Bci, double2 Bco);
void calcDivergence(double2* div, vfield u);
void calcVorticity(vfield vor, vfield u);
void forcing(double2* uz);

//Boundary conditions
void setBoundary(double* mesh);
void calcHsol(double* M, vfield hsolu_i, vfield hsolu_o, vfield hsolp_i, vfield hsolp_o,
              double d_implicit, double dt);
void imposeBoundarycond(vfield u, vfield hsolu_i, vfield hsolu_o, vfield hsolp_i,
                        vfield hsolp_o, double* M);

//Time integrator
void integrate(vfield u, vfield u_w, vfield rhs, vfield rhs_w,
               int nsteps, double dt,double* time);
void setIntegrator(size_p sizes);
void check_convergence(int* iter, double dt);
void new_step(double* dt);
void measurecorr(vfield up, vfield um);

//Statistics
void setStatistics(double* mesh);
void readTorq(double2* ut,double2* utdr,double time);

//Non-linear
void padForward(double2* aux,double2* u);
void padBackward(double2* u,double2* aux);
void calcDut(vfield u, vfield du);
void calcDuz(vfield u, vfield du);

//Check cublas
void cublasCheck(cublasStatus_t error, const char* function);                                                               

