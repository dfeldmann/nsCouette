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

/*static cufftHandle fft2_r2c; 
static cufftHandle fft2_c2r;*/ 

static cufftHandle fft2_r2c_zp; 
static cufftHandle fft2_c2r_zp;

static cufftHandle fft2_c2r;
static cufftHandle fft2_r2c;

static double2* aux_fft;

static size_t sizep;

void cufftCheck( cufftResult error, const char* function )
{
	if(error != CUFFT_SUCCESS)
	{
		printf("\n error  %s : %d \n", function, error);
		exit(1);
	}
		
	return;
}  

void setFft(size_p sizes)
{

    sizep = sizeof(double2)*NR*NZP*NTP;

    int nzp[2]={NTP,2*NZP-2};

    cufftCheck(cufftPlanMany( &fft2_r2c_zp,2,nzp,NULL,1,0,NULL,1,0,CUFFT_D2Z,NR/NSTEPS_CONV),"ALLOCATE_FFT3_R2C_ZP");
    cufftCheck(cufftPlanMany( &fft2_c2r_zp,2,nzp,NULL,1,0,NULL,1,0,CUFFT_Z2D,NR/NSTEPS_CONV),"ALLOCATE_FFT3_C2R_ZP");

    int nz[2]={NT,2*NZ-2};
    
    cufftCheck(cufftPlanMany( &fft2_c2r,2,nz,NULL,1,0,NULL,1,0,CUFFT_Z2D,NR),"ALLOCATE_FFT3_R2C_ZP");
    cufftCheck(cufftPlanMany( &fft2_r2c,2,nz,NULL,1,0,NULL,1,0,CUFFT_D2Z,NR),"ALLOCATE_FFT3_C2R_ZP");

    
    if(NR%NSTEPS_CONV!=0){
        printf("\nError tama?os:NR must be divisible by NSTEPS_CONV");exit(1);
    }

    CHECK_CUDART(cudaMalloc((void**)&aux_fft,sizep));
      
    return;
}

void fftDestroy(void)
{
  cufftDestroy(fft2_r2c_zp);
  cufftDestroy(fft2_c2r_zp);

  return;
}

void fftForward(double2* buffer)
{
//   cufftCheck(cufftExecD2Z(fft2_r2c_zp,(double*)buffer,(double2*)aux_fft),"forward transform_zp");
//   normalize(aux_fft,(double)NTP*(2*NZP-2));
//   CHECK_CUDART(cudaMemcpy(buffer,aux_fft,sizep,cudaMemcpyDeviceToDevice));

     cufftCheck(cufftExecD2Z(fft2_r2c_zp,(double*)buffer,(double2*)buffer),"backward transform_zp");
     normalize(buffer,(double)NTP*(2*NZP-2),NR*NTP*NZP);

    //Normalize
    
  return;
}

void fftBackward(double2* buffer)
{
//   cufftCheck(cufftExecZ2D(fft2_c2r_zp,(double2*)buffer,(double*)aux_fft),"backward transform_zp");
//   CHECK_CUDART(cudaMemcpy(buffer,aux_fft,sizep,cudaMemcpyDeviceToDevice));

  cufftCheck(cufftExecZ2D(fft2_c2r_zp,(double2*)buffer,(double*)buffer),"backward transform_zp");
  //normalize(buffer,(double)NTP*(2*NZP-2),NR*NTP*NZP);

  //   CHECK_CUDART(cudaMemcpy(buffer,aux_fft,sizep,cudaMemcpyDeviceToDevice));
  
  return;

}

void fftBackward_reduced(double2* buffer)
{
  
  cufftCheck(cufftExecZ2D(fft2_c2r,(double2*)buffer,(double*)buffer),"backward transform");
//   normalize(buffer,(double)NT*(2*NZ-2),NR*NT*NZ);
  
  return;

}

void fftForward_reduced(double2* buffer)
{
  
  cufftCheck(cufftExecD2Z(fft2_r2c,(double*)buffer,(double2*)buffer),"backward transform");
  normalize(buffer,(double)NT*(2*NZ-2),NR*NT*NZ);

  return;

}
