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

static double  Tinner, Touter;
static double2 dut_i, dut_o;
static double2 ut_i, ut_o;
static double  r_i, r_o;

void setStatistics(double* mesh){
  
  r_i=mesh[0];
  r_o=mesh[NR-1];
  
}

void readTorq(double2* ut,double2* utdr,double time){
  
  
  //Reads the derivatives of the zero mode of dut and ut
  CHECK_CUDART(cudaMemcpy(&dut_i,(double2*)utdr +      0*NT*NZ,sizeof(double2),cudaMemcpyDeviceToHost));
  CHECK_CUDART(cudaMemcpy(&dut_o,(double2*)utdr + (NR-1)*NT*NZ,sizeof(double2),cudaMemcpyDeviceToHost));

  CHECK_CUDART(cudaMemcpy(&ut_i,(double2*)ut +      0*NT*NZ,sizeof(double2),cudaMemcpyDeviceToHost));
  CHECK_CUDART(cudaMemcpy(&ut_o,(double2*)ut + (NR-1)*NT*NZ,sizeof(double2),cudaMemcpyDeviceToHost));
  
  Tinner=-PI2*r_i*(dut_i.x - ut_i.x/r_i);
  Touter=-PI2*r_o*(dut_o.x - ut_o.x/r_o);
  
  double C_1 = (REYNOLDS_INNER/r_o - REYNOLDS_OUTER/r_i)/(r_i/r_o - r_o/r_i);
  double C_2 = (REYNOLDS_INNER*r_o - REYNOLDS_OUTER*r_i)/(r_o/r_i - r_i/r_o);

  double ut_l_i = C_1*r_i + C_2/r_i;
  double ut_l_o = C_1*r_o + C_2/r_o;
  
  double dut_l_i = C_1 - C_2/(r_i*r_i);
  double dut_l_o = C_1 - C_2/(r_o*r_o);
  
  double Tinner_l=-PI2*r_i*(dut_l_i - ut_l_i/r_i);
  double Touter_l=-PI2*r_o*(dut_l_o - ut_l_o/r_o);
  
  //Torque normalized by the laminar torque
  
  printf("\nTorque(r_i,r_o)=(%e,%e)",Tinner/Tinner_l,Touter/Touter_l);
  
  FILE *fp;
  fp=fopen("./torque.dat","a");
  fprintf(fp,"\n%+e %+e %+e",time,Tinner/Tinner_l,Touter/Touter_l);
  fclose(fp);
  
}
