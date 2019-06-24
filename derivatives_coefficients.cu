// #include "TayC.h"
// 
// //Highly non-compact, non-optimised function to calculate the coefficients
// //of the FD with Lagrange pol. Just to make sure I understand.
// //x0, point to evaluate
// //x[n] , points of the function
// //n , length of the stencil
// //L1[n], coefficients of the first deriv
// //L2[n], coefficients of the secod deriv
// 
// void fd_weigths(double x0,double* x,int n,double* L1,double*L2){
//   
// //    double* L =(double*)malloc(n*sizeof(double));
// //    double* L1=(double*)malloc(n*sizeof(double));
// //    double* L2=(double*)malloc(n*sizeof(double));
// // 
// //     
// //     //Zero polinomia
// //     for(int i=0;i<n;i++){
// //       L[i]=1.0;
// //       for(int j=0;j<n;j++){
// //         if(i!=j){
// //           L[i]=L[i]*(x0-x[j])/(x[i]-x[j]);
// //         }
// //       }
// //     }
//         
//     //First polinomia
//     for(int i=0;i<n;i++){
//       L1[i]=0.0;
//       for(int j=0;j<n;j++){
//         if(i!=j){
//           double prod=1.0;
//           for(int m=0;m<n;m++){
//             if(m!=i && m!=j){
//               prod=prod*(x0-x[m])/(x[i]-x[m]);
//             }
//           }
//           L1[i]=L1[i]+prod/(x[i]-x[j]);
//         }
//       }
//     }
//   
//     //Second polinomia
//     for(int i=0;i<n;i++){
//       L2[i]=0.0;
//       for(int j=0;j<n;j++){
//         if(i!=j){
//           double sum=0.0;
//           for(int m=0;m<n;m++){
//             if(m!=i && m!=j){
//               double prod=1.0;
//               for(int l=0;l<n;l++){
//                 if(l!=i && l!=j && l!=m){
//                   prod=prod*(x0-x[l])/(x[i]-x[l]);
//                 }  
//               }  
//             sum=sum+prod/(x[i]-x[m]);
//             }
//           }
//           L2[i]=L2[i]+sum/(x[i]-x[j]);
//         }
//       }
//     }
//   
// //     printf("\n");
// //     for(int i=0;i<n;i++){
// //       printf("\n%e",L1[i]);
// //     }
// //     
// //     printf("\n");
// //     for(int i=0;i<n;i++){
// //       printf("\n%e",L2[i]);
// //     }
// //     
// //     free(L);
// //     free(L1);
// //     free(L2);
// /*    
//     c[0] = 1.0;
// //     printf("\nMn=%d",mn);
//     for(int i=1;i<n;i++){
// //        printf("\ni,m=%d,%d",i,m);
//        mn = min(i,m);//printf("\nmn=%d",mn);
//        c2 = 1.0;
//        c5 = c4;
//        c4 = x[i] - x0;
//        for(int j=0;j<i;j++){
// //           printf("j=%d,%d",j,i);exit(1);
//           c3 = x[i] - x[j];
//           c2 = c2*c3;
//           if(j==i-1){
//              for(int k = mn;k>0;k--){
// //                          printf("k=%d",k);exit(1);
//                 c[i+n*k] = c1*((double)k*c[i-1 + n*(k-1)]-c5*c[i-1+n*k])/c2;
//              }
//              c[i] = -c1*c5*c[i-1]/c2;
//           }
// //           printf("j=%d,%d",j,i);exit(1);
//           for(int k = mn;k>0;k--){
// //              printf("\nk,j,i,mn=%d,%d,%d,%d",k,j,i,mn);
//              c[j+n*k] = (c4*c[j+n*k]-(double)k*c[j+n*(k-1)])/c3;
//           }
// 
//           c[j] = c4*c[j]/c3;
//        }
// //        printf("k=%d",i);exit(1);
// 
//        c1 = c2;
//     }
// */
//     return;
// } 
// 
// void setCoefficientsDerivatives(double* mesh_points){
//   
//   
//   
//   
//   
// }
// 
// 
// 
// 
// void checkCoef(void){
//   
//   int n=9;
// 
//   double* L1=(double*)malloc(sizeof(double)*n);
//   double* L2=(double*)malloc(sizeof(double)*n);
//   double* x =(double*)malloc(sizeof(double)*n);
// 
//   double x0;
//  
//   for(int i=0;i<n;i++){
//   x[i]=(double)(i-4);
//   }
//   
//   x0=0;
//   printf("\nChecking coeff");
//   fd_weigths(x0,x,n,L1,L2);
// 
// //   x0=-4;
// //   n=4;
//   
//   for(int i=0;i<n;i++){
//   x[i]=(double)(i);
//   }
//   x0=0;
//   fd_weigths(x0,x,5,L1,L2);
//   
//   
//    for(int i=0;i<5;i++){
//      printf("\n%d,%f",i,L1[i]);
//    }
// 
//   for(int i=0;i<5;i++){
//      printf("\n%d,%f",i,L2[i]);
//    }
// }
