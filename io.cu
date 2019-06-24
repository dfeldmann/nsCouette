#include "TayC.h"

float* host_buffer;
float* host_1;
float* host_2;
float* host_3;

FILE *fp;

void writeBufferBinary(double* w,const char* file,size_t elsize, size_t elements)
{
  
  size_t size=elements;
  fp=fopen(file,"wb");
  if(fp==NULL){printf("\nerror escritura: %s",file); exit(1);}
  size_t fsize =fwrite( (unsigned char*)w,elsize,size,fp);
  if(fsize!=size){ printf("\nwriting error: %s",file); exit(1);}
  fclose(fp);
  
}

void readBufferBinary(double* w,const char* file,size_t elsize,size_t elements)
{

  
  size_t size=elements;
  fp=fopen(file,"rb");
  if(fp==NULL){printf("\nerror lectura: %s",file); exit(1);}
  size_t fsize =fread( (unsigned char*)w,elsize,size,fp);
  if(fsize!=size){ printf("\nreading error: %s",file); exit(1);}
  fclose(fp);
    
  return;

}

void writeBufferBinaryGPU(double* w,const char* file,size_t elsize,size_t elements)
{
    double* t_host=(double*)malloc(elsize*elements);
    CHECK_CUDART(cudaMemcpy(t_host,w,elsize*elements, cudaMemcpyDeviceToHost));
    writeBufferBinary(t_host,file,elsize,elements);
    free(t_host);
}
    
void writeBufferBinaryGPUreal(double* w,const char* file,size_t elsize,size_t elements)
{
    double* t_host=(double*)malloc(elsize*elements);
    double* t_dev;
    CHECK_CUDART(cudaMalloc(&t_dev,elsize*elements));
    copyBuffer((double2*)w,(double2*)t_dev);
    fftBackward_reduced((double2*)t_dev);
    CHECK_CUDART(cudaMemcpy(t_host,t_dev,elsize*elements, cudaMemcpyDeviceToHost));
    writeBufferBinary(t_host,file,elsize,elements);
    free(t_host);
    CHECK_CUDART(cudaFree(t_dev));
}

void writeCheckpoint(vfield u, double* grid_mesh,double* time,char* fileName){
  

  hsize_t dims[3] = {NR,NT,2*NZ};
  hsize_t dimg[1] = {NR};
  hsize_t dimh[1] = {1};

  size_t size_p=NR*NT*NZ*sizeof(double2);
  
  hid_t file_id;
  herr_t status;

  double2* aux_h=(double2*)malloc(size_p);
  
  file_id = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  CHECK_CUDART(cudaMemcpy(aux_h,u.r,size_p, cudaMemcpyDeviceToHost));
  status = H5LTmake_dataset(file_id, "/ur", 3, dims, H5T_IEEE_F64LE,aux_h);

  CHECK_CUDART(cudaMemcpy(aux_h,u.t,size_p, cudaMemcpyDeviceToHost));
  status = H5LTmake_dataset(file_id, "/ut", 3, dims, H5T_IEEE_F64LE,aux_h);
  
  CHECK_CUDART(cudaMemcpy(aux_h,u.z,size_p, cudaMemcpyDeviceToHost));
  status = H5LTmake_dataset(file_id, "/uz", 3, dims, H5T_IEEE_F64LE,aux_h);
  
  
  status = H5LTmake_dataset(file_id, "/mesh", 1, dims, H5T_IEEE_F64LE,grid_mesh);
  status = H5LTmake_dataset(file_id, "/time", 1, dimh, H5T_IEEE_F64LE,time);

  status = H5Fclose(file_id);

  free(aux_h);

}



void readCheckpoint(vfield u, double* grid_mesh,double* time,char* fileName){
  

  size_t size_p=NR*NT*NZ*sizeof(double2);
  
  hid_t file_id;
  herr_t status;

  double2* aux_h=(double2*)malloc(size_p);
  
  file_id = H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);

  status = H5LTread_dataset(file_id,"/ur",H5T_IEEE_F64LE,aux_h);
  CHECK_CUDART(cudaMemcpy(u.r,aux_h,size_p, cudaMemcpyHostToDevice));

  status = H5LTread_dataset(file_id,"/ut",H5T_IEEE_F64LE,aux_h);
  CHECK_CUDART(cudaMemcpy(u.t,aux_h,size_p, cudaMemcpyHostToDevice));
  
  status = H5LTread_dataset(file_id,"/uz",H5T_IEEE_F64LE,aux_h);
  CHECK_CUDART(cudaMemcpy(u.z,aux_h,size_p, cudaMemcpyHostToDevice));
  
  
  status = H5LTread_dataset(file_id, "/mesh", H5T_IEEE_F64LE,grid_mesh);
  status = H5LTread_dataset(file_id, "/time", H5T_IEEE_F64LE,time);

  status = H5Fclose(file_id);

  free(aux_h);
  
  
}

void truncate(double2* u_trun,double2* u_pad)
{
  
  for(int i=0;i<NR;i++){
    for(int j=0;j<NT;j++){
      for(int k=0;k<NZ-1;k++){
    
        u_trun[i*NT*(NZ-1)+j*(NZ-1)+k].x=u_pad[i*NT*NZ+j*NZ+k].x;
        u_trun[i*NT*(NZ-1)+j*(NZ-1)+k].y=u_pad[i*NT*NZ+j*NZ+k].y;

      }
    }
  }
  
}

void writeFieldVis(vfield u, double* grid_mesh,double* time,char* fileName)
{
  
  int NZr=2*(NZ-1);
  
  hsize_t dims[3] = {NR,NT,NZr};
  
  hsize_t dimr[1] = {NR};
  hsize_t dimt[1] = {NT};
  hsize_t dimz[1] = {NZr};

  hsize_t dimh[1] = {1};

  
  size_t size_p=NR*NT*NZ*sizeof(double2);
  
  hid_t file_id;
  hid_t group_f,group_g,group_s,group_v;
  hid_t dataset_id, dataspace_id;
  herr_t status;

  double2* aux_h     =(double2*)malloc(size_p);
  double2* aux_h_trun=(double2*)malloc(size_p);

  
  double2* aux_d;
  CHECK_CUDART(cudaMalloc(&aux_d,size_p));

  
  file_id = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  
  group_f = H5Gcreate2(file_id, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  group_g = H5Gcreate2(file_id, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  group_s = H5Gcreate2(file_id, "/setup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  group_v = H5Gcreate2(file_id, "/fields/velocity", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  //Write data set for grid
  
  dataspace_id = H5Screate_simple(1, dimr, NULL);
  dataset_id = H5Dcreate2(file_id,"/grid/r", H5T_IEEE_F64LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     grid_mesh);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);

  double* mesh_z=(double*)malloc(NZr*sizeof(double));
  double* mesh_t=(double*)malloc(NT*sizeof(double));

  for(int i=0;i<NT;i++){
      mesh_t[i]=LT*(double)i/(double)NT;
  }

  for(int i=0;i<NZr;i++){
      mesh_z[i]=LZ*(double)i/(double)NZr;
  } 
  
  
  dataspace_id = H5Screate_simple(1, dimt, NULL);
  dataset_id = H5Dcreate2(file_id,"/grid/th", H5T_IEEE_F64LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     mesh_z);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);

  dataspace_id = H5Screate_simple(1, dimz, NULL);
  dataset_id = H5Dcreate2(file_id,"/grid/z", H5T_IEEE_F64LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     mesh_z);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  
  
  //Write data set for fields
  
  dataspace_id = H5Screate_simple(3, dims, NULL);
  dataset_id = H5Dcreate2(file_id,"/fields/pressure", H5T_IEEE_F64LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     aux_h);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  
  copyBuffer((double2*)u.r,(double2*)aux_d);
  fftBackward_reduced((double2*)aux_d);
  CHECK_CUDART(cudaMemcpy(aux_h,aux_d,size_p,cudaMemcpyDeviceToHost));
  truncate(aux_h_trun,aux_h);
  
  dataspace_id = H5Screate_simple(3, dims, NULL);
  dataset_id = H5Dcreate2(file_id,"/fields/velocity/u_r", H5T_IEEE_F64LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     aux_h_trun);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);

  copyBuffer((double2*)u.t,(double2*)aux_d);
  fftBackward_reduced((double2*)aux_d);
  CHECK_CUDART(cudaMemcpy(aux_h,aux_d,size_p, cudaMemcpyDeviceToHost));
  truncate(aux_h_trun,aux_h);

  dataspace_id = H5Screate_simple(3, dims, NULL);
  dataset_id = H5Dcreate2(file_id,"/fields/velocity/u_th", H5T_IEEE_F64LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     aux_h_trun);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  
  copyBuffer((double2*)u.z,(double2*)aux_d);
  fftBackward_reduced((double2*)aux_d);
  CHECK_CUDART(cudaMemcpy(aux_h,aux_d,size_p, cudaMemcpyDeviceToHost));
  truncate(aux_h_trun,aux_h);

  dataspace_id = H5Screate_simple(3, dims, NULL);
  dataset_id = H5Dcreate2(file_id,"/fields/velocity/u_z", H5T_IEEE_F64LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     aux_h_trun);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  

  status = H5Gclose(group_f);
  status = H5Gclose(group_g);
  status = H5Gclose(group_s);
  status = H5Gclose(group_v);
 
  //Close file
  status = H5Fclose(file_id);

  free(aux_h);
  free(aux_h_trun);
  free(mesh_z);
  free(mesh_t);
  CHECK_CUDART(cudaFree(aux_d));
   
}
