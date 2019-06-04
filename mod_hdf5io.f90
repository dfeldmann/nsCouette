!==============================================
! Collection of driver routines for writing
!  HDF5 output and a corresponding XML file
!   with xdmf metadata
!
! SVN $HeadURL: $
! SVN $Id: $
! SVN $LastChangedDate: $
!==============================================

MODULE mod_hdf5io
  
  USE HDF5
  IMPLICIT NONE

  private
  public :: init_io,finalize_io,create_file,close_file,write_hdf,write_xdmf,datafile

  INTEGER :: mpierror       ! MPI error flag
  INTEGER :: mpi_comm,mpi_info,ierr,mpi_size,mpi_rank


  type :: datafile
     integer(kind=hid_t) :: id
     integer(HID_T) :: current_group
  end type datafile

  integer :: error

  interface write_hdf
     
     module procedure write_real,write_int,write_string,write_1dreal,write_3dreal_collective

  end interface write_hdf

CONTAINS
  

  subroutine create_file(filename,file)

    type(datafile), intent(out) :: file
    character(*), intent(in) :: filename

    integer(HID_T) :: plist_id,file_id

    ! Setup file access property list with parallel I/O access.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, mpi_comm, mpi_info, error)

    ! Create the file collectively.
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)

    CALL h5pclose_f(plist_id, error)

    file%current_group=0
    file%id = file_id

  end subroutine create_file


  subroutine close_file(file)

    type(datafile), intent(in) :: file
    CALL h5fclose_f(file%id, error)

  end subroutine close_file

  subroutine init_io(comm,info,size,rank)
    integer, intent(in) :: comm,info,size,rank

    mpi_comm=comm
    mpi_info=info
    mpi_rank=rank
    mpi_size=size

     ! Initialize FORTRAN predefined datatypes
    call h5open_f(error) 
    if (error .lt. 0) call abort_io(error,'err init_io')

  end subroutine init_io


  subroutine finalize_io

    ! Close FORTRAN predefined datatypes.
    CALL h5close_f(error)
    if (error .lt. 0) call abort_io(error,'err finalize_io')
    
  end subroutine finalize_io


  subroutine abort_io(error,message)

    integer, intent(in) :: error
    character(*), intent(in) :: message

    print *,error,message

    call MPI_finalize(ierr)
    stop 'panic'

  end subroutine abort_io



  subroutine write_real(file,data,dname)

    IMPLICIT NONE

    real(KIND=8), intent(in) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file

    integer :: error

    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id


    group_id = file%current_group

    ! Create scalar dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)

    !Create the attribute
    call h5acreate_f(group_id, dname, H5T_NATIVE_DOUBLE, dspace_id, attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, (/1_hsize_t/), error)
    call h5aclose_f(attr_id, error)


    CALL h5sclose_f(dspace_id, error)


  end subroutine write_real

  subroutine write_int(file,data,dname)

    IMPLICIT NONE

    integer(KIND=4), intent(in) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file

    integer :: error

    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id


    group_id = file%current_group

    ! Create scalar dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)

    !Create the attribute
    call h5acreate_f(group_id, dname, H5T_NATIVE_INTEGER, dspace_id, attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, data, (/1_hsize_t/), error)
    call h5aclose_f(attr_id, error)


    CALL h5sclose_f(dspace_id, error)


  end subroutine write_int

  subroutine write_string(file,data,dname)

    IMPLICIT NONE

    character(*), intent(in) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file

    integer :: error

    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id,dtype_id


    group_id = file%current_group

    ! Create scalar dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)

    ! Create the datatype
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    call h5tset_size_f(dtype_id, int(len(trim(data)), kind=size_t), error)


    !Create the attribute
    call h5acreate_f(group_id, dname,dtype_id , dspace_id, attr_id, error)
    call h5awrite_f(attr_id, dtype_id, data, (/1_hsize_t/), error)
    call h5aclose_f(attr_id, error)


    CALL h5sclose_f(dspace_id, error)
    call h5tclose_f(dtype_id, error)

  end subroutine write_string

  subroutine write_1dreal(file,data,dname)

    IMPLICIT NONE

    integer, parameter :: rank=1

    real(KIND=8), intent(in) :: data(:)
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file

    integer(HSIZE_T) :: dims(rank)
    integer(HSSIZE_T), DIMENSION(rank) :: offset 
    integer :: error

    integer(HID_T) :: filespace,memspace,dset_id,plist_id,group_id

    dims = shape(data)


    group_id = file%current_group



    ! Create the data space for the dataset. 
    CALL h5screate_simple_f(rank, dims, filespace, error)

    ! Create the dataset with default properties.
    CALL h5dcreate_f(group_id, dname, H5T_NATIVE_REAL, filespace, &
         dset_id, error)

    CALL h5sclose_f(filespace, error)

    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 


    CALL h5screate_simple_f(rank, dims, memspace, error) 
    !if (rank /=0 ) call h5sselect_none_f(memspace, error)


    ! Select all.
    CALL h5dget_space_f(dset_id, filespace, error)
    CALL h5sselect_all_f (filespace,error)

    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    ! Write the dataset collectively. 
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, real(data,kind=4), dims, error, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! Close dataspaces.
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    ! Close the dataset and property list.
    CALL h5dclose_f(dset_id, error)
    CALL h5pclose_f(plist_id, error)

  end subroutine write_1dreal

  subroutine write_3dreal_collective(file,data,dimsf,dname)

    IMPLICIT NONE

    integer, parameter :: rank=3

    real(KIND=8), intent(in) :: data(:,:)
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer(HSIZE_T), dimension(rank), intent(in) :: dimsf

    integer(HSIZE_T) :: dims(rank)
    integer(HSSIZE_T), DIMENSION(rank) :: offset 
    integer :: error

    integer(HID_T) :: filespace,memspace,dset_id,plist_id,group_id

    dims(1)=dimsf(1)
    dims(2)=dimsf(2)
    dims(3)=dimsf(3)/mpi_size

    group_id = file%current_group



    ! Create the data space for the dataset. 
    CALL h5screate_simple_f(rank, dimsf, filespace, error)

    ! Create the dataset with default properties.
    CALL h5dcreate_f(group_id, dname, H5T_NATIVE_REAL, filespace, &
         dset_id, error)

    CALL h5sclose_f(filespace, error)

    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 

    offset(1) = 0
    offset(2) = 0
    offset(3) = mpi_rank * dims(3) 

    CALL h5screate_simple_f(rank, dims, memspace, error) 

    ! Select hyperslab in the file.
    CALL h5dget_space_f(dset_id, filespace, error)
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, dims, error)

    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    ! Write the dataset collectively. 
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, real(data,kind=4), dims, error, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! Close dataspaces.
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)

    ! Close the dataset and property list.
    CALL h5dclose_f(dset_id, error)
    CALL h5pclose_f(plist_id, error)


  end subroutine write_3dreal_collective



  subroutine write_xdmf(filename,time,dimsf)

    integer, parameter :: lun=42
    character(*), intent(in) :: filename
    real(kind=8), intent(in) :: time
    integer(HSIZE_T), intent(in) :: dimsf(3)

    open(lun,file=filename//'.xmf',action="write")
    write(lun,'(A)') '<?xml version="1.0" ?>'
    write(lun,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(lun,'(A)') '<Xdmf Version="2.0">'
    write(lun,'(A)') '<Domain>'
    write(lun,'(A)') '<Grid Name="mesh" GridType="Uniform">'
    write(lun,'(A,3I6,A)') '<Topology TopologyType="3DRectMesh" Dimensions="',dimsf(3:1:-1),'"/>'
    write(lun,'(A)') '<Geometry GeometryType="VXVYVZ">'

    write(lun,'(A,I6,A)') ' <DataItem Dimensions="',dimsf(1),'" Name="th" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/grid/th'
    write(lun,'(A)') ' </DataItem>'

    write(lun,'(A,I6,A)') ' <DataItem Dimensions="',dimsf(2),'" Name="z" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/grid/z'
    write(lun,'(A)') ' </DataItem>'

    write(lun,'(A,I6,A)') ' <DataItem Dimensions="',dimsf(3),'" Name="r" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/grid/r'
    write(lun,'(A)') ' </DataItem>'

    write(lun,'(A)') '</Geometry>'
    write(lun,'(A,1E11.4,A)') '<Time Value="',time,'" />'
    write(lun,'(A)') '<Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/pressure'
    write(lun,'(A)') ' </DataItem>'
    write(lun,'(A)') '</Attribute>'


#ifdef TE_CODE
    write(lun,'(A)') '<Attribute Name="temperature" AttributeType="Scalar" Center="Node">'
    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/temperature'
    write(lun,'(A)') ' </DataItem>'
    write(lun,'(A)') '</Attribute>'
#endif /* TE_CODE */

    write(lun,'(A)') '<Attribute Name="u_r" AttributeType="Scalar" Center="Node">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/u_r'
    write(lun,'(A)') '  </DataItem>'
    write(lun,'(A)') '</Attribute>'

    write(lun,'(A)') '<Attribute Name="u_th" AttributeType="Scalar" Center="Node">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/u_th'
    write(lun,'(A)') '  </DataItem>'
    write(lun,'(A)') '</Attribute>'

    write(lun,'(A)') '<Attribute Name="u_z" AttributeType="Scalar" Center="Node">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/u_z'
    write(lun,'(A)') '  </DataItem>'
    write(lun,'(A)') '</Attribute>'




    write(lun,'(A)') '<Attribute Name="velocity" AttributeType="Vector" Center="Node">'
    write(lun,'(A,4I6,A)') '<DataItem ItemType="Function" Dimensions="',dimsf(3:1:-1),3,'" Function="JOIN($0 , $1, $2)">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/u_r'
    write(lun,'(A)') '  </DataItem>'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/u_th'
    write(lun,'(A)') '  </DataItem>'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(3:1:-1),'" NumberType="Float" Precision="4" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/u_z'
    write(lun,'(A)') '  </DataItem>'

    write(lun,'(A)') '</DataItem>'
    write(lun,'(A)') '</Attribute>'


    write(lun,'(A)') '</Grid>'
    write(lun,'(A)') '</Domain>'
    write(lun,'(A)') '</Xdmf>'
    close(lun)

  end subroutine write_xdmf

END MODULE mod_hdf5io
