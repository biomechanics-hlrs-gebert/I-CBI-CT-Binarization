!------------------------------------------------------------------------------
! MODULE: raw_binary
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief:
!> Module for reading/writing raw binary files parallely.
!> PureDat provides an extended funcitonality and calls a binary blob a
!> stream file (for example .int4.st).
!
!> @description
!> There's one routine for reading raw blobs for each datatype. 
!> And a switch for big/little endian. From the viewpoint of clean code, it is 
!> not well done. But creating subroutines for preparing a read raw and 
!> releasing/finishing read raw results in three calls and an allocate in the 
!> main program. Otherwise, it results in an mpi_read_raw_prepare subroutine 
!> with a specific data type, that itself must call the actual 
!> MPI_FILE_READ_ALL routine. Since the preparation is the biggest part, the 
!> calls won't be shorter significantly.
!------------------------------------------------------------------------------
MODULE raw_binary

USE ISO_FORTRAN_ENV
USE MPI
USE global_std
USE user_interaction

IMPLICIT NONE
INTERFACE mpi_read_raw
   MODULE PROCEDURE mpi_read_raw_rk4
   MODULE PROCEDURE mpi_read_raw_rk8
   MODULE PROCEDURE mpi_read_raw_ik4
   MODULE PROCEDURE mpi_read_raw_ik2
END INTERFACE mpi_read_raw

INTERFACE mpi_write_raw
   MODULE PROCEDURE mpi_write_raw_ik2
   MODULE PROCEDURE mpi_write_raw_ik4
END INTERFACE mpi_write_raw

INTERFACE ser_read_raw
   MODULE PROCEDURE ser_read_raw_ik2
   MODULE PROCEDURE ser_read_raw_ik4
   MODULE PROCEDURE ser_read_raw_ik8
   MODULE PROCEDURE ser_read_raw_rk4
   MODULE PROCEDURE ser_read_raw_rk8
END INTERFACE ser_read_raw

INTERFACE ser_write_raw
   MODULE PROCEDURE ser_write_raw_ik2
   MODULE PROCEDURE ser_write_raw_ik4
   MODULE PROCEDURE ser_write_raw_ik8
   MODULE PROCEDURE ser_write_raw_rk4
   MODULE PROCEDURE ser_write_raw_rk8
END INTERFACE ser_write_raw

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: get_rank_section
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Get the section address/number of a specific domain
!
!> @param[in] domain No of the control volume. 
!> @param[in] sections x/y/z mesh of domains
!> @param[out] rank_section position of domain in x/y/z mesh
!------------------------------------------------------------------------------  
SUBROUTINE get_rank_section(domain, sections, rank_section)
  
INTEGER(ik), INTENT(IN)  :: domain
INTEGER(ik), DIMENSION(3), INTENT(IN)  :: sections
INTEGER(ik), DIMENSION(3), INTENT(OUT) :: rank_section

INTEGER(ik) :: rank, yrmndr, zrmndr ! remainder

!------------------------------------------------------------------------------
! Power of 2 is handled here, because with the algorithm of CASE DEFAULT, Greedy suboptimality kicks in!  
! Have a look at the corresponding Matlab/Octave testing file!
! In example at size_mpi = 128 Processors, where CASE DEFAULT will deliver 125 Processors!
!------------------------------------------------------------------------------
IF ( domain == 0_ik ) THEN
   rank_section = [ 1, 1, 1 ]
ELSE
   rank = domain + 1 ! MPI starts at 0

   !------------------------------------------------------------------------------
   ! Calculate the rank_section out of my_rank and sections [ x, y, z ]
   ! Tested via Octave. Not fully implemented by 20210503
   !------------------------------------------------------------------------------
   zrmndr = MODULO(rank, sections(1)*sections(2))
   IF (zrmndr == 0_ik) THEN
      rank_section = [ sections(1), sections(2), (rank - zrmndr) / (sections(1)*sections(2)) ]
   ELSE
      rank_section(3) = (rank - zrmndr) / (sections(1) * sections(2)) 
      yrmndr = MODULO(zrmndr, sections(1))

      IF (yrmndr == 0_ik) THEN
         rank_section = [ sections(1), (zrmndr - yrmndr) / sections(1), rank_section(3)+1 ]
      ELSE
         rank_section = [ yrmndr, (zrmndr - yrmndr) / sections(1) + 1, rank_section(3) + 1 ]
      END IF
   END IF
END IF

END SUBROUTINE get_rank_section


!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw int2 data of a binary blob. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!> @param[in] dtrep Datarepresentation 
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_ik2(filename, disp, dims, subarray_dims, subarray_origin, subarray, dtrep)

CHARACTER(*), INTENT(IN) :: filename
INTEGER(MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(INT16), DIMENSION (:,:,:), ALLOCATABLE, INTENT(INOUT) :: subarray

!------------------------------------------------------------------------------  
! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
!------------------------------------------------------------------------------  
INTEGER(mik) :: ierr, type_subarray, fh
CHARACTER(scl) :: datarep
CHARACTER(*), INTENT(IN), OPTIONAL :: dtrep

datarep = 'NATIVE'

IF(PRESENT(dtrep)) THEN
   IF (dtrep /= "") datarep = TRIM(dtrep)
END IF 

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, mik), INT(subarray_dims, mik), &
   INT(subarray_origin, mik), MPI_ORDER_FORTRAN, MPI_INTEGER2, type_subarray, ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER2, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

IF(.NOT. ALLOCATED(subarray)) ALLOCATE(subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3)))

CALL MPI_FILE_READ(fh, subarray, INT(SIZE(subarray), mik), MPI_INTEGER2, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw int4 data of a binary blob. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!> @param[in] dtrep Whether the input file is big or little endian. 
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_ik4(filename, disp, dims, subarray_dims, subarray_origin, subarray, dtrep)

CHARACTER(*), INTENT(IN) :: filename
INTEGER(MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(INT32), DIMENSION (:,:,:), ALLOCATABLE, INTENT(INOUT) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(mik) :: ierr, type_subarray, fh
CHARACTER(scl) :: datarep
CHARACTER(*), INTENT(IN), OPTIONAL :: dtrep

datarep = 'NATIVE'

IF(PRESENT(dtrep)) THEN
   IF (dtrep /= "") datarep = TRIM(dtrep)
END IF 

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, mik), INT(subarray_dims, mik), & 
   INT(subarray_origin, mik), MPI_ORDER_FORTRAN, MPI_INTEGER4, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)


CALL MPI_FILE_READ(fh, subarray, INT(SIZE(subarray), mik), MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_ik4

!------------------------------------------------------------------------------
! SUBROUTINE: uik2_to_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Convert unsigned int2 to int 2 data. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @description
!> Fortran does not know this shit. Therefore a workaround...
!
!> @param[in] subarray_in Input data
!> @param[out] subarray_out Output data
!------------------------------------------------------------------------------  
SUBROUTINE uik2_to_ik2(subarray)

INTEGER(INT16), DIMENSION (:,:,:), INTENT(INOUT) :: subarray
INTEGER(ik) :: ii, jj, kk
INTEGER(ik), DIMENSION(3) :: shp

INTEGER(INT32), DIMENSION (:,:,:), ALLOCATABLE :: temp

!------------------------------------------------------------------------------  
! Storing the array with + 65536 will cut off the image.
! At least INT32 required. All of the required variables are INT32.
!------------------------------------------------------------------------------  
shp = SHAPE(subarray)

ALLOCATE(temp(shp(1), shp(2), shp(3)))
temp = 0

DO kk=1, shp(3)
DO jj=1, shp(2)
DO ii=1, shp(1)
   IF(subarray(ii,jj,kk) .LT. 0) THEN
      temp(ii,jj,kk) = INT(subarray(ii,jj,kk), INT32) + INT(65536, INT32)
   ELSE
      temp(ii,jj,kk) = INT(subarray(ii,jj,kk), INT32)
   END IF 
END DO
END DO
END DO

subarray = INT(temp - 32768, INT16)

DEALLOCATE(temp)
END SUBROUTINE uik2_to_ik2


!------------------------------------------------------------------------------
! SUBROUTINE: uik2_to_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Convert unsigned int2 to int 2 data. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @description
!> Fortran does not know this shit. Therefore a workaround...
!
!> @param[in] subarray_in Input data
!> @param[out] subarray_out Output data
!------------------------------------------------------------------------------  
SUBROUTINE uik2_to_ik4(subarray_in, subarray_out)

INTEGER(INT16), DIMENSION (:,:,:), INTENT(IN) :: subarray_in
INTEGER(INT32), DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT) :: subarray_out
INTEGER(ik) :: ii, jj, kk
INTEGER(ik), DIMENSION(3) :: shp

!------------------------------------------------------------------------------  
! Storing the array with + 65536 will cut off the image.
! At least INT32 required. All of the required variables are INT32.
!------------------------------------------------------------------------------  
INTEGER(INT32), PARAMETER :: conv_param=0, offset=65536

shp = SHAPE(subarray_in)

ALLOCATE(subarray_out(shp(1), shp(2), shp(3)))
subarray_out = INT(0, INT32)

subarray_out = INT(subarray_in, INT32)

DO kk=1, shp(3)
DO jj=1, shp(2)
DO ii=1, shp(1)
   IF(subarray_out(ii,jj,kk) .LT. conv_param) THEN
      subarray_out(ii,jj,kk) = subarray_out(ii,jj,kk) + offset
   END IF 
END DO
END DO
END DO


END SUBROUTINE uik2_to_ik4

!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_rk4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw real 4 (single precision) data of a binary blob. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!> @param[in] dtrep Whether the input file is big or little endian. 
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_rk4(filename, disp, dims, subarray_dims, subarray_origin, subarray, dtrep)

CHARACTER(*), INTENT(IN) :: filename
INTEGER(MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
REAL(REAL32), DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT) :: subarray
CHARACTER(*), INTENT(IN), OPTIONAL :: dtrep

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(mik) :: ierr, type_subarray, fh
CHARACTER(scl) :: datarep

datarep = 'NATIVE'

IF(PRESENT(dtrep)) THEN
   IF (dtrep /= "") datarep = TRIM(dtrep)
END IF 

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  
CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, mik), INT(subarray_dims, mik), & 
   INT(subarray_origin, mik), MPI_ORDER_FORTRAN, MPI_REAL, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_REAL, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

ALLOCATE(subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3)))

CALL MPI_FILE_READ_ALL(fh, subarray, INT(SIZE(subarray), mik), MPI_REAL, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_rk4

!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_rk8
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw real 8 (double precision) data of a binary blob. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @param[in] fh File handle
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!> @param[in] dtrep Whether the input file is big or little endian. 
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_rk8(filename, disp, dims, subarray_dims, subarray_origin, subarray, dtrep)

CHARACTER(*), INTENT(IN) :: filename
INTEGER(MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
REAL(REAL64), DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(mik) :: ierr, type_subarray, fh
CHARACTER(scl) :: datarep
CHARACTER(*), INTENT(IN), OPTIONAL :: dtrep

datarep = 'NATIVE'

IF(PRESENT(dtrep)) THEN
   IF (dtrep /= "") datarep = TRIM(dtrep)
END IF 

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, mik), INT(subarray_dims, mik), & 
   INT(subarray_origin, mik), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_DOUBLE_PRECISION, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

ALLOCATE(subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3)))

CALL MPI_FILE_READ_ALL(fh, subarray, INT(SIZE(subarray), mik), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_rk8


!------------------------------------------------------------------------------
! SUBROUTINE: mpi_write_raw_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Write raw binary data. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @param[in] fh File handle
!> @param[in] disp Length of the header (bytes) - position to write to
!> @param[in] filename File name
!> @param[in] dims Voxels per direction
!> @param[in] subarray_dims Voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the subarray
!> @param[in] subarray Scalar field / Image data
!> @param[in] dtrep Whether the input file is big or little endian. 
!------------------------------------------------------------------------------  
 SUBROUTINE mpi_write_raw_ik2 (filename, disp, dims, subarray_dims, subarray_origin, subarray, dtrep)
! type = 'int2', 'int4'
! IF type = uint2 - send an int4 and let it convert into int2 (!) Have a look at the src for details

CHARACTER(*), INTENT(IN) :: filename
INTEGER(MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(INT16), DIMENSION (:,:,:), INTENT(IN) :: subarray
CHARACTER(*), INTENT(IN), OPTIONAL :: dtrep

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(mik)  :: fh, ierr, type_subarray
CHARACTER(scl) :: datarep = 'NATIVE'

datarep = 'NATIVE'

IF(PRESENT(dtrep)) THEN
   IF (dtrep /= "") datarep = TRIM(dtrep)
END IF 

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, mik), INT(subarray_dims, mik), & 
   INT(subarray_origin, mik), MPI_ORDER_FORTRAN, MPI_INTEGER2, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER2, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

CALL MPI_FILE_WRITE(fh, subarray, INT(SIZE(subarray), mik), MPI_INTEGER2, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_write_raw_ik2


!------------------------------------------------------------------------------
! SUBROUTINE: mpi_write_raw_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Write raw binary data. @Datarepresenation: 
!> NATIVE=LittleEndian on file system, EXTERNAL32 --> BigEndian. 
!> Please check and test if you need this feature!! Depends on Hardware.
!
!> @description
!
!> @param[in] disp Length of the header (bytes) - position to write to
!> @param[in] filename File name
!> @param[in] dims Voxels per direction
!> @param[in] subarray_dims Voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the subarray
!> @param[in] subarray Scalar field / Image data
!> @param[in] dtrep Whether the input file is big or little endian. 
!------------------------------------------------------------------------------  
 SUBROUTINE mpi_write_raw_ik4 (filename, disp, dims, subarray_dims, subarray_origin, subarray, dtrep)
! type = 'int2', 'int4'
! IF type = uint2 - send an int4 and let it convert into int2 (!) Have a look at the src for details

CHARACTER(*), INTENT(IN) :: filename
INTEGER(MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(INT32), DIMENSION (:,:,:), INTENT(IN) :: subarray
CHARACTER(*), INTENT(IN), OPTIONAL :: dtrep

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(mik)  :: fh, ierr, type_subarray
CHARACTER(scl) :: datarep = 'NATIVE'

datarep = 'NATIVE'

IF(PRESENT(dtrep)) THEN
   IF (dtrep /= "") datarep = TRIM(dtrep)
END IF 

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), &
   MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, mik), INT(subarray_dims, mik), & 
   INT(subarray_origin, mik), MPI_ORDER_FORTRAN, MPI_INTEGER4, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, type_subarray, &
   TRIM(datarep), MPI_INFO_NULL, ierr)

CALL MPI_FILE_WRITE(fh, subarray, INT(SIZE(subarray), mik), &
   MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)
CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_write_raw_ik4


!------------------------------------------------------------------------------
! SUBROUTINE: ser_write_raw_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[in] array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_write_raw_ik2(fh, filename, array, representation)

INTEGER(INT16), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/ser_write_raw.aux.f90"

END SUBROUTINE ser_write_raw_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: ser_write_raw_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[in] array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_write_raw_ik4(fh, filename, array, representation)

INTEGER(INT32), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/ser_write_raw.aux.f90"

END SUBROUTINE ser_write_raw_ik4

!------------------------------------------------------------------------------
! SUBROUTINE: ser_write_raw_ik8
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[in] array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_write_raw_ik8(fh, filename, array, representation)

INTEGER(INT64), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/ser_write_raw.aux.f90"

END SUBROUTINE ser_write_raw_ik8

!------------------------------------------------------------------------------
! SUBROUTINE: ser_write_raw_rk4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[in] array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_write_raw_rk4(fh, filename, array, representation)

REAL(REAL32), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/ser_write_raw.aux.f90"

END SUBROUTINE ser_write_raw_rk4

!------------------------------------------------------------------------------
! SUBROUTINE: ser_write_raw_rk8
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[in] array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_write_raw_rk8(fh, filename, array, representation)

REAL(REAL64), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/ser_write_raw.aux.f90"

END SUBROUTINE ser_write_raw_rk8

!------------------------------------------------------------------------------
! SUBROUTINE: ser_read_raw_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Read raw binary data serially. Swap endianness if necessary.
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[out] Array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_read_raw_ik2(fh, filename, array, representation)

INTEGER(INT16), DIMENSION(:,:,:), INTENT(OUT) :: array
INCLUDE "./include_f90/ser_read_raw.aux.f90"

END SUBROUTINE ser_read_raw_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: ser_read_raw_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Read raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[out] Array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_read_raw_ik4(fh, filename, array, representation)

INTEGER(INT32), DIMENSION(:,:,:), INTENT(OUT) :: array
INCLUDE "./include_f90/ser_read_raw.aux.f90"

END SUBROUTINE ser_read_raw_ik4

!------------------------------------------------------------------------------
! SUBROUTINE: ser_read_raw_ik8
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Read raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[out] Array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_read_raw_ik8(fh, filename, array, representation)

INTEGER(INT64), DIMENSION(:,:,:), INTENT(OUT) :: array
INCLUDE "./include_f90/ser_read_raw.aux.f90"

END SUBROUTINE ser_read_raw_ik8

!------------------------------------------------------------------------------
! SUBROUTINE: ser_read_raw_rk4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Read raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[out] Array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_read_raw_rk4(fh, filename, array, representation)

REAL(REAL32), DIMENSION(:,:,:), INTENT(OUT) :: array
INCLUDE "./include_f90/ser_read_raw.aux.f90"

END SUBROUTINE ser_read_raw_rk4

!------------------------------------------------------------------------------
! SUBROUTINE: ser_read_raw_rk8
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Read raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[out] Array Raw data
!> @param[in] representation Optional swap of endianness
!------------------------------------------------------------------------------
SUBROUTINE ser_read_raw_rk8(fh, filename, array, representation)

REAL(REAL64), DIMENSION(:,:,:), INTENT(OUT) :: array
INCLUDE "./include_f90/ser_read_raw.aux.f90"

END SUBROUTINE ser_read_raw_rk8

END MODULE raw_binary

!------------------------------------------------------------------------------
! MODULE: vtk_meta_data
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief:
!> Module for reading/writing vtk-structured-point files
!------------------------------------------------------------------------------
MODULE vtk_meta_data

USE ISO_FORTRAN_ENV
USE global_std
USE strings
USE user_interaction
USE raw_binary

IMPLICIT NONE

INTERFACE write_ser_vtk
   MODULE PROCEDURE write_ser_vtk_ik2
   MODULE PROCEDURE write_ser_vtk_ik4
   MODULE PROCEDURE write_ser_vtk_rk8
END INTERFACE write_ser_vtk

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: write_vtk_struct_points_header
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write a *.vtk structured points header
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[in] type Data type
!> @param[in] spcng Distance between to voxels/scalars
!> @param[in] origin physical origin (mm)
!> @param[in] dims Amount of voxels/scalars per direction
!------------------------------------------------------------------------------
SUBROUTINE write_vtk_struct_points_header (fh, filename, type, spcng, origin, dims)

! It's HIGHLY recommended to check the existence of the output file prior to CALLing this
! Subroutine! Otherwise the program will crash. It's not double-checkd here, because this
! sequence often is placed at the very end of a program, which may run some time.

INTEGER  (ik), INTENT(IN) :: fh
CHARACTER(*) , INTENT(IN) :: filename
CHARACTER(*) , INTENT(IN) :: type
INTEGER(ik), INTENT(IN) :: dims(3)
REAL(rk), INTENT(IN) :: spcng(3), origin(3)

CHARACTER(scl) :: datatype=''
LOGICAL :: exist

INQUIRE(UNIT=fh, exist=exist)

IF(.NOT. exist) THEN
   OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='OLD', POSITION='APPEND')
END IF

SELECT CASE(TRIM(ADJUSTL(type)))
   CASE('uik2'); datatype = "unsigned_short"
   CASE('ik2'); datatype = "short"
   CASE('ik4'); datatype = "int"
   CASE('rk4'); datatype = "float"
   CASE('rk8'); datatype = "double"
   CASE DEFAULT
      CALL print_err_stop(fh, "No valid datatype given.", 1)
END SELECT

OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='NEW')

WRITE(fh,'(A)')          "# vtk DataFile Version 4.2" ! Compatibility issue
WRITE(fh,'(A)')          "vtk output"
WRITE(fh,'(A)')          "BINARY"
WRITE(fh,'(A)')          "DATASET STRUCTURED_POINTS"
WRITE(fh,'(A,3(I5))')    "DIMENSIONS", dims
WRITE(fh,'(A,3(F11.6))') "SPACING ", spcng
WRITE(fh,'(A,3(F11.6))') "ORIGIN ", origin
WRITE(fh,'(A, I0)')      "POINT_DATA ", PRODUCT(INT(dims, INT64))
WRITE(fh,'(A)')          "SCALARS DICOMImage "//TRIM(ADJUSTL(datatype))
WRITE(fh,'(A)')          "LOOKUP_TABLE default"

CLOSE(UNIT=fh)
END SUBROUTINE write_vtk_struct_points_header

!------------------------------------------------------------------------------
! SUBROUTINE: write_vtk_struct_points_footer
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write a *.vtk structured points footer
!
!> @param[in] fh File handle
!> @param[in] filename File name
!------------------------------------------------------------------------------
SUBROUTINE write_vtk_struct_points_footer (fh, filename)

INTEGER(ik), INTENT(IN) :: fh
CHARACTER(*) :: filename

LOGICAL :: opened

INQUIRE(UNIT=fh, opened=opened)

IF(.NOT. opened) THEN
   OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='OLD', POSITION='APPEND')
END IF

WRITE(fh ,'(A)')''
WRITE(fh ,'(A)')"METADATA"
WRITE(fh ,'(A)')"INFORMATION 0"
WRITE(fh ,'(A)')

CLOSE(UNIT=fh)
END SUBROUTINE write_vtk_struct_points_footer


!------------------------------------------------------------------------------
! SUBROUTINE: read_vtk_meta
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Read a *.vtk structured points header (footer not relevant)
!
!> @param[in] filename File name
!> @param[out] disp Displacement, length of header (bytes)
!> @param[out] dims Voxels per direction
!> @param[out] origin Physical (mm) origin
!> @param[out] spcng Physical distance between two voxels (mm)
!> @param[out] type Data type contained in binary blob
!------------------------------------------------------------------------------
SUBROUTINE read_vtk_meta(filename, disp, dims, origin, spcng, type)

CHARACTER(*), INTENT(IN) :: filename
INTEGER(ik), INTENT(OUT) :: disp
INTEGER(ik), INTENT(OUT) :: dims(3)
REAL   (rk), INTENT(OUT) :: spcng(3), origin(3)
CHARACTER(*), INTENT(OUT) :: type

!-- Initialize variables in case they're not used
INTEGER  (ik) :: ii=0, ntokens, fh

CHARACTER(mcl) :: line
CHARACTER(mcl) :: tokens(100)
CHARACTER(mcl), DIMENSION(3) :: token

!------------------------------------------------------------------------------
! Determine a new unit
!------------------------------------------------------------------------------  
fh = give_new_unit()
OPEN(UNIT=fh, FILE=TRIM(filename), STATUS="OLD")

disp=0

DO ii=1,10
   READ(fh,'(A)') line
   disp=disp+LEN(TRIM(line))+1_ik ! eol characters, white charcter
   
   CALL parse(str=line,delims=" ",args=tokens,nargs=ntokens)

   IF (ntokens > 0) THEN
   
      SELECT CASE(tokens(1))
         CASE('DIMENSIONS'); READ(tokens(2:4),'(I15)') dims(1:3)
         CASE('SPACING'); READ(tokens(2:4),'(F15.6)') spcng(1:3)  
         CASE('ORIGIN'); READ(tokens(2:4),'(F15.6)') origin(1:3)  
         CASE('DATASET')
            IF (tokens(2) /= "STRUCTURED_POINTS") THEN
               mssg = "The input file "//TRIM(filename)//" does not contain STRUCTURED_POINTS!"
               CALL print_err_stop(std_out, mssg, 1)
            END IF

         CASE('SCALARS')
            token(3) = tokens(3)

            SELECT CASE( TRIM( token(3) ) )
               CASE('float') ; type = 'rk4'
               CASE('double'); type = 'rk8'
               CASE('int')   ; type = 'ik4'
               CASE('short') ; type = 'ik2'
               CASE('unsigned_short'); type = 'uik2'
               CASE DEFAULT
                  WRITE(*,'(A)') "No valid type given in *.vtk File." 
            END SELECT
      END SELECT
   END IF !ntokens <0
END DO

CLOSE(fh)

END SUBROUTINE read_vtk_meta

!------------------------------------------------------------------------------
! SUBROUTINE: write_ser_vtk_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Basically a simple wrapper
!
!> @param[in] filename File name
!> @param[in] type Data type contained in binary blob
!> @param[in] spcng Physical distance between two voxels (mm)
!> @param[in] dims Voxels per direction
!> @param[in] origin Physical (mm) origin
!> @param[in] array Data
!------------------------------------------------------------------------------
SUBROUTINE write_ser_vtk_ik2(filename, type, spcng, dims, origin, array)

INTEGER(INT16), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/write_ser_vtk.aux.f90"

END SUBROUTINE write_ser_vtk_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: write_ser_vtk_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Basically a simple wrapper
!
!> @param[in] filename File name
!> @param[in] type Data type contained in binary blob
!> @param[in] spcng Physical distance between two voxels (mm)
!> @param[in] dims Voxels per direction
!> @param[in] origin Physical (mm) origin
!> @param[in] array Data
!------------------------------------------------------------------------------
SUBROUTINE write_ser_vtk_ik4(filename, type, spcng, dims, origin, array)

INTEGER(INT32), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/write_ser_vtk.aux.f90"

END SUBROUTINE write_ser_vtk_ik4

!------------------------------------------------------------------------------
! SUBROUTINE: write_ser_vtk_rk8
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Basically a simple wrapper
!
!> @param[in] filename File name
!> @param[in] type Data type contained in binary blob
!> @param[in] spcng Physical distance between two voxels (mm)
!> @param[in] dims Voxels per direction
!> @param[in] origin Physical (mm) origin
!> @param[in] array Data
!------------------------------------------------------------------------------
SUBROUTINE write_ser_vtk_rk8(filename, type, spcng, dims, origin, array)

REAL(REAL64), DIMENSION(:,:,:), INTENT(IN) :: array
INCLUDE "./include_f90/write_ser_vtk.aux.f90"

END SUBROUTINE write_ser_vtk_rk8

END MODULE vtk_meta_data
