!------------------------------------------------------------------------------
! MODULE: binarize
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @description: 
!> Routines to binarze a scalar field/a computed tomography image.
!> Ideally, the images already are filtered.
!------------------------------------------------------------------------------
MODULE binarize

USE ISO_FORTRAN_ENV
USE global_std

IMPLICIT NONE

INTERFACE ct_binarize
    MODULE PROCEDURE ct_binarize_ik2
    MODULE PROCEDURE ct_binarize_ik4 
END INTERFACE ct_binarize

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: ct_binarize_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Integer 2 routine to binarze a scalar field/a computed tomography image.
!
!> @param[inout] array Data to binarize
!> @param[in] in_lo_hi Input values of the CT image
!> @param[in] out_lo_hi Final values of binarized image
!------------------------------------------------------------------------------
SUBROUTINE ct_binarize_ik2(array, in_lo_hi, out_lo_hi)

INTEGER(KIND=INT16), DIMENSION(:,:,:), INTENT(INOUT) :: array
INTEGER(KIND=ik), DIMENSION(2), INTENT(IN) :: in_lo_hi, out_lo_hi

INTEGER(KIND=ik) :: ii, jj, kk

DO kk=1, SIZE(array, DIM=3)
DO jj=1, SIZE(array, DIM=2)
DO ii=1, SIZE(array, DIM=1)
    IF ((array(ii,jj,kk) >= in_lo_hi(1)) .AND. &
        (array(ii,jj,kk) <= in_lo_hi(2))) THEN
        array(ii,jj,kk) = out_lo_hi(1)
    ELSE
        array(ii,jj,kk) = out_lo_hi(2)
    END IF
END DO
END DO
END DO      
END SUBROUTINE ct_binarize_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: ct_binarize_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Integer 2 routine to binarze a scalar field/a computed tomography image.
!
!> @param[inout] array Data to binarize
!> @param[in] in_lo_hi Input values of the CT image
!> @param[in] out_lo_hi Final values of binarized image
!------------------------------------------------------------------------------
SUBROUTINE ct_binarize_ik4(array, in_lo_hi, out_lo_hi)

INTEGER(KIND=INT32), DIMENSION(:,:,:), INTENT(INOUT) :: array
INTEGER(KIND=ik), DIMENSION(2), INTENT(IN) :: in_lo_hi, out_lo_hi

INTEGER(KIND=ik) :: ii, jj, kk

DO kk=1, SIZE(array, DIM=3)
DO jj=1, SIZE(array, DIM=2)
DO ii=1, SIZE(array, DIM=1)
    IF ((array(ii,jj,kk) >= in_lo_hi(1)) .AND. &
        (array(ii,jj,kk) <= in_lo_hi(2))) THEN
        array(ii,jj,kk) = out_lo_hi(1)
    ELSE
        array(ii,jj,kk) = out_lo_hi(2)
    END IF
END DO
END DO
END DO      
END SUBROUTINE ct_binarize_ik4

END MODULE binarize


!>-------------------------------------------------------------------
!> Parallel Binarization of CT Images
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    25.05.2021
!> LastMod: 02.01.2022
!>-------------------------------------------------------------------
PROGRAM ctbinarization

USE ISO_FORTRAN_ENV
USE global_std
USE user_interaction
USE meta
USE MPI
USE raw_binary
USE binarize

IMPLICIT NONE

! Parameter
INTEGER(KIND=ik), PARAMETER :: debug = 2   ! Choose an even integer!!

CHARACTER(LEN=mcl), DIMENSION(:), ALLOCATABLE :: m_rry
CHARACTER(LEN=scl) :: type, binary, invert, restart, restart_cmd_arg
CHARACTER(LEN=  8) :: date
CHARACTER(LEN= 10) :: time

INTEGER(KIND=INT16), DIMENSION(:,:,:), ALLOCATABLE :: rry_ik2
INTEGER(KIND=INT32), DIMENSION(:,:,:), ALLOCATABLE :: rry_ik4
INTEGER(KIND=ik), DIMENSION(3) :: dims, rry_dims, sections, rank_section
INTEGER(KIND=ik), DIMENSION(3) :: remainder_per_dir, dims_reduced, subarray_origin
INTEGER(KIND=ik), DIMENSION(2) :: array_lo_hi
INTEGER(KIND=ik) :: in_lo, in_hi, img_max

REAL(KIND=rk) :: start, end

LOGICAL :: stp

! MPI variables
INTEGER(KIND=mik) :: ierr, my_rank, size_mpi

!------------------------------------------------------------------------------
! Invoke MPI 
!------------------------------------------------------------------------------
CALL mpi_init(ierr)
CALL print_err_stop(std_out, "MPI_INIT didn't succeed", INT(ierr, KIND=ik))

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't be retrieved", INT(ierr, KIND=ik))

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't be retrieved", INT(ierr, KIND=ik))

IF (size_mpi < 2) CALL print_err_stop(std_out, "At least two ranks required to execute this program.", 1)

!------------------------------------------------------------------------------
! Rank 0 -- Init (Master) Process and broadcast init parameters 
!------------------------------------------------------------------------------
IF (my_rank==0) THEN
    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    std_out = determine_stout()

    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL show_title()
 
    IF(debug >=0) WRITE(*,FMT_MSG) "Post mortem info probably in ./datasets/.temporary.std_out"

    CALL CPU_TIME(start)

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, stp, restart, restart_cmd_arg)
    IF(stp) GOTO 1001

    !------------------------------------------------------------------------------
    ! Check and open the input file; Modify the Meta-Filename / Basename
    ! Define the new application name first
    !------------------------------------------------------------------------------
    global_meta_prgrm_mstr_app = 'CTBI' 
    global_meta_program_keyword = 'CT_BINARIZATION'
    CALL meta_append(m_rry)

    !------------------------------------------------------------------------------
    ! Parse input
    !------------------------------------------------------------------------------
    WRITE(std_out, FMT_TXT) 'Reading data from *.meta file.'

    CALL meta_read(std_out, 'TYPE_RAW'  , m_rry, type)
    CALL meta_read(std_out, 'DIMENSIONS', m_rry, dims)

    CALL meta_read(std_out, 'BINARIZE_HI', m_rry, img_max)
    CALL meta_read(std_out, 'BIN_INVERT' , m_rry, invert)
    CALL meta_read(std_out, 'HU_THRSH_HI', m_rry, in_hi)
    CALL meta_read(std_out, 'HU_THRSH_LO', m_rry, in_lo)

    IF((type /= "ik2") .AND. (type /= "ik4")) THEN
        mssg = "Program only supports ik2 and ik4 for 'TYPE_RAW'"
        CALL print_err_stop(std_out, mssg, 1)
    END IF

END IF ! my_rank==0

!------------------------------------------------------------------------------
! Send required variables
!------------------------------------------------------------------------------
CALL MPI_BCAST(in%p_n_bsnm , INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(type  , INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(invert, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(img_max, 1_mik, MPI_INTEGER, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(dims   , 3_mik, MPI_INTEGER, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! Get dimensions for each domain. Every processor reveives its own domain.
! Therefore, each my_rank calculates its own address/dimensions/parameters.
!
! The decomposition of an image depends on the task to be done. Therefore, 
! these calculations are not encapsuled in a routine. However, this may
! follow in an update during refactoring the 3D scalar image filter.
!
! Allocation of subarray memory is done in the read_raw routines.
!------------------------------------------------------------------------------
sections=0
CALL MPI_DIMS_CREATE (size_mpi, 3_mik, sections, ierr)
CALL get_rank_section(my_rank, sections, rank_section)

remainder_per_dir = MODULO(dims, sections)

dims_reduced   = dims - remainder_per_dir

rry_dims  = (dims_reduced / sections)

subarray_origin = (rank_section-1_ik) * (rry_dims)

! Add the remainder to the last domains of each dimension
IF(rank_section(1) == sections(1)) rry_dims(1) = rry_dims(1) + remainder_per_dir(1)
IF(rank_section(2) == sections(2)) rry_dims(2) = rry_dims(2) + remainder_per_dir(2)
IF(rank_section(3) == sections(3)) rry_dims(3) = rry_dims(3) + remainder_per_dir(3)

IF(my_rank == 0) THEN
    WRITE(std_out, FMT_TXT) 'Reading binary information of *.raw file.'

    ! DEBUG INFORMATION
    IF (debug >= 0) THEN 
        CALL DATE_AND_TIME(date, time)
        
        WRITE(std_out, FMT_TXT) "Date: "//date//" [ccyymmdd]"
        WRITE(std_out, FMT_TXT) "Time: "//time//" [hhmmss.sss]"  
        WRITE(std_out, FMT_TXT_SEP)
        WRITE(std_out, FMT_MSG_AI0) "Debug Level:", debug
        WRITE(std_out, FMT_MSG_AI0) "Processors:", size_mpi  
        WRITE(std_out, FMT_MSG) "Calculation of domain sectioning:"
        WRITE(std_out, FMT_MSG)
        WRITE(std_out, FMT_MSG_A3I0) "sections: ", sections
        WRITE(std_out, FMT_MSG_A3I0) "dims: ", dims
        WRITE(std_out, FMT_MSG_A3I0) "dims_reduced: ", dims_reduced
        WRITE(std_out, FMT_MSG_A3I0) "subarray_origin: ", subarray_origin
        WRITE(std_out, FMT_MSG_SEP)
        WRITE(std_out, FMT_MSG)     "Binarization:"
        WRITE(std_out, FMT_MSG_AI0) "Threshold lo: ", in_lo
        WRITE(std_out, FMT_MSG_AI0) "Threshold hi: ", in_hi
        WRITE(std_out, FMT_MSG_AI0) "Effective array low: ", array_lo_hi(1)
        WRITE(std_out, FMT_MSG_AI0) "Effective array high: ", array_lo_hi(2)
        WRITE(std_out, FMT_MSG_SEP)
        
        FLUSH(std_out)
    END IF
END IF

!------------------------------------------------------------------------------
! Read binary part of the vtk file - basically a *.raw file
!------------------------------------------------------------------------------
SELECT CASE(type)
    CASE('ik2') 
        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//raw_suf, 0_8, dims, rry_dims, subarray_origin, rry_ik2)
    CASE('ik4') 
        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//raw_suf, 0_8, dims, rry_dims, subarray_origin, rry_ik4)
END SELECT

!------------------------------------------------------------------------------
! Initialize array low/high and invert image if requested.
!------------------------------------------------------------------------------
array_lo_hi = [ 0_ik, img_max ]
IF(invert == 'YES') array_lo_hi = [ img_max, 0_ik ]

!------------------------------------------------------------------------------
! Compute binarization
!------------------------------------------------------------------------------
IF(my_rank==0) WRITE(std_out, FMT_TXT) 'Binarizing image.'


SELECT CASE(type)
    CASE('ik2'); CALL ct_binarize(rry_ik2, [ in_lo, in_hi ], array_lo_hi)
    CASE('ik4'); CALL ct_binarize(rry_ik4, [ in_lo, in_hi ], array_lo_hi)
END SELECT

!------------------------------------------------------------------------------
! Write raw data
!------------------------------------------------------------------------------
IF(my_rank==0) WRITE(std_out, FMT_TXT) 'Writing binary information to *.raw file.'

SELECT CASE(type)
    CASE('ik2') 
        CALL mpi_write_raw(TRIM(out%p_n_bsnm)//raw_suf, 0_8, dims, rry_dims, subarray_origin, rry_ik2)
    CASE('ik4') 
        CALL mpi_write_raw(TRIM(out%p_n_bsnm)//raw_suf, 0_8, dims, rry_dims, subarray_origin, rry_ik4)
END SELECT


!------------------------------------------------------------------------------
! Jump to end for a more gracefully ending of the program in specific cases :-)
!------------------------------------------------------------------------------
1001 CONTINUE

!------------------------------------------------------------------------------
! Finish program
!------------------------------------------------------------------------------
IF(my_rank == 0) THEN
    CALL CPU_TIME(end)

    WRITE(std_out, FMT_TXT_AF0A) 'Finishing the program took', end-start,'seconds.'
    WRITE(std_out, FMT_TXT_SEP)

    CALL meta_signing(binary)
    CALL meta_close()

    IF (std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END IF ! (my_rank == 0)

Call MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", INT(ierr, KIND=ik))

END PROGRAM ctbinarization
