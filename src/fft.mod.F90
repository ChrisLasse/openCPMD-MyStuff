#include "cpmd_global.h"

MODULE fft

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_PTR, C_NULL_PTR
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE


  ! ==================================================================
  INTEGER, ALLOCATABLE :: mg(:,:)
  INTEGER, ALLOCATABLE, TARGET :: ms(:,:)
  INTEGER, ALLOCATABLE :: mz(:)

  INTEGER :: naux1,naux2

#if defined(_HAS_CUDA) || defined(_USE_SCRATCHLIBRARY)
  COMPLEX(real_8), POINTER __CONTIGUOUS, ASYNCHRONOUS :: xf(:,:)
  COMPLEX(real_8), POINTER __CONTIGUOUS, ASYNCHRONOUS :: yf(:,:)
#else
  COMPLEX(real_8), ALLOCATABLE, TARGET, ASYNCHRONOUS :: xf(:,:)
  COMPLEX(real_8), ALLOCATABLE, TARGET, ASYNCHRONOUS :: yf(:,:)
#endif

  ! ==================================================================
  INTEGER :: nr1m,nr2m,nr3m,kr1m,kr2m,kr3m,nhrm,ngrm,&
       kr2min,kr2max,kr3min,kr3max,maxrpt
  INTEGER, DIMENSION(:), ALLOCATABLE :: mxy !(maxcpu)
  ! ==================================================================
  ! ==   GATHER/SCATTER ARRAYS                                      ==
  ! ==================================================================
  INTEGER, ALLOCATABLE :: msp(:,:,:)

  ! ==================================================================
  ! NEW GENERAL PARALLEL FFT CODE
  ! ==================================================================
  INTEGER :: msrays,mfrays,llr1
  INTEGER :: qr1s,qr2s,qr3s,qr1,qr2,qr3
  INTEGER :: lr1s,lr2s,lr3s,lr1,lr2,lr3
  INTEGER :: qr2max,qr2min,qr3max,qr3min
  INTEGER :: lsrm,lfrm,lr1m,lmsq
  INTEGER :: jgw,jgws,jhg,jhgs
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: lrxpl !(0:maxcpu,2)
  INTEGER, DIMENSION(:), ALLOCATABLE :: sp5,sp8,sp9!(0:maxcpu)
  INTEGER, ALLOCATABLE :: msqs(:)
  INTEGER, ALLOCATABLE :: msqf(:)

  INTEGER, POINTER :: nzff(:)
  INTEGER, POINTER :: nzfs(:)
  INTEGER, POINTER :: inzf(:)
  INTEGER, POINTER :: inzs(:)
  INTEGER, POINTER :: inzh(:,:)

  ! ==================================================================
  ! POOL
  ! ==================================================================
  INTEGER, PARAMETER :: fftpoolsize=3
  INTEGER :: fftpool
  INTEGER :: fpoolv(28,fftpoolsize)
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: lrxpool!(0:maxcpu,2,fftpoolsize)
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: spm!(9,0:maxcpu,fftpoolsize)
  INTEGER :: lmsqmax,lnzf,lnzs
  INTEGER, ALLOCATABLE :: msqspool(:,:,:)
  INTEGER, ALLOCATABLE :: msqfpool(:,:,:)

  INTEGER, ALLOCATABLE, TARGET :: nzffp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: nzfsp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: inzfp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: inzsp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: inzhp(:,:,:)

  ! ==================================================================
  ! ==================================================================
  ! NEW SPARSE BATCH PARALLEL FFT CODE
  ! ==================================================================
  INTEGER :: a2a_msgsize, fft_batchsize, fft_numbatches, fft_residual, fft_total, fft_tune_num_it, fft_tune_max_it, fft_min_numbatches
  LOGICAL :: batch_fft
#ifdef _USE_SCRATCHLIBRARY
  COMPLEX(real_8), POINTER, SAVE __CONTIGUOUS, ASYNCHRONOUS   :: wfn_r(:,:),wfn_g(:,:)
#else
  COMPLEX(real_8), ALLOCATABLE, SAVE, TARGET, ASYNCHRONOUS  :: wfn_r(:,:),wfn_g(:,:)
#endif
  LOGICAL, ALLOCATABLE, ASYNCHRONOUS       :: locks_inv(:,:), locks_fw(:,:)
  REAL(real_8), ALLOCATABLE                :: fft_time_total(:)
  INTEGER, ALLOCATABLE                     :: fft_batchsizes(:)

  ! ==================================================================
  ! NEW GDISTRIBUTION
  ! ==================================================================
  INTEGER :: fft_numbuff
  TYPE FFT_TYPE_DESCRIPTOR

     INTEGER, ALLOCATABLE :: indx_map(:,:)
     INTEGER, ALLOCATABLE :: indx(:)
     INTEGER, ALLOCATABLE :: ind1(:)
     INTEGER, ALLOCATABLE :: ind2(:)
     INTEGER, ALLOCATABLE :: nr3_ranges(:,:)
     INTEGER, ALLOCATABLE :: nr3p(:)
     INTEGER, ALLOCATABLE :: nr3p_offset(:)
     INTEGER, ALLOCATABLE :: ir1w(:)
     INTEGER, ALLOCATABLE :: ir1p(:)
     INTEGER, ALLOCATABLE :: indw(:)
     INTEGER, ALLOCATABLE :: indp(:)
     INTEGER, ALLOCATABLE :: nsw(:)
     INTEGER, ALLOCATABLE :: nsp(:)
     INTEGER, ALLOCATABLE :: iss(:)
     INTEGER, ALLOCATABLE :: ismap(:)
     INTEGER, ALLOCATABLE :: cp_ngws(:)
     INTEGER, ALLOCATABLE :: cp_nstates(:)

     INTEGER, ALLOCATABLE :: thread_z_sticks(:,:,:,:)
     INTEGER, ALLOCATABLE :: thread_prepare_sticks(:,:)
     INTEGER, ALLOCATABLE :: thread_y_sticks(:,:)
     INTEGER, ALLOCATABLE :: thread_x_sticks(:,:)
     INTEGER, ALLOCATABLE :: thread_z_start(:,:,:,:)
     INTEGER, ALLOCATABLE :: thread_prepare_start(:,:)
     INTEGER, ALLOCATABLE :: thread_y_start(:,:)
     INTEGER, ALLOCATABLE :: thread_x_start(:,:)
     INTEGER, ALLOCATABLE :: thread_z_end(:,:,:,:)
     INTEGER, ALLOCATABLE :: thread_prepare_end(:,:)
     INTEGER, ALLOCATABLE :: thread_y_end(:,:)
     INTEGER, ALLOCATABLE :: thread_x_end(:,:)
     INTEGER, ALLOCATABLE :: thread_ngms(:)
     INTEGER, ALLOCATABLE :: thread_ngms_start(:)
     INTEGER, ALLOCATABLE :: thread_ngms_end(:)
     INTEGER, ALLOCATABLE :: cg_thread_ngms(:,:)
     INTEGER, ALLOCATABLE :: cg_thread_ngms_start(:,:)
     INTEGER, ALLOCATABLE :: cg_thread_ngms_end(:,:)
     INTEGER, ALLOCATABLE :: thread_rspace(:)
     INTEGER, ALLOCATABLE :: thread_rspace_start(:)
     INTEGER, ALLOCATABLE :: thread_rspace_end(:)

     INTEGER, ALLOCATABLE :: map_set_psi(:,:)
     INTEGER, ALLOCATABLE :: map_transpose_y2x(:,:)
     INTEGER, ALLOCATABLE :: map_transpose_x2y(:,:)
     INTEGER :: zero_transpose_y2x_start( 2 )
     INTEGER :: zero_transpose_y2x_end( 2 )
     INTEGER, ALLOCATABLE :: map_y2z(:,:)
     INTEGER, ALLOCATABLE :: map_z2y_wave(:,:)
     INTEGER, ALLOCATABLE :: map_z2y_pot(:)
     INTEGER, ALLOCATABLE :: zero_z2y_start(:,:)
     INTEGER, ALLOCATABLE :: zero_z2y_end(:,:)

     INTEGER :: nr1w
     INTEGER :: nr1p
     INTEGER :: my_nr3p
     INTEGER :: nr3px
     INTEGER :: nhg
     INTEGER :: ngw
     INTEGER :: nwst
     INTEGER :: npst
     INTEGER :: max_ngw
     INTEGER :: max_nstates

     INTEGER :: which ! 1 -> wave sticks ; 2 -> pot sticks
     INTEGER :: which_wave ! 1 -> rho ; -> 2 -> vpsi
     LOGICAL :: no_comm

     DOUBLE PRECISION :: tscale

     INTEGER :: small_chunks( 2 )
     INTEGER :: big_chunks( 2 )
     LOGICAL :: do_comm( 2 )
     INTEGER :: comm_sendrecv(2,2)
     INTEGER, ALLOCATABLE :: c2_com_num(:,:)
     INTEGER, ALLOCATABLE :: c2_com_recv(:,:)
     INTEGER, ALLOCATABLE :: s4_coms(:)

  END TYPE
  Type( FFT_TYPE_DESCRIPTOR ) :: tfft

END MODULE fft

