#include "cpmd_global.h"
#if defined(__FFT_HASNT_THREADED_COPIES)
#define HASNT_THREADED_COPIES .TRUE.
#else
#define HASNT_THREADED_COPIES .FALSE.
#endif

#if defined(__FFT_HAS_LOW_LEVEL_TIMERS)
#define HAS_LOW_LEVEL_TIMERS .TRUE.
#else
#define HAS_LOW_LEVEL_TIMERS .FALSE.
#endif

#if defined(__FFT_HAS_SPECIAL_COPY)
#define HAS_SPECIAL_COPY __FFT_HAS_SPECIAL_COPY
#else
#define HAS_SPECIAL_COPY 0
#endif

#ifndef __COLLAPSE2
#if defined(__FFT_HAS_OMP_COLLAPSE)
#define __COLLAPSE2 collapse(2)
#else
#define __COLLAPSE2
#endif
#endif

MODULE fftutil_utils
  USE cppt,                            ONLY: nzh_r,&
                                             indz_r
  USE fft,                             ONLY: FFT_TYPE_DESCRIPTOR,&
                                             fft_batchsize,&
                                             fft_numbuff
  USE fft_maxfft,                      ONLY: maxfft
  USE fftkernal_utils,                 ONLY: z2y_fillsend,&
                                             z2y_fillrecv,&
                                             y2z_fillsend,&
                                             y2z_fillrecv,&
                                             y2z_voidrecv
  USE fftnew_utils,                    ONLY: locks_omp,&
                                             locks_omp_big,&
                                             locks_calc_1,&
                                             locks_calc_2
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE mltfft_utils,                    ONLY: mltfft_fftw_threadsafe
  USE mp_interface,                    ONLY: mp_all2all,&
                                             mp_startall,&
                                             mp_waitall
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: type_cast,&
                                             reshape_inplace
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             parap,&
                                             spar,&
                                             parm,&
                                             cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr_no_omp,&
                                             zsctr_no_omp
  USE zeroing_utils,                   ONLY: zeroing

  USE, INTRINSIC :: iso_fortran_env

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phase
  PUBLIC :: putz
  PUBLIC :: getz
  PUBLIC :: unpack_y2x
  PUBLIC :: pack_y2x
  PUBLIC :: fft_comm
  PUBLIC :: pack_x2y
  PUBLIC :: unpack_x2y
  PUBLIC :: phasen
!TK special routines for batched fft
  PUBLIC :: getz_n
  PUBLIC :: putz_n
  PUBLIC :: pack_x2y_n
  PUBLIC :: unpack_x2y_n
  PUBLIC :: pack_y2x_n
  PUBLIC :: unpack_y2x_n
!TK
!CLR special routines for new_gdist FFT
  PUBLIC :: set_psi_new_gdist
  PUBLIC :: fft_comm_preinitialized
  PUBLIC :: fft_comm_ALL2ALL
  PUBLIC :: invfft_z_section
  PUBLIC :: invfft_y_section
  PUBLIC :: invfft_x_section
  PUBLIC :: fwfft_x_section
  PUBLIC :: fwfft_y_section
  PUBLIC :: fwfft_z_section
  REAL(real_8), PARAMETER :: scal = 1.0
!CLR
CONTAINS


  ! ==================================================================
  SUBROUTINE phase(f)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: f(fpar%kr1,fpar%kr2s,fpar%kr3s)

    INTEGER                                  :: ii, ijk1, ijk2, isub, j, k

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
    ijk2=parap%nrxpl(parai%mepos,2)-parap%nrxpl(parai%mepos,1)+1 ! jh-mb
    !$omp parallel do private (J,K,II,IJK1) shared(IJK2) __COLLAPSE2
    DO k=1,spar%nr3s
       DO j=1,spar%nr2s
          ii=k+j+parap%nrxpl(parai%mepos,1)
          ijk1=MOD(ii,2)+1
          DO ii=ijk1,ijk2,2
             f(ii,j,k)=-f(ii,j,k)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phase
  ! ==================================================================
  SUBROUTINE putz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: a(krmax-krmin+1,*)
    INTEGER                                  :: kr, m
    COMPLEX(real_8)                          :: b(kr,m)

    CHARACTER(*), PARAMETER                  :: procedureN = 'putz'

    INTEGER                                  :: isub, n

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL zeroing(b)!,m*kr)
    CALL matmov(n,m,a,n,b(krmin,1),kr)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz
  ! ==================================================================
  SUBROUTINE getz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: b(krmax-krmin+1,*)
    INTEGER                                  :: kr
    COMPLEX(real_8)                          :: a(kr,*)
    INTEGER                                  :: m

    CHARACTER(*), PARAMETER                  :: procedureN = 'getz'

    INTEGER                                  :: isub, n

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL matmov(n,m,a(krmin,1),kr,b,n)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE getz
  ! ==================================================================
  ! ==================================================================
  ! CODE FOR NEW FFT ROUTINES
  ! ==================================================================


  SUBROUTINE unpack_y2x(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_Y2X'
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ..Pack the data for sending
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy_s(2*nrx,yf4(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf4(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,yf(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_y2x
  ! ==================================================================
  SUBROUTINE pack_y2x(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft, mproc, &
                                                sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, &
                                                mxrp, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_Y2X'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Pack the data for sending
    ! IF(HAS_LOW_LEVEL_TIMERS) CALL TISET(procedureN,ISUB1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL cgthr_z(mxrp,yf(jj),xf4(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf4(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default

          !$omp parallel do shared(SP8,NRX,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zgthr_no_omp(mxrp,yf(jj),xf(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_y2x
  ! ==================================================================
  SUBROUTINE fft_comm(xf,yf,lda,tr4a2a, comm )
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8), TARGET                  :: xf(*), yf(*)
    INTEGER, INTENT(IN)                      :: lda
    LOGICAL, INTENT(IN)                      :: tr4a2a
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)                      :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif

    CHARACTER(*), PARAMETER                  :: procedureN = 'fft_comm'

    COMPLEX(real_4), POINTER                 :: xf4(:), yf4(:)
    INTEGER                                  :: isub1, isub2

! Variables
! ==--------------------------------------------------------------==
! ..All to all communication

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    IF (tr4a2a) THEN
       ! is that needed on P?
       CALL type_cast(xf, maxfft, xf4)
       CALL type_cast(yf, maxfft, yf4)       
       CALL mp_all2all( xf4, yf4, lda, comm )
    ELSE
       CALL mp_all2all( xf, yf, lda, comm )
    ENDIF
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fft_comm
  ! ==================================================================
  SUBROUTINE pack_x2y(xf,yf,nrays,lda,jrxpl,sp5,&
       maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Prepare data for sending
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL scopy_d(2*nrx,xf(nrs),1,yf4(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf4(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,xf(nrs),1,yf(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_x2y
  ! ==================================================================
  SUBROUTINE unpack_x2y(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    COMPLEX(real_8)                          :: xf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft
    COMPLEX(real_8)                          :: yf(maxfft)
    INTEGER                                  :: mproc, sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, mxrp

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Unpacking the data
    CALL zeroing(yf)!,maxfft)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_c(mxrp,xf4(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf4(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_no_omp(mxrp,xf(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_x2y
  ! ==================================================================
  SUBROUTINE phasen(f,kr1,kr2s,kr3s,n1u,n1o,nr2s,nr3s)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: kr1, kr2s, kr3s
    COMPLEX(real_8)                          :: f(kr1,kr2s,kr3s)
    INTEGER                                  :: n1u, n1o, nr2s, nr3s

    INTEGER                                  :: i, ii, ijk, isub, j, k
    REAL(real_8), DIMENSION(2)               :: pf = (/1._real_8,-1._real_8/)

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
    !$omp parallel do default(none) __COLLAPSE2 &
    !$omp             private(K,J,I,II,IJK) &
    !$omp             shared(F,PF,NR3S,NR2S,N1U,N1O)
    DO k=1,nr3s
       DO j=1,nr2s
          DO i=n1u,n1o
             ii=i-n1u+1
             ijk=MOD(k+j+i+1,2)+1
             f(ii,j,k)=f(ii,j,k)*pf(ijk)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phasen
  ! ==================================================================

  !TK
  ! Rewritten routines taking care of batches of states
  ! no support for A2A in single precision
  ! Author:
  ! Tobias Kloeffel, CCC,FAU Erlangen-Nuernberg tobias.kloeffel@fau.de
  ! Gerald Mathias, LRZ, Garching Gerald.Mathias@lrz.de
  ! Bernd Meyer, CCC, FAU Erlangen-Nuernberg bernd.meyer@fau.de

  SUBROUTINE putz_n(a,b,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    COMPLEX(real_8),INTENT(IN)               :: a(krmax-krmin+1,kr1,nperbatch,*)
    COMPLEX(real_8),INTENT(INOUT)            :: b(kr,kr1,kr2s,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'putz_n'
    REAL(real_8), POINTER __CONTIGUOUS       :: b_r(:,:,:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: a_r(:,:,:,:)
    INTEGER                                  :: isub1,isub2

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    
    CALL reshape_inplace(a,(/(krmax-krmin+1)*2,kr1,nperbatch,kr2s/),a_r)
    CALL reshape_inplace(b,(/kr*2,kr1,kr2s,nperbatch/),b_r)
    CALL putz_n_r(a_r,b_r,krmin,krmax,kr,kr1,kr2s,nperbatch)

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz_n

  SUBROUTINE getz_n(a,b,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    COMPLEX(real_8),INTENT(OUT)              :: b(krmax-krmin+1,kr1,nperbatch,*)
    COMPLEX(real_8),INTENT(IN)               :: a(kr,kr1,kr2s,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'getz_n'

    INTEGER                                  :: isub1,isub2,n,is,i,j,k,n1

! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    n=krmax-krmin+1
    n1=krmin-1
    !$omp parallel private (i,is,k,j) proc_bind(close)
    DO is=1,nperbatch
       !$omp do
       DO i=1,kr2s
          DO k=1,kr1
             !$omp simd
             DO j=1,n
                b(j,k,is,i)=a(n1+j,k,i,is)
             END DO
          END DO
       END DO
       !$omp end do nowait
    END DO
    !$omp end parallel

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE getz_n

  SUBROUTINE pack_y2x_n(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       tr4a2a,nperbatch,offset_in)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(OUT)              :: xf(*)
    COMPLEX(real_8),INTENT(IN)               :: yf(*)
    INTEGER,INTENT(IN)                       :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), mproc, &
                                                sp8(0:mproc-1), maxfft
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN),OPTIONAL              :: nperbatch,offset_in

    INTEGER                                  :: i, ii, ip, isub1, isub2, &
                                                jj, k, mxrp, is, nstate, offset

    CHARACTER(*),PARAMETER :: procedureN='PACK_Y2X_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    IF(PRESENT(nperbatch))THEN
       nstate=nperbatch
    ELSE
       nstate=1
    END IF
    IF(PRESENT(offset_in))THEN
       offset=offset_in
    ELSE
       offset=0
    END IF

    !$omp parallel private(is,ip,mxrp,i,ii,jj,k) proc_bind(close)
    DO is=1,nstate
       !$omp do
       DO ip=0,mproc-1
          mxrp = sp8(ip)
          DO i=1,lr1
             ii = ip*lda*nstate + (i-1)*mxrp + (is-1)*lda
             jj = (i-1)*m + (is-1)*offset
             DO k=1,mxrp
                xf(ii+k) = yf(jj+msp(k,ip+1))
             ENDDO
          ENDDO
       ENDDO
       !$omp end do nowait
    end do
    !$omp end parallel

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_y2x_n

  SUBROUTINE unpack_y2x_n(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc,tr4a2a,count)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN)               :: yf(*)
    COMPLEX(real_8),INTENT(OUT)              :: xf(*)
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN)                       :: m, nrays, lda, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1), maxfft
    INTEGER,INTENT(IN),OPTIONAL              :: count

    INTEGER                                  :: ip, ipp, isub1, isub2, k, nrs, &
                                                nrx, i,is,nstate

    CHARACTER(*),PARAMETER :: procedureN='UNPACK_Y2X_n'
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ..Pack the data for sending
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    If(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    END IF

    !$omp parallel do private(ip,i,is,nrs,ipp,k) proc_bind(close)
    DO ip=0,mproc-1
       DO i=1,sp5(ip)
          DO is=1,nstate
             nrs = (jrxpl(ip)-1)*nrays*nstate + (i-1)*nrays*nstate +(is-1)*nrays
             ipp = ip*lda*nstate + (i-1)*nrays + (is-1)*lda
             !$omp simd
             DO k=1,nrays
                xf(nrs+k)=yf(ipp+k)
             END DO
          END DO
       END DO
    END DO

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_y2x_n
  ! ==================================================================

  SUBROUTINE pack_x2y_n(xf,yf,nrays,lda,jrxpl,sp5,maxfft,mproc,tr4a2a,count)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN)               :: xf(*)
    COMPLEX(real_8),INTENT(OUT)              :: yf(*)
    INTEGER,INTENT(IN)                       :: nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN),OPTIONAL              :: count

    INTEGER                                  :: ip, ipp, isub1, isub2, &
                                                k, nrs, nrx, i,is,nstate
    CHARACTER(*),PARAMETER :: procedureN='PACK_X2Y_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    IF(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    END IF

    !$omp parallel do private(ip,i,is,nrs,ipp,k) proc_bind(close)
    DO ip=0,mproc-1
       DO i=1,sp5(ip)
          DO is=1,nstate
             nrs = (jrxpl(ip)-1)*nrays*nstate + (i-1)*nrays*nstate +(is-1)*nrays
             ipp = ip*lda*nstate + (i-1)*nrays + (is-1)*lda
             !$omp simd
             DO k=1,nrays
                yf(ipp+k)=xf(nrs+k)
             END DO
          END DO
       END DO
    END DO

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_x2y_n

  ! ==================================================================
  SUBROUTINE unpack_x2y_n(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a,&
       count,offset_in)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    COMPLEX(real_8),INTENT(IN)               :: xf(*)
    COMPLEX(real_8),INTENT(INOUT)            :: yf(*)
    INTEGER, INTENT(IN)                      :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft,&
                                                mproc, sp8(0:mproc-1)
    LOGICAL, INTENT(IN)                      :: tr4a2a
    INTEGER, INTENT(IN), OPTIONAL            :: count,offset_in
    INTEGER                                  :: isub1, isub2, nstate, offset
    REAL(real_8), POINTER __CONTIGUOUS       :: yf_r(:),xf_r(:)
    
    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_X2Y_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    ! ..Unpacking the data
    IF(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    ENDIF
    IF(PRESENT(offset_in))THEN
       offset=offset_in
    ELSE
       offset=0
    END IF
    CALL reshape_inplace(yf,(/nstate*offset*2/),yf_r)
    CALL reshape_inplace(xf,(/nstate*offset*2/),xf_r)
    CALL unpack_x2y_n_r(xf_r,yf_r,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       offset,nstate)
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_x2y_n
  ! ==================================================================
  SUBROUTINE unpack_x2y_n_r(xf_r,yf_r,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       offset,nstate)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    real(real_8),INTENT(IN)               :: xf_r(*)
    real(real_8),INTENT(OUT)              :: yf_r(*)
    INTEGER, INTENT(IN)                      :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft,&
                                                mproc, sp8(0:mproc-1)
    INTEGER, INTENT(IN)                    :: offset,nstate
    INTEGER                                  :: i, ii, ip, isub1, jj, k, mxrp,is,kk

    !$omp parallel private(is,ip,mxrp,i,ii,jj,k,kk) proc_bind(close)
    !    call zero(yf_r,nstate*offset*2)
    !$omp do
    do is=1,nstate*offset*2
       yf_r(is)=0._real_8
    end do
    DO is=1,nstate
       !$omp do schedule (static)
       DO ip=0,mproc-1
          mxrp = sp8(ip)
          DO i=1,lr1
             ii = ip*lda*nstate + (i-1)*mxrp + (is-1)*lda
             jj = (i-1)*m + (is-1)*offset
             DO k=1,mxrp
                do kk=-1,0
                   yf_r((jj+msp(k,ip+1))*2+kk) = xf_r((ii+k)*2+kk)
                end do
             END DO
          END DO
       END DO
       !$omp end do nowait
    END DO
    !$omp end parallel

  END SUBROUTINE unpack_x2y_n_r
  SUBROUTINE putz_n_r(a_r,b_r,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    REAL(real_8),INTENT(IN)               :: a_r((krmax-krmin+1)*2,kr1,nperbatch,*)
    REAL(real_8),INTENT(OUT)              :: b_r(kr*2,kr1,kr2s,*)

    INTEGER                                  :: isub,n,n1,n2,n3,n4,is,i,j,k,krmin_loc,kr_loc

    
    n=krmax-krmin+1
    n1=(krmin-1)*2
    n2=(krmin-1+n)*2
    n3=(krmin+n)*2-1
    n4=(krmin-1)*2
    kr_loc=kr*2
    krmin_loc=krmin*2-1

    !$omp parallel private (i,is,k,j) proc_bind(close)
    DO is=1,nperbatch
       !$omp do
       DO i=1,kr2s
          DO k=1,kr1
             !$omp simd
             DO j=1,n1
                b_r(j,k,i,is)=0.0_real_8
             END DO
             !$omp simd
             DO j=krmin_loc,n2
                b_r(j,k,i,is)=a_r(j-n4,k,is,i)
             END DO
             !$omp simd
             DO j=n3,kr_loc
                b_r(j,k,i,is)=0.0_real_8
             END DO
          END DO
       END DO
       !$omp end do nowait
    END DO
    !$omp end parallel

    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz_n_r

  subroutine zero(in,len)
    real(real_8), intent(out) __CONTIGUOUS :: in(:)
    integer, intent(in) :: len
    integer :: i
    !$omp do simd
    do i=1,len
       in(i)=0.0_real_8
    end do
  end subroutine zero
  subroutine zero_noomp(in,len)
    real(real_8), intent(out) __CONTIGUOUS :: in(:)
    integer, intent(in) :: len
    integer :: i
    !$omp simd
    do i=1,len
       in(i)=0.0_real_8
    end do
  end subroutine zero_noomp
  !TK
  !CLR
  SUBROUTINE set_psi_new_gdist( tfft, psi, aux, remswitch, mythread, last_single, counter )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) ::tfft
    INTEGER, INTENT(IN) :: remswitch, mythread, counter
    COMPLEX(real_8), INTENT(IN)  :: psi ( : , : )
    COMPLEX(real_8), INTENT(OUT)  :: aux ( fpar%kr3s , * )
    LOGICAL, INTENT(IN) :: last_single

    INTEGER :: j, i, iter, l, lter, f, fter
    INTEGER :: offset, offset2, offset3, offset4, offset5, offset6
    CHARACTER(*), PARAMETER :: procedureN = 'set_psi_new_gdist'

    IF( .not. last_single ) THEN

       DO i = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
          iter = mod( i-1, tfft%nsw(parai%me+1) ) + 1
          offset  = ( iter - 1 ) * fpar%kr3s
          offset2 = 2 * ( ( (i-1) / tfft%nsw(parai%me+1) ) + 1 )

          DO j = tfft%map_set_psi(1,iter), tfft%map_set_psi(2,iter)
             aux( j, i ) = (0.d0, 0.d0)
          ENDDO
          DO j = tfft%map_set_psi(3,iter), tfft%map_set_psi(4,iter)
             aux( j, i ) = conjg( psi( indz_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset + j ), offset2 ) )
          ENDDO
          DO j = tfft%map_set_psi(5,iter), tfft%map_set_psi(6,iter)
             aux( j, i ) = psi( nzh_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset + j ), offset2 )
          ENDDO

          DO j = tfft%map_set_psi(7,iter), tfft%map_set_psi(8,iter)
             aux( j, i ) = (0.d0, 0.d0)
          ENDDO

          DO j = tfft%map_set_psi(9,iter), tfft%map_set_psi(10,iter)
             aux( j, i ) = psi( nzh_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset + j ), offset2 )
          ENDDO
          DO j = tfft%map_set_psi(11,iter), tfft%map_set_psi(12,iter)
             aux( j, i ) = conjg( psi( indz_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset + j ), offset2 ) )
          ENDDO
          DO j = tfft%map_set_psi(13,iter), tfft%map_set_psi(14,iter)
             aux( j, i ) = (0.d0, 0.d0)
          ENDDO

       ENDDO

    ELSE

       DO l = tfft%thread_prepare_start( mythread+1, 1 ), tfft%thread_prepare_end( mythread+1, 1 )
          lter = mod( l-1, tfft%nsw(parai%me+1) ) + 1
          offset3  = ( lter - 1 ) * fpar%kr3s
          offset4 = 2 * ( ( (l-1) / tfft%nsw(parai%me+1) ) + 1 )

          DO j = tfft%map_set_psi(1,lter), tfft%map_set_psi(2,lter)
             aux( j, i ) = (0.d0, 0.d0)
          ENDDO
          DO j = tfft%map_set_psi(3,lter), tfft%map_set_psi(4,lter)
             aux( j, l ) = conjg( psi( indz_r( offset3 + j ), offset4 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset3 + j ), offset4 ) )
          ENDDO
          DO j = tfft%map_set_psi(5,lter), tfft%map_set_psi(6,lter)
             aux( j, l ) = psi( nzh_r( offset3 + j ), offset4 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset3 + j ), offset4 )
          ENDDO

          DO j = tfft%map_set_psi(7,lter), tfft%map_set_psi(8,lter)
             aux( j, l ) = (0.d0, 0.d0)
          ENDDO

          DO j = tfft%map_set_psi(9,lter), tfft%map_set_psi(10,lter)
             aux( j, l ) = psi( nzh_r( offset3 + j ), offset4 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset3 + j ), offset4 )
          ENDDO
          DO j = tfft%map_set_psi(11,lter), tfft%map_set_psi(12,lter)
             aux( j, l ) = conjg( psi( indz_r( offset3 + j ), offset4 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset3 + j ), offset4 ) )
          ENDDO
          DO j = tfft%map_set_psi(13,lter), tfft%map_set_psi(14,lter)
             aux( j, i ) = (0.d0, 0.d0)
          ENDDO

       ENDDO

       DO f = tfft%thread_prepare_start( mythread+1, 2 ), tfft%thread_prepare_end( mythread+1, 2 )
          fter = mod( f-1, tfft%nsw(parai%me+1) ) + 1
          offset5  = ( fter - 1 ) * fpar%kr3s
          offset6 = 2 * ( ( (f-1) / tfft%nsw(parai%me+1) ) + 1 ) - 1

          DO j = tfft%map_set_psi(1,fter), tfft%map_set_psi(2,fter)
             aux( j, i ) = (0.d0, 0.d0)
          ENDDO
          DO j = tfft%map_set_psi(3,fter), tfft%map_set_psi(4,fter)
             aux( j, f ) = conjg( psi( indz_r( offset5 + j ), offset6 ) )
          ENDDO
          DO j = tfft%map_set_psi(5,fter), tfft%map_set_psi(6,fter)
             aux( j, f ) = psi( nzh_r( offset5 + j ), offset6 )
          ENDDO

          DO j = tfft%map_set_psi(7,fter), tfft%map_set_psi(8,fter)
             aux( j, f ) = (0.d0, 0.d0)
          ENDDO

          DO j = tfft%map_set_psi(9,fter), tfft%map_set_psi(10,fter)
             aux( j, f ) = psi( nzh_r( offset5 + j ), offset6 )
          ENDDO
          DO j = tfft%map_set_psi(11,fter), tfft%map_set_psi(12,fter)
             aux( j, f ) = conjg( psi( indz_r( offset5 + j ), offset6 ) )
          ENDDO
          DO j = tfft%map_set_psi(13,fter), tfft%map_set_psi(14,fter)
             aux( j, i ) = (0.d0, 0.d0)
          ENDDO

       ENDDO

       !$  locks_omp( mythread+1, counter, 4 ) = .false.
       !$omp flush( locks_omp )
       !$  DO WHILE( ANY( locks_omp( :, counter, 4 ) ) )
       !$omp flush( locks_omp )
       !$  END DO

    END IF

  END SUBROUTINE set_psi_new_gdist

  SUBROUTINE fft_comm_preinitialized( tfft, remswitch, work_buffer, which )
    IMPLICIT NONE

    INTEGER, INTENT(IN)                         :: remswitch, work_buffer, which
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT)    :: tfft

    CHARACTER(*), PARAMETER :: procedureN = 'fft_com_preinit'

    INTEGER :: ierr, isub, isub4

!    IF( cntl%fft_tune_batchsize ) THEN
!       CALL tiset(procedureN//'_tuning',isub4)
!    ELSE
!       CALL tiset(procedureN,isub)
!    END IF

    CALL MP_STARTALL( tfft%comm_sendrecv(1,tfft%which)+tfft%comm_sendrecv(2,tfft%which), parai%sendrecv_handle(:,work_buffer,remswitch,which) )

    CALL MP_WAITALL( tfft%comm_sendrecv(1,tfft%which)+tfft%comm_sendrecv(2,tfft%which), parai%sendrecv_handle(:,work_buffer,remswitch,which) )

!    IF( cntl%fft_tune_batchsize ) THEN
!       CALL tihalt(procedureN//'_tuning',isub4)
!    ELSE
!       CALL tihalt(procedureN,isub)
!    END IF

  END SUBROUTINE fft_comm_preinitialized

  SUBROUTINE fft_comm_ALL2ALL( tfft, remswitch, work_buffer, which, comm_send, comm_recv, sendsize )
    USE mpi_f08
    IMPLICIT NONE

    INTEGER, INTENT(IN)                         :: remswitch, work_buffer, which, sendsize
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT)    :: tfft
!    COMPLEX(real_8), INTENT(INOUT)              :: comm_send(*), comm_recv(*)
    COMPLEX(real_8), INTENT(INOUT)              :: comm_send(:), comm_recv(:)

    CHARACTER(*), PARAMETER :: procedureN = 'fft_comm_ALL2ALL'

    INTEGER :: ierr, isub, isub4

!    IF( cntl%fft_tune_batchsize ) THEN
!       CALL tiset(procedureN//'_tuning',isub4)
!    ELSE
!       CALL tiset(procedureN,isub)
!    END IF

    CALL MPI_ALLTOALL( comm_send, sendsize, MPI_DOUBLE_COMPLEX, comm_recv, sendsize, MPI_DOUBLE_COMPLEX, parai%allgrp )

!    IF( cntl%fft_tune_batchsize ) THEN
!       CALL tihalt(procedureN//'_tuning',isub4)
!    ELSE
!       CALL tihalt(procedureN,isub)
!    END IF

  END SUBROUTINE fft_comm_ALL2ALL 

  SUBROUTINE invfft_z_section( tfft, aux, comm_mem_send, comm_mem_recv, batch_size, remswitch, mythread, nss, current )
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: batch_size, remswitch, mythread, current
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    COMPLEX(real_8), INTENT(INOUT) :: comm_mem_send( * ), comm_mem_recv( * )
    COMPLEX(real_8), INTENT(INOUT)  :: aux ( fpar%kr3s , * )
    INTEGER, INTENT(IN) :: nss(*)

    INTEGER :: l, m, j, k, i
    INTEGER :: offset, kdest, ierr
    CHARACTER(*), PARAMETER :: procedureN = 'invfft_z_section'

  !------------------------------------------------------
  !------------z-FFT Start-------------------------------

    CALL mltfft_fftw_threadsafe('n','n',aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     fpar%kr3s, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     fpar%kr3s, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     fpar%kr3s, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which),-1,scal,.FALSE.,mythread,parai%ncpus_FFT)

  !-------------z-FFT End--------------------------------
  !------------------------------------------------------

  !------------------------------------------------------
  !-----------pack_z2y Start-----------------------------

    CALL z2y_fillsend( tfft, aux, comm_mem_send, batch_size, remswitch, mythread, nss )

    IF( tfft%which .eq. 1 .and. .not. cntl%fft_distmem ) THEN

       !In theory, locks could be made faster by checking each iset individually for non-remainder cases
       !$omp flush( locks_calc_1 )
       !$  DO WHILE( ANY(locks_calc_1( :, 1+current:fft_batchsize+current ) ) )
       !$omp flush( locks_calc_1 )
       !$  END DO

    END IF

    CALL z2y_fillrecv( tfft, aux, comm_mem_recv, batch_size, remswitch, mythread, nss )

  !------------pack_z2y End------------------------------
  !------------------------------------------------------

  END SUBROUTINE invfft_z_section

  SUBROUTINE invfft_y_section( tfft, comm_mem_recv, aux, map_z2y, mythread, my_nr1s, ispec, counter )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) ::tfft
    INTEGER, INTENT(IN) :: mythread, my_nr1s, ispec, counter
    COMPLEX(real_8), INTENT(IN)  :: comm_mem_recv( * )
    COMPLEX(real_8), INTENT(INOUT) :: aux( fpar%kr2s , * )
    INTEGER, INTENT(IN) :: map_z2y( * )

    INTEGER :: i, k, offset, iter
    CHARACTER(*), PARAMETER :: procedureN = 'invfft_y_section'

  !------------------------------------------------------
  !----------unpack_z2y Start----------------------------

    offset = (ispec-1)*tfft%big_chunks( tfft%which )

    DO i = tfft%thread_y_start( mythread+1, tfft%which ), tfft%thread_y_end( mythread+1, tfft%which )
       iter = mod( i-1, my_nr1s ) + 1
       DO k = tfft%map_z2y_bounds( iter, 1, tfft%which ), tfft%map_z2y_bounds( iter, 2, tfft%which ) - 1
          aux( k, i ) = (0.0_real_8,0.0_real_8)
       END DO
       DO k = tfft%map_z2y_bounds( iter, 5, tfft%which ), tfft%map_z2y_bounds( iter, 6, tfft%which ) - 1
          aux( k, i ) = comm_mem_recv( map_z2y( (i-1) * fpar%kr2s + k ) + offset )
       END DO
       DO k = tfft%map_z2y_bounds( iter, 3, tfft%which ), tfft%map_z2y_bounds( iter, 4, tfft%which ) - 1
          aux( k, i ) = (0.0_real_8,0.0_real_8)
       END DO
       DO k = tfft%map_z2y_bounds( iter, 7, tfft%which ), tfft%map_z2y_bounds( iter, 8, tfft%which ) - 1
          aux( k, i ) = comm_mem_recv( map_z2y( (i-1) * fpar%kr2s + k ) + offset )
       END DO
    END DO

  !-----------unpack_z2y End-----------------------------
  !------------------------------------------------------

    IF( tfft%which .eq. 1 ) THEN
       !$  locks_omp_big( mythread+1, ispec, counter, 5 ) = .false.
       !$omp flush( locks_omp_big )
    END IF

  !------------------------------------------------------
  !------------y-FFT Start-------------------------------

    CALL mltfft_fftw_threadsafe('n','n',aux( : , tfft%thread_y_start( mythread+1, tfft%which ) : tfft%thread_y_end( mythread+1, tfft%which ) ), &
                     fpar%kr2s, tfft%thread_y_sticks(mythread+1,tfft%which), &
                     aux( : , tfft%thread_y_start( mythread+1, tfft%which ) : tfft%thread_y_end( mythread+1, tfft%which ) ), &
                     fpar%kr2s, tfft%thread_y_sticks(mythread+1,tfft%which), &
                     fpar%kr2s, tfft%thread_y_sticks(mythread+1,tfft%which),-1,scal,.FALSE.,mythread,parai%ncpus_FFT)

  !-------------y-FFT End--------------------------------
  !------------------------------------------------------

  END SUBROUTINE invfft_y_section

  SUBROUTINE invfft_x_section( tfft, aux2, aux_r, mythread, my_nr1s )
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: mythread, my_nr1s
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    COMPLEX(real_8), INTENT(INOUT) :: aux2( * )
    COMPLEX(real_8), INTENT(INOUT) :: aux_r( : )

    CHARACTER(*), PARAMETER :: procedureN = 'invfft_x_section'

    INTEGER :: i, k, offset

    Call First_Part_x_section( aux_r )

    Call Second_Part_x_section( aux_r )

    CONTAINS

      SUBROUTINE First_Part_x_section( aux_r )

        IMPLICIT NONE
        COMPLEX(real_8), INTENT(INOUT) :: aux_r( * )

        !------------------------------------------------------
        !---------transpose y2x Start--------------------------

          DO i = tfft%thread_x_start( mythread+1, tfft%which ), tfft%thread_x_end( mythread+1, tfft%which )
             offset = (i-1) * fpar%kr1s
             DO k = 1, tfft%zero_transpose_y2x_start( tfft%which ) - 1
                aux_r( offset + k ) = aux2( tfft%map_transpose_y2x( offset + k, tfft%which ) )
             END DO
             DO k = tfft%zero_transpose_y2x_start( tfft%which ), tfft%zero_transpose_y2x_end( tfft%which )
                aux_r( offset + k ) = (0.0_real_8, 0.0_real_8)
             END DO
             DO k = tfft%zero_transpose_y2x_end( tfft%which ) + 1, fpar%kr1s
                aux_r( offset + k ) = aux2( tfft%map_transpose_y2x( offset + k, tfft%which ) )
             END DO
          END DO

        !----------transpose y2x End---------------------------
        !------------------------------------------------------

      END SUBROUTINE First_Part_x_section

      SUBROUTINE Second_Part_x_section( aux_r )

        Implicit NONE
        COMPLEX(real_8), INTENT(INOUT) :: aux_r( fpar%kr1s , * )

        !------------------------------------------------------
        !------------x-FFT Start-------------------------------

          CALL mltfft_fftw_threadsafe('n','n',aux_r( : , tfft%thread_x_start( mythread+1, tfft%which ) : tfft%thread_x_end( mythread+1, tfft%which ) ), &
                           fpar%kr1s, tfft%thread_x_sticks(mythread+1, tfft%which), &
                           aux_r( :, tfft%thread_x_start( mythread+1, tfft%which ) : tfft%thread_x_end( mythread+1, tfft%which ) ), &
                           fpar%kr1s, tfft%thread_x_sticks(mythread+1, tfft%which), &
                           fpar%kr1s, tfft%thread_x_sticks(mythread+1, tfft%which),-1,scal,.FALSE.,mythread,parai%ncpus_FFT)

        !-------------x-FFT End--------------------------------
        !------------------------------------------------------

      END SUBROUTINE Second_Part_x_section

  END SUBROUTINE invfft_x_section

  SUBROUTINE fwfft_x_section( tfft, aux, mythread )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: mythread
    COMPLEX(real_8), INTENT(INOUT) :: aux( fpar%kr1s , * )

    CHARACTER(*), PARAMETER :: procedureN = 'fwfft_x_section'

  !------------------------------------------------------
  !------------x-FFT Start-------------------------------

    CALL mltfft_fftw_threadsafe('n','n',aux( : , tfft%thread_x_start( mythread+1, tfft%which ) : tfft%thread_x_end( mythread+1, tfft%which ) ), &
                     fpar%kr1s, tfft%thread_x_sticks(mythread+1, tfft%which), &
                     aux( : , tfft%thread_x_start( mythread+1,  tfft%which ) : tfft%thread_x_end( mythread+1, tfft%which ) ), &
                     fpar%kr1s, tfft%thread_x_sticks(mythread+1, tfft%which), &
                     fpar%kr1s, tfft%thread_x_sticks(mythread+1, tfft%which),1,scal,.FALSE.,mythread,parai%ncpus_FFT)

  !-------------x-FFT End--------------------------------
  !------------------------------------------------------

  END SUBROUTINE fwfft_x_section

  SUBROUTINE fwfft_y_section( tfft, aux, aux2_r, comm_mem_send, comm_mem_recv, map_y2z, batch_size, ispec, counter, mythread )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: counter, batch_size, mythread, ispec
    COMPLEX(real_8), INTENT(INOUT)  :: comm_mem_send( * )
    COMPLEX(real_8), INTENT(INOUT)  :: comm_mem_recv( * )
    COMPLEX(real_8), INTENT(INOUT) :: aux( * )
    COMPLEX(real_8), INTENT(INOUT) :: aux2_r( : )
    INTEGER, INTENT(IN) :: map_y2z( * )

    INTEGER :: l, m, i, offset, j, k, ibatch, jter, offset2
    CHARACTER(*), PARAMETER :: procedureN = 'fwfft_y_section'

    Call First_Part_y_section( aux2_r )

    Call Second_Part_y_section( aux2_r )

    !$  locks_omp_big( mythread+1, ispec, counter, 6 ) = .false.
    !$omp flush( locks_omp_big )
    !$  DO WHILE( ANY( locks_omp_big( :, ispec, counter, 6 ) ) )
    !$omp flush( locks_omp_big )
    !$  END DO

    Call Third_Part_y_section( aux2_r )

    CONTAINS

      SUBROUTINE First_Part_y_section( aux2 )
        Implicit NONE
        COMPLEX(real_8), INTENT(INOUT) :: aux2( * )

        !------------------------------------------------------
        !---------transpose x2y Start--------------------------

          DO i = tfft%thread_y_start( mythread+1, tfft%which ), tfft%thread_y_end( mythread+1, tfft%which )
             DO j = 1, fpar%kr2s
                aux2( j + (i-1) * fpar%kr2s ) = aux( tfft%map_transpose_x2y( j + (i-1) * fpar%kr2s, tfft%which ) )
             ENDDO
          ENDDO

        !----------transpose x2y End---------------------------
        !------------------------------------------------------

      END SUBROUTINE First_Part_y_section

      SUBROUTINE Second_Part_y_section( aux2 )
        Implicit NONE
        COMPLEX(real_8), INTENT(INOUT) :: aux2( fpar%kr2s, * )

        !------------------------------------------------------
        !------------y-FFT Start-------------------------------

          CALL mltfft_fftw_threadsafe('n','n',aux2( : , tfft%thread_y_start( mythread+1, tfft%which ) : tfft%thread_y_end( mythread+1, tfft%which ) ), &
                           fpar%kr2s, tfft%thread_y_sticks(mythread+1,tfft%which), &
                           aux2( : , tfft%thread_y_start( mythread+1, tfft%which ) : tfft%thread_y_end( mythread+1, tfft%which ) ), &
                           fpar%kr2s, tfft%thread_y_sticks(mythread+1,tfft%which), &
                           fpar%kr2s, tfft%thread_y_sticks(mythread+1,tfft%which),1,scal,.FALSE.,mythread,parai%ncpus_FFT)

        !-------------y-FFT End--------------------------------
        !------------------------------------------------------

      END SUBROUTINE Second_Part_y_section

      SUBROUTINE Third_Part_y_section( aux2 )

        IMPLICIT NONE
        COMPLEX(real_8), INTENT(INOUT) :: aux2  ( * )

        !------------------------------------------------------
        !-----------pack_y2z Start-----------------------------

          CALL y2z_fillsend( tfft, aux2, comm_mem_send, batch_size, mythread, ispec, map_y2z )

          IF( tfft%which .eq. 1 .and. .not. cntl%fft_distmem ) THEN

             !In theory, locks could be made faster by checking each iset individually for non-remainder cases
             !$omp flush( locks_calc_2 )
             !$  DO WHILE( ANY(locks_calc_2(:,1+(counter-1)*fft_batchsize:ispec+(counter-1)*fft_batchsize ) ) )
             !$omp flush( locks_calc_2 )
             !$  END DO

          END IF

          CALL y2z_fillrecv( tfft, aux2, comm_mem_recv, batch_size, mythread, ispec, map_y2z )

        !------------pack_y2z End------------------------------
        !------------------------------------------------------

      END SUBROUTINE Third_Part_y_section

  END SUBROUTINE fwfft_y_section

  SUBROUTINE fwfft_z_section( tfft, comm_mem_recv, aux, counter, batch_size, remswitch, mythread, nss, factor_in )
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: counter, batch_size, remswitch, mythread
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    COMPLEX(real_8), INTENT(IN)  :: comm_mem_recv( * )
    COMPLEX(real_8), INTENT(INOUT)  :: aux ( fpar%kr3s , * )
    INTEGER, INTENT(IN) :: nss( * )
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: factor_in

    DOUBLE PRECISION :: factor
    INTEGER :: j, l, k, i, m
    INTEGER :: offset, kfrom, ierr
    CHARACTER(*), PARAMETER :: procedureN = 'fwfft_z_section'

    IF( present( factor_in ) ) THEN
       factor = factor_in
    ELSE
       factor = 1
    END IF

  !------------------------------------------------------
  !----------unpack_y2z Start----------------------------

    CALL y2z_voidrecv( tfft, comm_mem_recv, aux, batch_size, remswitch, mythread, factor, nss )

  !-----------unpack_y2z End-----------------------------
  !------------------------------------------------------

    IF( tfft%which .eq. 1 .and. .not. cntl%fft_distmem ) THEN
       !$  locks_omp( mythread+1, counter, 5 ) = .false.
       !$omp flush( locks_omp )
       IF( parai%ncpus_FFT .eq. 1 .or. .not. ANY( locks_omp( :, counter, 5 ) ) ) THEN
          IF( cntl%krwfn ) THEN
          !$   locks_calc_2( parai%node_me+1, 1+(counter+fft_numbuff-1)*fft_batchsize:batch_size+(counter+fft_numbuff-1)*fft_batchsize ) = .false.
          !$omp flush( locks_calc_2 )
          ELSE
          !$   locks_calc_1( parai%node_me+1, 1+(counter+fft_numbuff-1)*fft_batchsize:fft_batchsize+(counter+fft_numbuff-1)*fft_batchsize ) = .false.
          !$omp flush( locks_calc_1 )
          END IF
       END IF
    END IF

  !------------------------------------------------------
  !------------z-FFT Start-------------------------------

    CALL mltfft_fftw_threadsafe('n','n',aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     fpar%kr3s, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     fpar%kr3s, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     fpar%kr3s, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which),1,scal,.FALSE.,mythread,parai%ncpus_FFT)

  !-------------z-FFT End--------------------------------
  !------------------------------------------------------

  END SUBROUTINE fwfft_z_section
!CLR
END MODULE fftutil_utils
