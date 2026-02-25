MODULE fftnew_utils
  USE cell,                            ONLY: cell_com
  USE cnst,                            ONLY: pi
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             indzs,&
                                             inyh,&
                                             nzh,&
                                             nzhs
  USE cuda_types,                      ONLY: cuda_memory_t
  USE cuda_utils,                      ONLY: cuda_memcpy_host_to_device
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: &
       fftpool, fftpoolsize, fpoolv, inzf, inzfp, inzh, inzhp, inzs, inzsp, &
       jgw, jgws, jhg, jhgs, kr2max, kr2min, kr3max, kr3min, lfrm, llr1, &
       lmsq, lmsqmax, lnzf, lnzs, lr1, lr1m, lr1s, lr2, lr2s, lr3, lr3s, &
       lrxpl, lrxpool, lsrm, mfrays, mg, msp, msqf, msqfpool, msqs, msqspool, &
       msrays, mz, ngrm, nhrm, nr1m, nzff, nzffp, nzfs, nzfsp, qr1, qr1s, &
       qr2, qr2max, qr2min, qr2s, qr3, qr3max, qr3min, qr3s, sp5, sp8, sp9, &
       spm, FFT_TYPE_DESCRIPTOR, fft_batchsize, fft_residual, fft_numbatches, &
       fft_numbuff
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE fftchk_utils,                    ONLY: fftchk
  USE kinds,                           ONLY: real_8
  USE loadpa_utils,                    ONLY: leadim
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_win_alloc_shared_mem,&
                                             mp_send_init_COMPLEX,&
                                             mp_recv_init_COMPLEX
  USE parac,                           ONLY: parai,&
                                             paral
  USE store_types,                     ONLY: restart1 
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar,&
                                             cntl,&
                                             cnti
  USE utils,                           ONLY: icopy
  USE zeroing_utils,                   ONLY: zeroing

  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setfftn
  !public :: rmfftnset
  PUBLIC :: addfftnset
  !public :: setrays
  PUBLIC :: fft_new_gdist_batch_setup
  PUBLIC :: fft_new_gdist_setup
  PUBLIC :: Make_Manual_Maps
  PUBLIC :: Make_z2y_Maps
  PUBLIC :: Pre_Initialize_C2_Com

  COMPLEX(real_8), POINTER, SAVE, CONTIGUOUS :: comm_send(:,:)
  PUBLIC :: comm_send
  COMPLEX(real_8), POINTER, SAVE, CONTIGUOUS :: comm_recv(:,:)
  PUBLIC :: comm_recv
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_inv(:,:)
  PUBLIC :: locks_calc_inv
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_com_inv(:,:)
  PUBLIC :: locks_com_inv
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_fw(:,:)
  PUBLIC :: locks_calc_fw
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_com_fw(:,:)
  PUBLIC :: locks_com_fw
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_cc_invfw(:,:,:)
  PUBLIC :: locks_cc_invfw
  LOGICAL, ALLOCATABLE, SAVE :: locks_omp(:,:,:)
  PUBLIC :: locks_omp
  LOGICAL, ALLOCATABLE, SAVE :: locks_omp_big(:,:,:,:)
  PUBLIC :: locks_omp_big
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_1(:,:)
  PUBLIC :: locks_calc_1
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_2(:,:)
  PUBLIC :: locks_calc_2


CONTAINS

  ! ==================================================================
  SUBROUTINE setfftn(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    INTEGER                                  :: i_device, ip, ipx
    TYPE(cuda_memory_t), POINTER             :: inzs_d, lrxpl_d, msqf_d, &
                                                msqs_d, nzfs_d, sp5_d, sp8_d, &
                                                sp9_d

    IF (ipool.EQ.0) THEN
       msrays   = parai%ngrays
       mfrays   = parai%nhrays
       llr1     = fpar%nnr1
       qr1s     = fpar%kr1s
       qr2s     = fpar%kr2s
       qr3s     = fpar%kr3s
       qr1      = fpar%kr1
       lr1s     = spar%nr1s
       lr2s     = spar%nr2s
       lr3s     = spar%nr3s
       lr1      = parm%nr1
       qr2max   = kr2max
       qr2min   = kr2min
       qr3max   = kr3max
       qr3min   = kr3min
       lsrm     = ngrm
       lfrm     = nhrm
       lr1m     = nr1m
       lmsq     = nhrm
       maxfftn  = maxfft
       jgw      = ncpw%ngw
       jgws     = spar%ngws
       jhg      = ncpw%nhg
       jhgs     = spar%nhgs
       nzff => nzh
       inzf => indz
       nzfs => nzhs
       inzs => indzs
       inzh => inyh

       DO ip=0,parai%nproc-1
          lrxpl(ip,1)    = parap%nrxpl(ip,1)
          lrxpl(ip,2)    = parap%nrxpl(ip,2)
          sp5(ip)        = parap%sparm(5,ip)
          sp8(ip)        = parap%sparm(8,ip)
          sp9(ip)        = parap%sparm(9,ip)
          ipx=lmsq*ip
          IF (nhrm.GT.0) THEN
             CALL icopy(nhrm,msp(1,1,ip+1),1,msqf(ipx+1),1)
             CALL icopy(nhrm,msp(1,2,ip+1),1,msqs(ipx+1),1)
          ENDIF
       ENDDO
    ELSEIF (ipool.GT.0 .AND. ipool.LE.fftpool) THEN
       ! LOAD FROM POOL
       msrays   = fpoolv( 1,ipool)
       mfrays   = fpoolv( 2,ipool)
       llr1     = fpoolv( 3,ipool)
       qr1s     = fpoolv( 4,ipool)
       qr2s     = fpoolv( 5,ipool)
       qr3s     = fpoolv( 6,ipool)
       qr1      = fpoolv( 7,ipool)
       qr2      = fpoolv( 8,ipool)
       qr3      = fpoolv( 9,ipool)
       lr1s     = fpoolv(10,ipool)
       lr2s     = fpoolv(11,ipool)
       lr3s     = fpoolv(12,ipool)
       lr1      = fpoolv(13,ipool)
       lr2      = fpoolv(14,ipool)
       lr3      = fpoolv(15,ipool)
       qr2max   = fpoolv(16,ipool)
       qr2min   = fpoolv(17,ipool)
       qr3max   = fpoolv(18,ipool)
       qr3min   = fpoolv(19,ipool)
       lsrm     = fpoolv(20,ipool)
       lfrm     = fpoolv(21,ipool)
       lr1m     = fpoolv(22,ipool)
       lmsq     = fpoolv(23,ipool)
       maxfftn  = fpoolv(24,ipool)
       jgw      = fpoolv(25,ipool)
       jgws     = fpoolv(26,ipool)
       jhg      = fpoolv(27,ipool)
       jhgs     = fpoolv(28,ipool)
       nzff => nzffp(:, ipool)
       inzf => inzfp(:, ipool)
       nzfs => nzfsp(:, ipool)
       inzs => inzsp(:, ipool)
       inzh => inzhp(:, :, ipool)

       DO ip=0,parai%nproc-1
          lrxpl(ip,1)    = lrxpool(ip,1,ipool)
          lrxpl(ip,2)    = lrxpool(ip,2,ipool)
          sp5(ip)        = spm(5,ip,ipool)
          sp8(ip)        = spm(8,ip,ipool)
          sp9(ip)        = spm(9,ip,ipool)
          ipx=lmsq*ip
          IF (lmsqmax.GT.0) THEN
             CALL icopy(lmsq,msqfpool(1,ip+1,ipool),1,msqf(ipx+1),1)
             CALL icopy(lmsq,msqspool(1,ip+1,ipool),1,msqs(ipx+1),1)
          ENDIF
       ENDDO
    ELSE
       CALL stopgm("SETFFTN","FFTPOOL NOT DEFINED",& 
            __LINE__,__FILE__)
    ENDIF

    !vw copy FFT arrays to GPU memory
    IF( cp_cuda_env%use_fft ) THEN
       DO i_device = 1, cp_cuda_env%fft_n_devices_per_task
          CALL cp_cufft_device_get_ptrs ( cp_cufft, i_device, sp5_d=sp5_d, sp8_d=sp8_d, sp9_d=sp9_d, &
               & msqs_d=msqs_d, msqf_d=msqf_d, lrxpl_d=lrxpl_d, nzfs_d=nzfs_d, inzs_d=inzs_d )
          CALL cuda_memcpy_host_to_device ( sp5, sp5_d )
          CALL cuda_memcpy_host_to_device ( sp8, sp8_d )
          CALL cuda_memcpy_host_to_device ( sp9, sp9_d )
          CALL cuda_memcpy_host_to_device ( lrxpl, lrxpl_d )
          CALL cuda_memcpy_host_to_device ( msqf, msqf_d )
          CALL cuda_memcpy_host_to_device ( msqs, msqs_d )
          CALL cuda_memcpy_host_to_device ( nzfs, nzfs_d )
          CALL cuda_memcpy_host_to_device ( inzs, inzs_d )
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE setfftn
  ! ==================================================================
  SUBROUTINE rmfftnset(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    INTEGER                                  :: i, ip, j

    IF (ipool.GT.0 .AND. ipool.LE.fftpool) THEN
       DO ip=ipool+1,fftpool
          i=ip
          j=ip-1
          CALL icopy(28,fpoolv(1,i),1,fpoolv(1,j),1)
          CALL icopy(lnzf,nzffp(1,i),1,nzffp(1,j),1)
          CALL icopy(lnzf,inzfp(1,i),1,inzfp(1,j),1)
          CALL icopy(lnzs,nzfsp(1,i),1,nzfsp(1,j),1)
          CALL icopy(lnzs,inzsp(1,i),1,inzsp(1,j),1)
          CALL icopy(3*lnzf,inzhp(1,1,i),1,inzhp(1,1,j),1)
          CALL icopy(2*SIZE(lrxpool,1),lrxpool(0,1,i),1,lrxpool(0,1,j),1)
          CALL icopy(9*SIZE(lrxpool,1),spm(1,0,i),1,spm(1,0,j),1)
          IF (lmsqmax.GT.0) THEN
             CALL icopy(lmsqmax*parai%nproc,msqfpool(1,1,i),1,&
                  msqfpool(1,1,j),1)
             CALL icopy(lmsqmax*parai%nproc,msqspool(1,1,i),1,&
                  msqspool(1,1,j),1)
          ENDIF
       ENDDO
       CALL zeroing(fpoolv(:,fftpool))!,28)
       CALL zeroing(nzffp(:,fftpool))!,lnzf)
       CALL zeroing(inzfp(:,fftpool))!,lnzf)
       CALL zeroing(nzfsp(:,fftpool))!,lnzf)
       CALL zeroing(inzsp(:,fftpool))!,lnzf)
       CALL zeroing(inzhp(:,:,fftpool))!,3*lnzf)
       CALL zeroing(lrxpool(:,:,fftpool))!,2*maxcpu+2)
       CALL zeroing(spm(:,:,fftpool))!,9*maxcpu+9)
       IF (lmsqmax.GT.0) THEN
          CALL zeroing(msqfpool(:,:,fftpool))!,lmsqmax*parai%nproc)
          CALL zeroing(msqspool(:,:,fftpool))!,lmsqmax*parai%nproc)
       ENDIF
       fftpool=fftpool-1
    ELSE
       CALL stopgm("RMFFTNSET","FFTPOOL NOT DEFINED",& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rmfftnset
  ! ==================================================================
  SUBROUTINE addfftnset(ecutf,ecuts,ipool)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ecutf, ecuts
    INTEGER                                  :: ipool

    CHARACTER(*), PARAMETER                  :: procedureN = 'addfftnset'

    INTEGER                                  :: i, ierr, ig, j, k, l, lh1, &
                                                lh2, lh3, nh1, nh2, nh3
    INTEGER, SAVE                            :: icount = 0
    REAL(real_8)                             :: aa1, aa2, aa3, rr, xpaim, &
                                                xplanes, xpnow

! ==--------------------------------------------------------------==

    IF (icount.EQ.0) THEN
       lmsqmax=nhrm
       CALL zeroing(lrxpool)!,2*fftpoolsize*(maxcpu+1))
       CALL zeroing(spm)!,9*fftpoolsize*(maxcpu+1))
       l=(lmsqmax*parai%nproc)
       ALLOCATE(msqf(l),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(msqs(l),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       l=(lmsqmax*parai%nproc*fftpoolsize)
       ALLOCATE(msqfpool(lmsqmax,parai%nproc,fftpoolsize),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(msqspool(lmsqmax,parai%nproc,fftpoolsize),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       l=(ncpw%nhg*fftpoolsize)
       ALLOCATE(nzffp(lnzf,l/lnzf),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (lnzs > 0) THEN
          ALLOCATE(nzfsp(lnzs,l/lnzs),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(nzfsp(1,l),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       ALLOCATE(inzfp(lnzf,l/lnzf),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (lnzs > 0) THEN
          ALLOCATE(inzsp(lnzs,l/lnzs),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(inzsp(1,l),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       l=(3*ncpw%nhg*fftpoolsize)
       ALLOCATE(inzhp(3,lnzf,l/(3*lnzf)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    icount=1
    IF (ecutf.LT.0._real_8 .OR. ecuts.LT.0._real_8) RETURN
    ! 
    ! IF(TKPNT) CALL STOPGM("ADDFFTNSET","NOT IMPLEMENTED")
    ! 
    fftpool=fftpool+1
    IF (fftpool.GT.fftpoolsize) THEN
       CALL stopgm("ADDFFTNSET","TOO MANY ENTRIES IN POOL",& 
            __LINE__,__FILE__)
    ENDIF
    ipool=fftpool
    ! 
    jhg=ncpw%nhg
    DO ig=2,ncpw%nhg
       IF (parm%tpiba2*hg(ig).GT.ecutf) THEN
          jhg=ig-1
          GOTO 100
       ENDIF
    ENDDO
100 CONTINUE
    rr=REAL(jhg,kind=real_8)
    CALL mp_sum(rr,parai%allgrp)
    jhgs=NINT(rr)
    IF (jhgs.GT.spar%nhgs) CALL stopgm("ADDFFTNSET","JHGS TOO LARGE",& 
         __LINE__,__FILE__)
    jgw=ncpw%nhg
    DO ig=2,ncpw%nhg
       IF (parm%tpiba2*hg(ig).GT.ecuts) THEN
          jgw=ig-1
          GOTO 101
       ENDIF
    ENDDO
101 CONTINUE
    rr=REAL(jgw,kind=real_8)
    CALL mp_sum(rr,parai%allgrp)
    jgws=NINT(rr)
    IF (jgws.GT.spar%ngws) CALL stopgm("ADDFFTNSET","JGWS TOO LARGE",& 
         __LINE__,__FILE__)
    ! 
    aa1=parm%alat
    aa2=parm%alat*cell_com%celldm(2)
    aa3=parm%alat*cell_com%celldm(3)
    lr1s=NINT(aa1/pi*SQRT(ecutf)+0.5_real_8)
    lr2s=NINT(aa2/pi*SQRT(ecutf)+0.5_real_8)
    lr3s=NINT(aa3/pi*SQRT(ecutf)+0.5_real_8)
    ! 
    lr1s=fftchk(lr1s,2)
    lr2s=fftchk(lr2s,2)
    lr3s=fftchk(lr3s,2)
    CALL leadim(lr1s,lr2s,lr3s,qr1s,qr2s,qr3s)
    ! 
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    lh1=lr1s/2+1
    lh2=lr2s/2+1
    lh3=lr3s/2+1
    !$omp parallel do private(IG,I,J,K) shared(NH1,NH2,NH3,LH1,LH2,LH3)
    DO ig=1,jhg
       i=inyh(1,ig)-nh1
       j=inyh(2,ig)-nh2
       k=inyh(3,ig)-nh3
       inzhp(1,ig,ipool)=lh1+i
       inzhp(2,ig,ipool)=lh2+j
       inzhp(3,ig,ipool)=lh3+k
    ENDDO
    inzh => inzhp(:,:,ipool)
    ! 
    CALL zeroing(lrxpool(:,:,ipool))!,2*(maxcpu+1))
    xplanes=REAL(lr1s,kind=real_8)
    xpnow=0.0_real_8
    DO i=parai%nproc,1,-1
       xpaim = xpnow + xplanes/parai%nproc
       lrxpool(i-1,1,ipool)=NINT(xpnow)+1
       lrxpool(i-1,2,ipool)=NINT(xpaim)
       IF (NINT(xpaim).GT.lr1s) lrxpool(i-1,2,ipool)=lr1s
       IF (i.EQ.1) lrxpool(i-1,2,ipool)=lr1s
       xpnow = xpaim
    ENDDO
    lr1=lrxpool(parai%mepos,2,ipool)-lrxpool(parai%mepos,1,ipool)+1
    CALL leadim(lr1,lr2s,lr3s,qr1,qr2s,qr3s)
    lr2=lr2s
    lr3=lr3s
    qr2=qr2s
    qr3=qr3s
    llr1=qr1*qr2*qr3
    ! 
    CALL setrays(ipool)
    ! 
    maxfftn = maxfft
    ! 
    fpoolv( 1,ipool) = msrays
    fpoolv( 2,ipool) = mfrays
    fpoolv( 3,ipool) = llr1
    fpoolv( 4,ipool) = qr1s
    fpoolv( 5,ipool) = qr2s
    fpoolv( 6,ipool) = qr3s
    fpoolv( 7,ipool) = qr1
    fpoolv( 8,ipool) = qr2
    fpoolv( 9,ipool) = qr3
    fpoolv(10,ipool) = lr1s
    fpoolv(11,ipool) = lr2s
    fpoolv(12,ipool) = lr3s
    fpoolv(13,ipool) = lr1
    fpoolv(14,ipool) = lr2
    fpoolv(15,ipool) = lr3
    fpoolv(16,ipool) = qr2max
    fpoolv(17,ipool) = qr2min
    fpoolv(18,ipool) = qr3max
    fpoolv(19,ipool) = qr3min
    fpoolv(20,ipool) = lsrm
    fpoolv(21,ipool) = lfrm
    fpoolv(22,ipool) = lr1m
    fpoolv(23,ipool) = lmsq
    fpoolv(24,ipool) = maxfftn
    fpoolv(25,ipool) = jgw
    fpoolv(26,ipool) = jgws
    fpoolv(27,ipool) = jhg
    fpoolv(28,ipool) = jhgs
    ! 
    IF (paral%io_parent) THEN
       WRITE(6,*)
       WRITE(6,'(A,T50,A,I4)') ' ADD NEW FFT SET ',' SET NUMBER ',&
            ipooL
       WRITE(6,'(A,T51,3I5)') ' REAL SPACE GRID ',lr1s,lr2s,lr3S
       WRITE(6,'(A,T20,A,F6.0,T44,A,I10)') ' SPARSE FFT SETUP: ',&
            'CUTOFF [Ry]:',ecuts,'PLANE WAVES:',jgwS
       WRITE(6,'(A,T20,A,F6.0,T44,A,I10)') ' FULL FFT SETUP  : ',&
            'CUTOFF [Ry]:',ecutf,'PLANE WAVES:',jhgS
       WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE addfftnset
  ! ==================================================================
  SUBROUTINE setrays(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    CHARACTER(*), PARAMETER                  :: procedureN = 'setrays'

    INTEGER :: i, ierr, ig, ij, img, iny1, iny2, iny3, ip, ipro, ixf, j, &
      jgwl, jhgl, jj, jmg, msglen, mxrp, nh1, nh2, nh3, ny1, ny2, ny3, qr1m
    INTEGER, ALLOCATABLE                     :: mq(:), my(:)

! Variables
! ==--------------------------------------------------------------==
! GATHER ARRAY FOR FFT ALONG X
! SPARSITY FOR FFT ALONG Y

    ALLOCATE(mg(qr2s,qr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mz((2*qr3s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(my((2*qr2s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(mg)!,qr2s*qr3s)
    CALL zeroing(mz)!,2*qr3s)
    CALL zeroing(my)!,2*qr2s)
    nh1=lr1s/2+1
    nh2=lr2s/2+1
    nh3=lr3s/2+1
    DO ig=1,jgw
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       mg(ny2,ny3)=mg(ny2,ny3)+1
       mg(iny2,iny3)=mg(iny2,iny3)+1
       my(ny2)=my(ny2)+1
       my(iny2)=my(iny2)+1
       mz(ny3)=mz(ny3)+1
       mz(iny3)=mz(iny3)+1
    ENDDO
    CALL mp_sum(my,my(lr2s+1:),lr2s,parai%allgrp)
    CALL icopy(lr2s,my(lr2s+1),1,my,1)
    CALL mp_sum(mz,mz(lr3s+1:),lr3s,parai%allgrp)
    CALL icopy(lr3s,mz(lr3s+1),1,mz,1)
    qr2min=1
    DO i=1,qr2s
       IF (my(i).NE.0) THEN
          qr2min=i
          GOTO 50
       ENDIF
    ENDDO
50  CONTINUE
    qr2max=qr2s
    DO i=qr2s,1,-1
       IF (my(i).NE.0) THEN
          qr2max=i
          GOTO 51
       ENDIF
    ENDDO
51  CONTINUE
    qr3min=1
    DO i=1,qr3
       IF (mz(i).NE.0) THEN
          qr3min=i
          GOTO 52
       ENDIF
    ENDDO
52  CONTINUE
    qr3max=1
    DO i=qr3,1,-1
       IF (mz(i).NE.0) THEN
          qr3max=i
          GOTO 53
       ENDIF
    ENDDO
53  CONTINUE
    ! ==--------------------------------------------------------------==
    img=0
    DO j=1,qr3s
       DO i=1,qr2s
          IF (mg(i,j).NE.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    msrays=img
    DO ig=jgw+1,jhg
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       jmg=mg(ny2,ny3)
       IF (jmg.EQ.0) mg(ny2,ny3)=-1
       jmg=mg(iny2,iny3)
       IF (jmg.EQ.0) mg(iny2,iny3)=-1
    ENDDO
    DO j=1,qr3s
       DO i=1,qr2s
          IF (mg(i,j).LT.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    mfrays=img
    ! 
    jhgl=0
    jgwl=0
    CALL zeroing(spm(:,:,ipool))!,9*maxcpu+9)
    spm(1,parai%mepos,ipool)=jhg
    spm(2,parai%mepos,ipool)=jhgl
    spm(3,parai%mepos,ipool)=jgw
    spm(4,parai%mepos,ipool)=jgwl
    spm(5,parai%mepos,ipool)=lr1
    spm(6,parai%mepos,ipool)=lr2
    spm(7,parai%mepos,ipool)=lr3
    spm(8,parai%mepos,ipool)=mfrays
    spm(9,parai%mepos,ipool)=msrays
    ! 
    nzff => nzffp(:,ipool)
    inzf => inzfp(:,ipool)
    nzfs => nzfsp(:,ipool)
    inzs => inzsp(:,ipool)
    DO i=0,parai%nproc-1
       ipro=i
       CALL mp_bcast(spm(:,i,ipool),9,ipro,parai%allgrp)
    ENDDO
    ! MAXIMUM OF LR1, MSRAYS AND MFRAYS FOR MP_INDEX
    lr1m = 0
    lfrm = 0
    lsrm = 0
    DO i=0,parai%nproc-1
       lr1m = MAX(lr1m,spm(5,i,ipool))
       lfrm = MAX(lfrm,spm(8,i,ipool))
       lsrm = MAX(lsrm,spm(9,i,ipool))
    ENDDO
    qr1m=MAX(lr1m+MOD(lr1m+1,2),qr1)
    lmsq=MAX(lfrm,lsrm)
    ! SCATTER ARRAY FOR FFT ALONG X
    ALLOCATE(mq(lmsq*2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(mq)!,2*lmsq)
    DO i=1,qr2s
       DO j=1,qr3s
          ij=mg(i,j)
          IF (ij.GT.0) THEN
             mq(ij)=i
             mq(lmsq+ij)=j
          ENDIF
       ENDDO
    ENDDO
    ! CONCATENATE GATHER/SCATTER ARRAYS
    msglen = lfrm * 8/2
    CALL my_concat(mq(1),msqs,msglen,parai%allgrp)
    CALL my_concat(mq(lmsq+1),msqf,msglen,parai%allgrp)
    ! TRANSLATE I,J TO A SINGLE G/S INDEX
    DO ip=0,parai%nproc-1
       mxrp=parap%sparm(8,ip)
       DO ixf=1,mxrp
          jj=ixf+ip*lmsq
          i=msqs(jj)
          j=msqf(jj)
          msqfpool(ixf,ip+1,ipool)=i+(j-1)*qr2s
          IF (ixf.LE.parap%sparm(9,ip))&
               msqspool(ixf,ip+1,ipool)=i+(j-qr3min)*qr2s
       ENDDO
    ENDDO
    ! REDEFINE NZH AND INDZ FOR COMPRESSED STORAGE
    !$omp parallel do private(IG,NY1,NY2,NY3,INY1,INY2,INY3)
    DO ig=1,jhg
       ny1=inzh(1,ig)
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny1=-ny1+2*nh1
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       nzff(ig)=ny1 + (mg(ny2,ny3)-1)*qr1s
       inzf(ig)=iny1 + (mg(iny2,iny3)-1)*qr1s
    ENDDO
    !$omp parallel do private(IG)
    DO ig=1,jgw
       nzfs(ig)=nzff(ig)
       inzs(ig)=inzf(ig)
    ENDDO
    DEALLOCATE(mg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(my,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mq,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE setrays
  ! ==================================================================
  SUBROUTINE fft_new_gdist_batch_setup( tfft, nstate, sendsize, sendsize_rem, spin )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN)  :: nstate
    INTEGER, INTENT(OUT) :: sendsize, sendsize_rem
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: spin(:)

    INTEGER :: ierr, Com_in_locks, sendsize_pot, irun, i, j
    INTEGER, SAVE :: remember_batch = 0
    LOGICAL, SAVE :: first, DEBUG_shared_mem = .false.
    TYPE(C_PTR) :: baseptr( 0:parai%node_nproc-1 )
    INTEGER :: arrayshape(3,4), needed_size(4)
    CHARACTER(*), PARAMETER                  :: procedureN = 'fft_new_gdist_batch_setup'
    COMPLEX(real_8), SAVE, POINTER, CONTIGUOUS   :: Big_Com_Pointer(:,:,:)
    LOGICAL,         SAVE, POINTER, CONTIGUOUS   :: Big_1Log_Pointer(:,:,:)
    LOGICAL,         SAVE, POINTER, CONTIGUOUS   :: Big_2Log_Pointer(:,:,:)
    LOGICAL,         SAVE, POINTER, CONTIGUOUS   :: Big_3Log_Pointer(:,:,:)
    LOGICAL :: war(4)

    IF( parai%cp_nproc .gt. 1 ) tfft%do_comm = .true.
    fft_numbuff = 3
    IF( cntl%krwfn ) fft_numbuff = 2
    IF( .not. ( cntl%overlapp_comm_comp .and. fft_numbatches .gt. 1 ) ) fft_numbuff = 1

    IF( remember_batch .ne. fft_batchsize ) THEN

       remember_batch = fft_batchsize

       IF( ALLOCATED( tfft%map_z2y_wave ) )        DEALLOCATE( tfft%map_z2y_wave )
       ALLOCATE( tfft%map_z2y_wave( tfft%my_nr3p * tfft%nr1w * fpar%kr2s * fft_batchsize, 2 ) )
       CALL Make_z2y_Maps( tfft, tfft%map_z2y_wave(:,1), fft_batchsize, tfft%ir1w, tfft%nsw, tfft%nr1w, tfft%small_chunks(1), tfft%big_chunks(1), tfft%map_z2y_bounds(:,:,1) ) 
       IF( fft_residual .ne. 0 ) THEN
          CALL Make_z2y_Maps( tfft, tfft%map_z2y_wave(:,2), fft_residual, tfft%ir1w, tfft%nsw, tfft%nr1w, tfft%small_chunks(1), tfft%big_chunks(1) )
       END IF

       sendsize     = MAXVAL( tfft%nr3p ) * MAXVAL ( tfft%nsw ) * fft_batchsize
       sendsize_rem = MAXVAL( tfft%nr3p ) * MAXVAL ( tfft%nsw ) * fft_residual
       sendsize_pot = MAXVAL( tfft%nr3p ) * MAXVAL(  tfft%nsp )

       arrayshape(1,1) = MAX( sendsize*parai%cp_nproc, sendsize_pot*parai%cp_nproc )
       arrayshape(2,1) = fft_numbuff
       arrayshape(3,1) = 2
       IF( associated( Big_Com_Pointer ) ) DEALLOCATE( Big_Com_Pointer )
       ALLOCATE( Big_Com_Pointer( arrayshape(1,1), arrayshape(2,1), arrayshape(3,1) ) )
       comm_send => Big_Com_Pointer(:,:,1)
       comm_recv => Big_Com_Pointer(:,:,2)

       IF( associated( locks_cc_invfw ) ) DEALLOCATE( locks_cc_invfw )
       ALLOCATE( locks_cc_invfw( 1, ( nstate / fft_batchsize ) + 1, 4 ) )

       IF( allocated( locks_omp ) ) DEALLOCATE( locks_omp )
       ALLOCATE( locks_omp( parai%ncpus_FFT, fft_numbatches+3, 20 ) )

       IF( allocated( locks_omp_big ) ) DEALLOCATE( locks_omp_big )
       ALLOCATE( locks_omp_big( parai%ncpus_FFT, fft_batchsize, fft_numbatches+3, 20 ) )

       CALL Make_Manual_Maps( tfft, fft_batchsize, fft_residual, tfft%nsw, tfft%nr1w, tfft%ngw, tfft%which, nstate )

       first = .true.

    END IF

    IF( first .or. .not. allocated( spin ) ) THEN
       IF( allocated( spin ) ) DEALLOCATE( spin )
       ALLOCATE( spin( 2 ), STAT=ierr )
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       first = .false.
    END IF


  END SUBROUTINE fft_new_gdist_batch_setup

  SUBROUTINE fft_new_gdist_setup( tfft, nss, nr1s, ngs )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: ngs
    INTEGER, INTENT(IN) :: nss(:), nr1s

    INTEGER :: sendsize
    LOGICAL, SAVE :: first = .true.
    TYPE(C_PTR) :: baseptr( 0:parai%node_nproc-1 )
    INTEGER :: arrayshape(3)           
    COMPLEX(real_8), SAVE, POINTER, CONTIGUOUS   :: Big_Pointer(:,:,:)

    IF( parai%cp_nproc .gt. 1 ) tfft%do_comm = .true.
    IF( first .and. .not. restart1%rwf ) THEN

       first = .false.

       arrayshape(2) = 1
       arrayshape(3) = 2
       sendsize = MAXVAL ( tfft%nr3p ) * MAXVAL( nss )
       arrayshape(1) = sendsize*parai%cp_nproc
       ALLOCATE( Big_Pointer( arrayshape(1), arrayshape(2), arrayshape(3) ) )
       comm_send => Big_Pointer(:,:,1)
       comm_recv => Big_Pointer(:,:,2)

       CALL Make_Manual_Maps( tfft, 1, 0, nss, nr1s, ngs, tfft%which, 0 )

       IF( .not. allocated( locks_omp ) ) ALLOCATE( locks_omp( parai%ncpus_FFT, 1, 20 ) )
       !$ locks_omp = .true.
       IF( .not. allocated( locks_omp_big ) ) ALLOCATE( locks_omp_big( parai%ncpus_FFT, 1, 1, 20 ) )
       !$ locks_omp_big = .true.

    ELSE IF( first .and. restart1%rwf ) THEN

       first = .false.

       CALL Make_Manual_Maps( tfft, 1, 0, nss, nr1s, ngs, tfft%which, 0 )

    END IF

  END SUBROUTINE fft_new_gdist_setup

  ! ==================================================================
  SUBROUTINE Make_Manual_Maps( tfft, batch_size, rem_size, nss, nr1s, ngs, which, nstate )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: batch_size, rem_size, nstate
    INTEGER, INTENT(IN) :: nss(:), nr1s, ngs

    CHARACTER(*), PARAMETER                     :: procedureN = 'Make_Manual_Maps'

    INTEGER :: i, j, overlap_cor, eff_nthreads, ierr, which

    ierr = 0

    IF( cntl%overlapp_comm_comp .and. parai%ncpus_FFT .gt. 1 .and. tfft%do_comm(which) .and. which .eq. 1 ) THEN
       eff_nthreads = parai%ncpus_FFT - 1
    ELSE
       eff_nthreads = parai%ncpus_FFT
    END IF
    overlap_cor = parai%ncpus_FFT - eff_nthreads

  ! z things

    DO j = 1, parai%nproc

       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_z_sticks( i, 1, j, which ) = ( nss( j ) * batch_size ) / eff_nthreads
       ENDDO
       DO i = 1+overlap_cor, mod( nss( j ) * batch_size, eff_nthreads ) + overlap_cor
          tfft%thread_z_sticks( i, 1, j, which ) = tfft%thread_z_sticks( i, 1, j, which ) + 1
       ENDDO

       tfft%thread_z_start( 1+overlap_cor, 1, j, which ) = 1
       DO i = 2+overlap_cor, parai%ncpus_FFT
          tfft%thread_z_start( i, 1, j, which ) = tfft%thread_z_start( i-1, 1, j, which ) + tfft%thread_z_sticks( i-1, 1, j, which )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_z_end( i, 1, j, which ) = tfft%thread_z_start( i, 1, j, which ) + tfft%thread_z_sticks( i, 1, j, which ) - 1
       ENDDO

       IF( rem_size .ne. 0 ) THEN
          DO i = 1+overlap_cor, parai%ncpus_FFT
             tfft%thread_z_sticks( i, 2, j, which ) = ( nss( j ) * rem_size ) / eff_nthreads
          ENDDO
          DO i = 1+overlap_cor, mod( nss( j ) * rem_size, eff_nthreads ) + overlap_cor
             tfft%thread_z_sticks( i, 2, j, which ) = tfft%thread_z_sticks( i, 2, j, which ) + 1
          ENDDO

          tfft%thread_z_start( 1+overlap_cor, 2, j, which ) = 1
          DO i = 2+overlap_cor, parai%ncpus_FFT
             tfft%thread_z_start( i, 2, j, which ) = tfft%thread_z_start( i-1, 2, j, which ) + tfft%thread_z_sticks( i-1, 2, j, which )
          ENDDO
          DO i = 1+overlap_cor, parai%ncpus_FFT
             tfft%thread_z_end( i, 2, j, which ) = tfft%thread_z_start( i, 2, j, which ) + tfft%thread_z_sticks( i, 2, j, which ) - 1
          ENDDO
       END IF

    ENDDO

    DO j = 1, parai%nproc

       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_z_sticks( i, 3, j, which ) = ( nss( j ) ) / eff_nthreads
       ENDDO
       DO i = 1+overlap_cor, mod( nss( j ), eff_nthreads ) + overlap_cor
          tfft%thread_z_sticks( i, 3, j, which ) = tfft%thread_z_sticks( i, 3, j, which ) + 1
       ENDDO

       tfft%thread_z_start( 1+overlap_cor, 3, j, which ) = 1
       DO i = 2+overlap_cor, parai%ncpus_FFT
          tfft%thread_z_start( i, 3, j, which ) = tfft%thread_z_start( i-1, 3, j, which ) + tfft%thread_z_sticks( i-1, 3, j, which )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_z_end( i, 3, j, which ) = tfft%thread_z_start( i, 3, j, which ) + tfft%thread_z_sticks( i, 3, j, which ) - 1
       ENDDO

    ENDDO

    IF( mod( nstate, 2 ) .ne. 0 ) THEN

       IF( rem_size .eq. 0 ) THEN

          IF( batch_size .gt. 1 ) THEN

             DO i = 1+overlap_cor, parai%ncpus_FFT
                tfft%thread_prepare_sticks( i, 1 ) = ( nss( parai%me+1 ) * (batch_size-1) ) / eff_nthreads
             ENDDO
             DO i = 1+overlap_cor, mod( nss( parai%me+1 ) * (batch_size-1), eff_nthreads ) + overlap_cor
                tfft%thread_prepare_sticks( i, 1 ) = tfft%thread_prepare_sticks( i, 1 ) + 1
             ENDDO

             tfft%thread_prepare_start( 1+overlap_cor, 1 ) = 1
             DO i = 2+overlap_cor, parai%ncpus_FFT
                tfft%thread_prepare_start( i, 1 ) = tfft%thread_prepare_start( i-1, 1 ) + tfft%thread_prepare_sticks( i-1, 1 )
             ENDDO
             DO i = 1+overlap_cor, parai%ncpus_FFT
                tfft%thread_prepare_end( i, 1 ) = tfft%thread_prepare_start( i, 1 ) + tfft%thread_prepare_sticks( i, 1 ) - 1
             ENDDO

          ELSE

             tfft%thread_prepare_sticks(:,1) = 0
             tfft%thread_prepare_start(:,1)  = 1
             tfft%thread_prepare_end(:,1)    = 0

          END IF

       ELSE

          IF( rem_size .gt. 1 ) THEN

             DO i = 1+overlap_cor, parai%ncpus_FFT
                tfft%thread_prepare_sticks( i, 1 ) = ( nss( parai%me+1 ) * (rem_size-1) ) / eff_nthreads
             ENDDO
             DO i = 1+overlap_cor, mod( nss( parai%me+1 ) * (rem_size-1), eff_nthreads ) + overlap_cor
                tfft%thread_prepare_sticks( i, 1 ) = tfft%thread_prepare_sticks( i, 1 ) + 1
             ENDDO

             tfft%thread_prepare_start( 1+overlap_cor, 1 ) = 1
             DO i = 2+overlap_cor, parai%ncpus_FFT
                tfft%thread_prepare_start( i, 1 ) = tfft%thread_prepare_start( i-1, 1 ) + tfft%thread_prepare_sticks( i-1, 1 )
             ENDDO
             DO i = 1+overlap_cor, parai%ncpus_FFT
                tfft%thread_prepare_end( i, 1 ) = tfft%thread_prepare_start( i, 1 ) + tfft%thread_prepare_sticks( i, 1 ) - 1
             ENDDO

          ELSE

             tfft%thread_prepare_sticks(:,1) = 0
             tfft%thread_prepare_start(:,1)  = 1
             tfft%thread_prepare_end(:,1)    = 0

          END IF

       END IF

       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_prepare_sticks( i, 2 ) = nss( parai%me+1 ) / eff_nthreads
       ENDDO
       DO i = 1+overlap_cor, mod( nss( parai%me+1 ) , eff_nthreads ) + overlap_cor
          tfft%thread_prepare_sticks( i, 2 ) = tfft%thread_prepare_sticks( i, 2 ) + 1
       ENDDO

       tfft%thread_prepare_start( 1+overlap_cor, 2 ) = tfft%thread_prepare_end( parai%ncpus_FFT, 1 ) + 1
       DO i = 2+overlap_cor, parai%ncpus_FFT
          tfft%thread_prepare_start( i, 2 ) = tfft%thread_prepare_start( i-1, 2 ) + tfft%thread_prepare_sticks( i-1, 2 )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_prepare_end( i, 2 ) = tfft%thread_prepare_start( i, 2 ) + tfft%thread_prepare_sticks( i, 2 ) - 1
       ENDDO

    END IF

  ! y things

    DO i = 1+overlap_cor, parai%ncpus_FFT
       tfft%thread_y_sticks( i, which ) = ( nr1s * tfft%my_nr3p ) / eff_nthreads
    ENDDO
    DO i = 1+overlap_cor, mod( nr1s * tfft%my_nr3p, eff_nthreads ) + overlap_cor
       tfft%thread_y_sticks( i, which ) = tfft%thread_y_sticks( i, which ) + 1
    ENDDO

    tfft%thread_y_start( 1+overlap_cor, which ) = 1
    DO i = 2+overlap_cor, parai%ncpus_FFT
       tfft%thread_y_start( i, which ) = tfft%thread_y_start( i-1, which ) + tfft%thread_y_sticks( i-1, which )
    ENDDO
    DO i = 1+overlap_cor, parai%ncpus_FFT
       tfft%thread_y_end( i, which ) = tfft%thread_y_start( i, which ) + tfft%thread_y_sticks( i, which ) - 1
    ENDDO

  ! x things

    DO i = 1+overlap_cor, parai%ncpus_FFT
       tfft%thread_x_sticks( i, which ) = ( tfft%my_nr3p * fpar%kr2s ) / eff_nthreads
    ENDDO
    DO i = 1+overlap_cor, mod( tfft%my_nr3p * fpar%kr2s, eff_nthreads ) + overlap_cor
       tfft%thread_x_sticks( i, which ) = tfft%thread_x_sticks( i, which ) + 1
    ENDDO

    tfft%thread_x_start( 1+overlap_cor, which ) = 1
    DO i = 2+overlap_cor, parai%ncpus_FFT
       tfft%thread_x_start( i, which ) = tfft%thread_x_start( i-1, which ) + tfft%thread_x_sticks( i-1, which )
    ENDDO
    DO i = 1+overlap_cor, parai%ncpus_FFT
       tfft%thread_x_end( i, which ) = tfft%thread_x_start( i, which ) + tfft%thread_x_sticks( i, which ) - 1
    ENDDO

    IF( tfft%which .eq. 1 ) THEN

     ! gspace things

       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_ngms( i ) = ( ngs ) / eff_nthreads
       ENDDO
       DO i = 1+overlap_cor, mod( ngs, eff_nthreads ) + overlap_cor
          tfft%thread_ngms( i ) = tfft%thread_ngms( i ) + 1
       ENDDO

       tfft%thread_ngms_start( 1+overlap_cor ) = 1
       DO i = 2+overlap_cor, parai%ncpus_FFT
          tfft%thread_ngms_start( i ) = tfft%thread_ngms_start( i-1 ) + tfft%thread_ngms( i-1 )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_ngms_end( i ) = tfft%thread_ngms_start( i ) + tfft%thread_ngms( i ) - 1
       ENDDO

     ! rspace things

       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_rspace( i ) =  tfft%thread_x_sticks( i, 1 ) * fpar%kr1s
       ENDDO

       tfft%thread_rspace_start( 1+overlap_cor ) = 1
       DO i = 2+overlap_cor, parai%ncpus_FFT
          tfft%thread_rspace_start( i ) = tfft%thread_rspace_start( i-1 ) + tfft%thread_rspace( i-1 )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus_FFT
          tfft%thread_rspace_end( i ) = tfft%thread_rspace_start( i ) + tfft%thread_rspace( i ) - 1
       ENDDO

    END IF

  END SUBROUTINE Make_Manual_Maps
  ! ==================================================================
  SUBROUTINE Make_z2y_Maps( tfft, map_z2y, batch_size, ir1s, nss, my_nr1s, small_chunks, big_chunks, map_z2y_bounds )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN)  :: batch_size, my_nr1s, small_chunks, big_chunks
    INTEGER, INTENT(OUT) :: map_z2y(:)
    INTEGER, INTENT(IN)  :: ir1s(:), nss(:)
    INTEGER, OPTIONAL, INTENT(OUT) :: map_z2y_bounds(:,:)

    LOGICAL :: l_map( tfft%my_nr3p * my_nr1s * fpar%kr2s )
    LOGICAL :: first, start
    INTEGER :: j, l, i, k
    INTEGER :: offset, m, m1, m2, ip, pos
    INTEGER :: ierr

    l_map = .false.
    map_z2y = 0

    !$omp parallel private( j, l, ip, offset, i, m, m1, m2, pos, k)
    ip = 0
    DO j = 1, parai%nnode
       DO l = 1, parai%node_nproc_overview( j )
          ip = ip + 1
          offset = (ip-1) * batch_size * small_chunks
          !$omp do
          DO i = 1, nss( ip )
             m = tfft%ismap( i + tfft%iss(ip) ) !number of current pencil
             m1 = mod ( m-1, fpar%kr1s ) + 1     !coordinate of pencil
             m2 = (m-1)/fpar%kr1s + 1            !other coordinate of pencil
             pos = m2 + ( ir1s(m1) - 1 ) * fpar%kr2s
             DO k = 1, tfft%my_nr3p
                l_map( pos ) = .true.
                map_z2y( pos ) = k + offset + tfft%nr3px * (i-1)
                pos = pos + fpar%kr2s * my_nr1s
             ENDDO
          ENDDO
          !$omp end do
       ENDDO
    ENDDO
    !$omp end parallel

    IF( present( map_z2y_bounds ) ) THEN

       map_z2y_bounds = 0

       !$omp parallel do private( i, start, j, k )
       DO i = 1, my_nr1s
          start = .true.
          k = 1
          DO j = 1, fpar%kr2s

             IF( j .eq. 1 ) THEN
                IF( l_map( (i-1)*fpar%kr2s + j ) .eqv. .true. ) THEN
                   k = k + 1
                END IF
             END IF
             IF( l_map( (i-1)*fpar%kr2s + j ) .eqv. .false. ) THEN
                IF( start .eqv. .true. ) THEN
                   map_z2y_bounds(i,2*k-1) = j
                   start = .false.
                END IF
             ELSE
                IF( start .eqv. .false. ) THEN
                   map_z2y_bounds(i,2*k) = j
                   k = k + 1
                   start = .true.
                END IF
             END IF
             IF( j .eq. fpar%kr2s ) THEN
                IF( start .eqv. .false. ) THEN
                   map_z2y_bounds(i,2*k) = j + 1
                END IF
             END IF

          ENDDO
       ENDDO
       !$omp end parallel do

       !$omp parallel do private( i, start, j, k )
       DO i = 1, my_nr1s
          start = .true.
          k = 3
          DO j = 1, fpar%kr2s

             IF( j .eq. 1 ) THEN
                IF( l_map( (i-1)*fpar%kr2s + j ) .eqv. .false. ) THEN
                   k = k + 1
                END IF
             END IF
             IF( l_map( (i-1)*fpar%kr2s + j ) .eqv. .true. ) THEN
                IF( start .eqv. .true. ) THEN
                   map_z2y_bounds(i,2*k-1) = j
                   start = .false.
                END IF
             ELSE
                IF( start .eqv. .false. ) THEN
                   map_z2y_bounds(i,2*k) = j
                   k = k + 1
                   start = .true.
                END IF
             END IF
             IF( j .eq. fpar%kr2s ) THEN
                IF( start .eqv. .false. ) THEN
                   map_z2y_bounds(i,2*k) = j + 1
                END IF
             END IF

          ENDDO
       ENDDO
       !$omp end parallel do

    END IF

  END SUBROUTINE Make_z2y_Maps
  ! ==================================================================
  SUBROUTINE Pre_Initialize_C2_Com( c2, nstate, las, nstate_local, my_start, c2_com_num, c2_com_recv, fft_batchsize, fft_residual, fft_numbatches, cp_nstates, s4_coms, my_ngw )
    IMPLICIT NONE

    COMPLEX(real_8), CONTIGUOUS, INTENT(IN) :: c2(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: c2_com_num(:,:), c2_com_recv(:,:)
    INTEGER, INTENT(IN) :: nstate, las, nstate_local, my_start, my_ngw
    INTEGER, INTENT(IN) :: fft_batchsize, fft_residual, fft_numbatches
    INTEGER, INTENT(OUT):: cp_nstates(:)

    INTEGER :: ngw( parai%nproc, parai%cp_nogrp )
    INTEGER :: start( parai%nproc, parai%cp_nogrp )
    INTEGER :: nstate_end( parai%cp_nogrp )
    INTEGER :: com_matrix( parai%nproc, parai%cp_nogrp )
    INTEGER :: counter, current, i, j, k, jter, kter, bsize, remove, next
    INTEGER :: kter_save( fft_numbatches+1 )
    INTEGER :: fft_info( 4, parai%cp_nogrp ), which_batch( parai%cp_nogrp )

    INTEGER, ALLOCATABLE :: baseoffset(:), state_per_Com(:,:), offset_per_Com(:,:), s4_coms(:)
    INTEGER :: inb, icp

    ngw = 0
    start = 0
    nstate_end = 0
    com_matrix = 0
    fft_info = 0

    remove = 0
    IF( mod(nstate_local,2) .ne. 0 ) remove = 1

    IF( parai%me .eq. 0 ) THEN
       fft_info( 1, parai%cp_inter_me+1 ) = fft_numbatches
       fft_info( 2, parai%cp_inter_me+1 ) = fft_batchsize
       fft_info( 3, parai%cp_inter_me+1 ) = fft_residual
       fft_info( 4, parai%cp_inter_me+1 ) = remove
    END IF

    com_matrix( parai%me+1, parai%cp_inter_me+1 ) = parai%cp_me
    IF( parai%me .eq. 0 ) nstate_end( parai%cp_inter_me+1 ) = las

    CALL cp_grp_get_sizes( ngw_l=ngw( parai%me+1, parai%cp_inter_me+1 ), first_g=start( parai%me+1, parai%cp_inter_me+1 ) )

    CALL MP_SUM( ngw         , parai%nproc*parai%cp_nogrp, parai%cp_grp )
    CALL MP_SUM( start       , parai%nproc*parai%cp_nogrp, parai%cp_grp )
    CALL MP_SUM( nstate_end  ,             parai%cp_nogrp, parai%cp_grp )
    CALL MP_SUM( com_matrix  , parai%nproc*parai%cp_nogrp, parai%cp_grp )
    CALL MP_SUM( fft_info    ,           4*parai%cp_nogrp, parai%cp_grp )

    IF( allocated( parai%c2_send_handle ) ) DEALLOCATE( parai%c2_send_handle )
    ALLOCATE( parai%c2_send_handle( fft_batchsize*2*(parai%cp_nogrp-1) ,fft_numbatches+1 ) )
   
    IF( .not. allocated( parai%c2_send_handle_counter ) ) ALLOCATE( parai%c2_send_handle_counter( nstate ) )
    IF( .not. allocated( parai%c2_recv_handle_counter ) ) ALLOCATE( parai%c2_recv_handle_counter( nstate ) )
   
    IF( allocated( c2_com_num ) ) DEALLOCATE( c2_com_num )
    ALLOCATE( c2_com_num( fft_numbatches+1, 4 ) )
    
    c2_com_num = 0
   
    IF( allocated( parai%c2_recv_handle ) ) DEALLOCATE( parai%c2_recv_handle )
    ALLOCATE( parai%c2_recv_handle( MAXVAL( fft_info( 2, : ) )*2*parai%cp_nogrp , MAXVAL( fft_info( 1, : ) )+1 ) )
   
    counter = 0
    kter = 0
    kter_save = 0
    DO i = 1, parai%cp_nogrp
       IF( i .eq. parai%cp_inter_me+1 ) CYCLE
       jter = 0
       DO j = 1, fft_numbatches+1
          bsize = fft_batchsize*2
          kter = kter_save(j)
          IF( j .eq. fft_numbatches+1 ) bsize = fft_residual*2
          IF( bsize .eq. 0 ) CYCLE
          IF( ( fft_residual .eq. 0 .and. j .eq. fft_numbatches ) .or. &
              ( j .eq. fft_numbatches+1 ) ) bsize = bsize - remove
          DO k = 1, bsize
             kter = kter + 1
             jter = jter + 1
             counter = counter + 1
             CALL mp_send_init_complex( c2(:,my_start+jter), start( parai%me+1, i )-1, ngw( parai%me+1, i ), com_matrix( parai%me+1, i ), parai%cp_me, parai%cp_grp, parai%c2_send_handle( kter, j ) )
             parai%c2_send_handle_counter( counter ) = parai%c2_send_handle( kter, j )
          ENDDO
          kter_save(j) = kter
       ENDDO
    ENDDO
    c2_com_num(1,2) = counter
    c2_com_num(:,1) = fft_batchsize*2
    IF( fft_residual .ne. 0 ) THEN
       c2_com_num(fft_numbatches+1,1) = fft_residual*2 - remove
    ELSE
       c2_com_num(fft_numbatches,1) = c2_com_num(fft_numbatches,1) - remove
    END IF
    c2_com_num(:,1) = c2_com_num(:,1) * (parai%cp_nogrp-1)
   
    counter = 0
    current = 1
    next = 0
    which_batch = 1
    DO i = 1, nstate
       IF( i .gt. nstate_end(current) ) THEN
          current = current + 1
          next = 1
       ELSE
          next = next + 1
       END IF
       IF( current .eq. parai%cp_inter_me+1 ) CYCLE
       counter = counter + 1
       c2_com_num( which_batch( current ), 3 ) = c2_com_num( which_batch( current ), 3 ) + 1
       IF( next .gt. fft_info( 2, current )*2 ) THEN
          next = 1
          c2_com_num( which_batch( current ), 3 ) = c2_com_num( which_batch( current ), 3 ) - 1
          which_batch( current ) = which_batch( current ) + 1
          c2_com_num( which_batch( current ), 3 ) = c2_com_num( which_batch( current ), 3 ) + 1
       END IF
   
   
       CALL mp_recv_init_complex( c2(:,i), start( parai%me+1, parai%cp_inter_me+1 )-1, ngw( parai%me+1, parai%cp_inter_me+1 ), &
                                  com_matrix( parai%me+1, current ), parai%cp_grp, parai%c2_recv_handle_counter( counter ) )
       parai%c2_recv_handle( c2_com_num( which_batch( current ), 3 ), which_batch( current ) ) = parai%c2_recv_handle_counter( counter )
   
    ENDDO
    c2_com_num(2,2) = counter

    IF( allocated( parai%c2_comb_handle ) ) DEALLOCATE( parai%c2_comb_handle )
    ALLOCATE( parai%c2_comb_handle( c2_com_num(1,1)+c2_com_num(1,3) , fft_numbatches+1  ) )
    
    DO i = 1, fft_numbatches+1
       c2_com_num(i,4) = c2_com_num(i,1) + c2_com_num(i,3)
       DO j = 1, c2_com_num(i,1)
          parai%c2_comb_handle( j, i ) = parai%c2_send_handle( j, i )
       ENDDO
       DO j = 1, c2_com_num(i,3)
          parai%c2_comb_handle( j + c2_com_num(i,1), i ) = parai%c2_recv_handle( j, i )
       ENDDO
    ENDDO

    IF( allocated( parai%c2_comb_handle_counter ) ) DEALLOCATE( parai%c2_comb_handle_counter )
    ALLOCATE( parai%c2_comb_handle_counter( nstate * 2 ) )
    
    c2_com_num(3,2) = c2_com_num(1,2) + c2_com_num(2,2)
    DO j = 1, c2_com_num(1,2)
       parai%c2_comb_handle_counter( j ) = parai%c2_send_handle_counter( j )
    ENDDO
    DO j = 1, c2_com_num(2,2)
       parai%c2_comb_handle_counter( j + c2_com_num(1,2) ) = parai%c2_recv_handle_counter( j )
    ENDDO

  END SUBROUTINE Pre_Initialize_C2_Com

END MODULE fftnew_utils
