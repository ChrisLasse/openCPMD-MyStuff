#include "cpmd_global.h"

MODULE loadpa_utils
  USE cppt,                            ONLY: hg,&
                                             inyh
  USE distribution_utils,              ONLY: dist_entity2,&
                                             dist_entity
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: tfft,&
                                             FFT_TYPE_DESCRIPTOR
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: epsg,&
                                             epsgx,&
                                             gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE mp_interface,                    ONLY: mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE prmem_utils,                     ONLY: prmem
  USE sort_utils,                      ONLY: sort2
  USE sphe,                            ONLY: gcutwmax
  USE system,                          ONLY: &
       fpar, iatpe, iatpe_cp, iatpt, ipept, ipept_cp, natpe_cp,  &
       mapgp, natpe, ncpw, nkpt, norbpe, parap, parm, spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: loadpa
  PUBLIC :: leadim

CONTAINS

  ! ==================================================================
  SUBROUTINE loadpa
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'loadpa'

    INTEGER :: i, i0, ia, iat, icpu, ierr, ig, ihrays, ii, img, in1, in2, &
      in3, iorb, ip, ipp, ir, is, isub, isub2, isub3, isub4, izrays, izpl, j, &
      j1, j2, jmax, jmin, k, kmax, kmin, mspace, nh1, nh2, nh3, nthreads,&
      first, last, offset
    INTEGER, ALLOCATABLE                     :: ihray(:,:), ixray(:,:), &
                                                mgpa(:,:), ind(:), pre_inyh(:,:)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: thread_buff
    LOGICAL                                  :: oldstatus
    REAL(real_8)                             :: g2, sign, t
    REAL(real_8), PARAMETER                  :: eps8=1.0E-8_real_8

! ==--------------------------------------------------------------==
! ==  DISTRIBUTION OF PARALLEL WORK                               ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    nthreads = 1
    !$    nthreads = omp_get_max_threads( )
    ALLOCATE(thread_buff(2*nthreads),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (paral%io_parent)WRITE(6,'(/," ",16("PARA"))')
    mspace = (fpar%kr2s*fpar%kr3s)/2 + 1
    ALLOCATE(ixray(fpar%kr1s,fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ihray(fpar%kr1s,fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mgpa(fpar%kr2s,fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(iatpe(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(iatpe_cp(ions1%nat,0:parai%cp_nogrp-1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ipept_cp(2,0:parai%nproc-1,0:parai%cp_nogrp-1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL mp_sync(parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE ATOMS
    ! ==--------------------------------------------------------------==
    CALL dist_entity(ions1%nat,parai%nproc,ipept)

    CALL mm_dim(mm_go_mm,oldstatus)
    ALLOCATE(iatpt(2,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          iatpt(1,iat)=ia
          iatpt(2,iat)=is
       ENDDO
    ENDDO
    CALL mm_dim(mm_revert,oldstatus)

    DO i=0,parai%nproc-1
       DO j=ipept(1,i),ipept(2,i)
          iatpe(j)=i
       ENDDO
    ENDDO
    natpe=ipept(2,parai%mepos)-ipept(1,parai%mepos)+1

    DO i=0,parai%nproc-1
       do j=0,parai%cp_nogrp-1
          CALL part_1d_get_blk_bounds((ipept(2,i)-ipept(1,i)+1),j,parai%cp_nogrp,first,last)
          ipept_cp(1,i,j)=ipept(1,i)+first-1
          ipept_cp(2,i,j)=ipept(1,i)+last-1
       end do
    ENDDO

    iatpe_cp=-1 !mepos .ge. 0

    DO i=0,parai%nproc-1
       do k=0,parai%cp_nogrp-1
          DO j=ipept_cp(1,i,k),ipept_cp(2,i,k)
             iatpe_cp(j,k)=i
          end do
       ENDDO
    ENDDO

    natpe_cp=&
         ipept_cp(2,parai%mepos,parai%cp_inter_me)&
         -ipept_cp(1,parai%mepos,parai%cp_inter_me)+1

    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE ORBITALS
    ! ==--------------------------------------------------------------==
    CALL dist_entity2(crge%n,parai%nproc,parap%nst12,nblocal=norbpe,iloc=parai%me)
    ! ==--------------------------------------------------------------==

!    ALLOCATE(tfft%nr3_ranges(0:parai%nproc-1,2),STAT=ierr)
!    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
!         __LINE__,__FILE__)
    ALLOCATE(tfft%indx_map(fpar%kr1s,fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%indx(fpar%kr1s*fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%ind1(fpar%kr1s*fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%ind2(fpar%kr1s*fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    tfft%indx_map = 0
    tfft%indx = 0
    tfft%ind1 = 0
    tfft%ind2 = 0

    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE REAL SPACE XY-PLANES
    ! ==--------------------------------------------------------------==
!CLR: Currently not used, maybe activated again at a later time
!    CALL dist_entity2(spar%nr3s,parai%nproc,tfft%nr3_ranges)
    CALL zeroing(parap%nrzpl)!,2*(maxcpu+1))
    IF (isos1%tclust.AND.isos3%ps_type.EQ.1) THEN
!CLR: Whatever case this is: it will probably not work if triggered
       ! DISTRIBUTE REAL SPACE XY-PLANES
       CALL dist_entity2(2*spar%nr3s,parai%nproc,parap%nrzpl)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE G-VECTORS AND MARK ASSOCIATED RAYS
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN//'_c',isub4)
    CALL zfft(ixray,gcutwmax)
    CALL zfft(ihray,gvec_com%gcut)

    CALL Distribute_Sticks( fpar%kr1s, fpar%kr2s, ixray )
    CALL Distribute_Sticks( fpar%kr1s, fpar%kr2s, ihray, ixray )

    CALL czfft(ixray,gcutwmax)
    CALL czfft(ihray,gvec_com%gcut)

    CALL tihalt(procedureN//'_c',isub4)
    ! ==--------------------------------------------------------------==
    ! ALLOCATE THE GLOBAL DATA ARRAYS
    ! ==--------------------------------------------------------------==
    ncpw%nhg=INT(REAL(spar%nhgs,kind=real_8)/REAL(parai%nproc,kind=real_8)*1.2_real_8)
    ncpw%nhg=MIN(ncpw%nhg,spar%nhgs)
    ncpw%nhg=MAX(ncpw%nhg,100)
    ncpw%ngw=INT(REAL(spar%ngws,kind=real_8)/REAL(parai%nproc,kind=real_8)*1.2_real_8)
    ncpw%ngw=MIN(ncpw%ngw,spar%ngws)
    ncpw%ngw=MAX(ncpw%ngw,100)

    ALLOCATE(hg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(pre_inyh(3,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(inyh(3,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mapgp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(hg)!,nhg)
    CALL zeroing(inyh)!,3*nhg)
    CALL zeroing(mapgp)!,nhg)
    ! ==--------------------------------------------------------------==
    ! SYSTEM PARAMETER PER CPU
    ! ==--------------------------------------------------------------==
    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    ! ==--------------------------------------------------------------==
    ! TO GET UNIQUE ORDERING OF G-VECTORS, BREAK THE SYMMETRY
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN//'_a',isub2)
    ! EPSG=EPSGX * ACOS(-1._real_8) *real(NHGS,kind=real_8)
    ! EPSG=GCUT*EPSGX
    epsg=SQRT(REAL(spar%nhgs-1,kind=real_8))*gvec_com%gcut*epsgx
    ig=0
    ncpw%ngw=0
    ncpw%nhg=0
    DO i=0,parm%nr3-1
       jmin=-parm%nr1+1
       jmax=parm%nr1-1
       IF (i.EQ.0) THEN
          jmin=0
       ENDIF
       DO j=jmin,jmax
          kmin=-parm%nr2+1
          kmax=parm%nr2-1
          IF (i.EQ.0.AND.j.EQ.0) THEN
             kmin=0
          ENDIF
          DO k=kmin,kmax
             g2=0._real_8
             DO ir=1,3
                t=REAL(i,kind=real_8)*gvec_com%b3(ir)+REAL(j,kind=real_8)*gvec_com%b1(ir)+REAL(k,kind=real_8)*gvec_com%b2(ir)
                g2=g2+t*t
             ENDDO
             IF (compare_lt(g2,gvec_com%gcut)) THEN
!             IF (g2.LT.gvec_com%gcut) THEN
                ig=ig+1
                in1=nh1+i
                in2=nh2+j
                in3=nh3+k
                IF (compare_lt(g2,gcutwmax)) THEN
!                IF (g2.LT.gcutwmax) THEN
                   icpu=ixray(in1,in2)
                   IF (-icpu.EQ.parai%mepos+1) THEN
                      ncpw%ngw=ncpw%ngw+1
                   ENDIF
                   sign=-1._real_8
                ELSE
                   sign=1._real_8
                ENDIF
                icpu=ihray(in1,in2)
                IF (-icpu.EQ.parai%mepos+1) THEN
                   ncpw%nhg=ncpw%nhg+1
                   ! HG(NHG)=G2+SQRT(real(IG-1,kind=real_8))*EPSG*SIGN
                   ! G2*EPSGX*SIGN is the epsilon added (>= epsilon(G2))
                   ! SQRT(FLOAT(IG-1)) is to break the symmetry
                   hg(ncpw%nhg)=g2
                   pre_inyh(1,ncpw%nhg)=in1
                   pre_inyh(2,ncpw%nhg)=in2
                   pre_inyh(3,ncpw%nhg)=in3
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    tfft%nhg = ncpw%nhg
    tfft%ngw = ncpw%ngw

    CALL zeroing(mgpa)!,kr2s*kr3s)
    parai%nhrays=0
    parai%ngrays=0
    DO j1=1,fpar%kr1s
       DO j2=1,fpar%kr2s
          IF (ihray(j2,j1).EQ.-(parai%mepos+1)) parai%nhrays=parai%nhrays+1
          IF (ixray(j2,j1).EQ.-(parai%mepos+1)) THEN
             parai%ngrays=parai%ngrays+1
             mgpa(j2,j1)=parai%ngrays
          ENDIF
       ENDDO
    ENDDO
    img=parai%ngrays
    DO j1=1,fpar%kr1s
       DO j2=1,fpar%kr2s
          IF (ihray(j2,j1).EQ.-(parai%mepos+1).AND.mgpa(j2,j1).EQ.0) THEN
             img=img+1
             mgpa(j2,j1)=img
          ENDIF
       ENDDO
    ENDDO
    ncpw%nhgl=spar%nhgls
    ncpw%ngwl=spar%ngwls
    parm%nr1 =parap%nrxpl(parai%mepos,2)-parap%nrxpl(parai%mepos,1)+1
    parm%nr2 =spar%nr2s
    parm%nr3 =spar%nr3s
    parap%sparm(1,parai%mepos)=ncpw%nhg
    parap%sparm(2,parai%mepos)=ncpw%nhgl
    CALL tihalt(procedureN//'_a',isub2)
    ! ==--------------------------------------------------------------==
    ! IF K POINTS, WE NEED TO DOUBLE DIMENSIONS.
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       nkpt%ngwk=ncpw%ngw*2
       nkpt%nhgk=ncpw%nhg*2
    ELSE
       nkpt%ngwk=ncpw%ngw
       nkpt%nhgk=ncpw%nhg
    ENDIF
    parap%sparm(3,parai%mepos)=ncpw%ngw
    parap%sparm(4,parai%mepos)=ncpw%ngwl
    parap%sparm(5,parai%mepos)=parm%nr1
    parap%sparm(6,parai%mepos)=parm%nr2
    parap%sparm(7,parai%mepos)=parm%nr3
    parap%sparm(8,parai%mepos)=parai%nhrays
    parap%sparm(9,parai%mepos)=parai%ngrays
    CALL my_allgather_i(parap%sparm,SIZE(parap%sparm,1),parai%allgrp)

    ! ==--------------------------------------------------------------==
    ! SORTING OF G-VECTORS
    ! ==--------------------------------------------------------------==

    ALLOCATE(ind(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ind = 0
    CALL hpsort_eps( ncpw%nhg, hg, ind, eps8 )

    DO i = 1, ncpw%nhg
       inyh(:,i) = pre_inyh(:,ind(i))
    ENDDO

    DO ig=1,ncpw%nhg
       mapgp(ig)=ig
    ENDDO
    ! ==--------------------------------------------------------------==
    geq0=.FALSE.
    i0=0
    DO ig=1,ncpw%ngw
       IF (hg(ig).LT.1.e-5_real_8) THEN
          geq0=.TRUE.
          i0=parai%mepos
       ENDIF
    ENDDO
    CALL mp_sum(i0,parai%igeq0,parai%allgrp)
    IF (paral%io_parent) THEN
       WRITE(6,'(13X,A,I5)') '   G=0 COMPONENT ON PROCESSOR : ',parai%igeq0
       WRITE(6,'(" ",16("PARA"),/)')
       CALL prmem(procedureN)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SETUP ARRAYS NEEDED FOR NEW GDISTRIBUTION FFT
    ! ==--------------------------------------------------------------==

    ALLOCATE( parai%cp_overview( 4, parai%nproc ) )
    parai%cp_overview = 0
    parai%cp_overview( 1, parai%me+1 ) = parai%cp_me
    parai%cp_overview( 2, parai%me+1 ) = parai%node_me
    parai%cp_overview( 3, parai%me+1 ) = parai%node_nproc
    Call mp_sum(parai%cp_overview,parai%nproc*4,parai%allgrp)
    offset = 0
    DO i = 1, parai%nproc
       parai%cp_overview( 4, i ) = offset
       IF( parai%cp_overview( 2, i ) + 1 .eq. parai%cp_overview( 3, i ) ) offset = offset + 1
    ENDDO

    parai%nnode = MAXVAL(parai%cp_overview(4,:))+1
    parai%my_node = parai%cp_overview( 4, parai%me+1 )
    parai%max_node_nproc = MAXVAL(parai%cp_overview(3,:))

    ALLOCATE( parai%node_nproc_overview( parai%nnode ) )
    DO i = 1, parai%nnode
       DO j = 1, parai%nproc
          IF( parai%cp_overview(4,j) .eq. i-1 ) parai%node_nproc_overview( i ) = parai%cp_overview(3,j)
       ENDDO
    ENDDO

    ALLOCATE( parai%node_grpindx( parai%node_nproc ) )
    DO i = 1, parai%nproc
       IF( parai%my_node .eq. parai%cp_overview(4,i) ) parai%node_grpindx( parai%cp_overview(2,i)+1 ) = i - 1
    ENDDO

    ALLOCATE( tfft%cp_ngws( parai%cp_nogrp ) )
    ALLOCATE( tfft%cp_nstates( parai%cp_nogrp ) )

    !For now here; need to have a better spot for this eventually!
    parai%ncpus_FFT = parai%ncpus

    ALLOCATE( tfft%thread_z_sticks( parai%ncpus_FFT, 3, parai%nproc, 2 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_z_start( parai%ncpus_FFT, 3, parai%nproc, 2 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_z_end( parai%ncpus_FFT, 3, parai%nproc, 2 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_prepare_sticks( parai%ncpus_FFT, 2 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_prepare_start( parai%ncpus_FFT, 2 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_prepare_end( parai%ncpus_FFT, 2 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_y_sticks( parai%ncpus_FFT, 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_y_start( parai%ncpus_FFT, 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_y_end( parai%ncpus_FFT, 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_x_sticks( parai%ncpus_FFT, 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_x_start( parai%ncpus_FFT, 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_x_end( parai%ncpus_FFT, 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_ngms( parai%ncpus_FFT ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_ngms_start( parai%ncpus_FFT ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_ngms_end( parai%ncpus_FFT ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_rspace( parai%ncpus_FFT ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_rspace_start( parai%ncpus_FFT ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE( tfft%thread_rspace_end( parai%ncpus_FFT ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    CALL SetupArrays( ihray, ixray )

    ! ==--------------------------------------------------------------==
    ! LEADING DIMENSIONS OF REAL SPACE ARRAYS
    ! ==--------------------------------------------------------------==
    CALL leadim(parm%nr1,parm%nr2,parm%nr3,fpar%kr1,fpar%kr2,fpar%kr3)
    fpar%kr1 = tfft%nr3p( parai%me+1 )
    fpar%nnr1=fpar%kr1*fpar%kr2s*fpar%kr3s
    fpar%nng1=tfft%nsp( parai%me+1 ) * fpar%kr3s
    DEALLOCATE(thread_buff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    IF (paral%io_parent) THEN
       WRITE(6,'(A,A)') '  NCPU     NGW',&
            '     NHG  PLANES  GXRAYS  HXRAYS ORBITALS Z-PLANES'
       DO i=0,parai%nproc-1
          iorb=parap%nst12(i,2)-parap%nst12(i,1)+1
          izpl=parap%nrzpl(i,2)-parap%nrzpl(i,1)+1
          WRITE(6,'(I6,8I8)') i,parap%sparm(3,i),parap%sparm(1,i),parap%sparm(5,i),&
               parap%sparm(9,i),parap%sparm(8,i),iorb,izpL,tfft%nr3p(i+1)
          ! IF(SPARM(3,I).LE.0) CALL stopgm(procedureN,
          ! *            'NGW .LE. 0')
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! DEALLOCATE DATA ARRAYS
    ! ==--------------------------------------------------------------==
    DEALLOCATE(mgpa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ixray,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ihray,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE loadpa
  ! ==================================================================
  SUBROUTINE iraymax(n1,n2,ir,i1,i2,thread_buff,nthreads)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n1, n2, ir(n1,n2), i1, i2, &
                                                nthreads, &
                                                thread_buff(2,0:nthreads-1)

    INTEGER                                  :: im, ithread, j1, j2

! Input
! Output
! ==--------------------------------------------------------------==

    thread_buff(:,:) = 0
    !$omp parallel &
    !$omp          private(J1,J2,IM,ithread) &
    !$omp          shared(N1,N2,thread_buff)
    im=-2**30
    ithread = 0
    !$    ithread = omp_get_thread_num ( )
    !$omp do __COLLAPSE2
    DO j2=1,n2
       DO j1=1,n1
          IF (ir(j1,j2).GE.im) THEN
             im=ir(j1,j2)
             thread_buff(1,ithread) = j1
             thread_buff(2,ithread) = j2
          ENDIF
       ENDDO
    ENDDO
    !$omp enddo
    !$omp end parallel

    ! finally find the max over threads
    i1=thread_buff(1,0)
    i2=thread_buff(2,0)
    im=ir(i1,i2)
    DO ithread=1,nthreads-1
       j1=thread_buff(1,ithread)
       j2=thread_buff(2,ithread)
       ! protect if a thread doesnt have an entry
       IF (j1>0.AND.j2>0) THEN
          IF (ir(j1,j2).GE.im) THEN
             im=ir(j1,j2)
             i1=j1
             i2=j2
          ENDIF
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE iraymax
  ! ==================================================================
  FUNCTION ixypr(nfull,kr1,kr2,kr3)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nfull, kr1, kr2, kr3, ixypr

    INTEGER                                  :: i2, i3, n2

! ==--------------------------------------------------------------==

    i3=nfull/(kr1*kr2) + 1
    n2=nfull - (i3-1)*kr1*kr2
    i2=n2/kr1 + 1
    ixypr = (i3 - 1)*kr2 + i2
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION ixypr
  ! ==================================================================
  SUBROUTINE leadim(nr1,nr2,nr3,kr1,kr2,kr3)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr1, nr2, nr3, kr1, kr2, kr3

! ==--------------------------------------------------------------==
! to align things properly for the fft we need that
! KR1=NR1+MOD(NR1,2)
! KR2=NR2+MOD(NR2,2)
! KR3=NR3+MOD(NR3,2)
! instead off that
    !TK kr[123]=nr[123] seems to be much faster on Intel architectures...
#ifdef _INTEL_MKL
    kr1=nr1
    kr2=nr2
    kr3=nr3
#else
    kr1=nr1+MOD(nr1+1,2)
    kr2=nr2+MOD(nr2+1,2)
    kr3=nr3+MOD(nr3+1,2)
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE leadim
  ! ==================================================================
  SUBROUTINE zfft(iray,gvcut)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iray(fpar%kr1,fpar%kr2)
    REAL(real_8)                             :: gvcut

    INTEGER                                  :: i, id1, id2, id3, in1, in2, &
                                                in3, ir, ix1, ix2, j, jmax, &
                                                jmin, k, kmax, kmin, nh1, &
                                                nh2, nh3
    REAL(real_8)                             :: g2, t

! ==--------------------------------------------------------------==

    CALL zeroing(iray)!,kr2*kr3)
    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    DO i=0,parm%nr3-1
       jmin=-parm%nr1+1
       jmax=parm%nr1-1
       IF (i.EQ.0) THEN
          jmin=0
       ENDIF
       DO j=jmin,jmax
          kmin=-parm%nr2+1
          kmax=parm%nr2-1
          IF (i.EQ.0.AND.j.EQ.0) THEN
             kmin=0
          ENDIF
          DO k=kmin,kmax
             g2=0._real_8
             DO ir=1,3
                t=REAL(i,kind=real_8)*gvec_com%b1(ir)+REAL(j,kind=real_8)*gvec_com%b2(ir)+REAL(k,kind=real_8)*gvec_com%b3(ir)
                g2=g2+t*t
             ENDDO
!             IF (g2.LT.gvcut) THEN
             IF (compare_lt(g2,gvcut)) THEN
                in1=nh1+i
                in2=nh2+j
                in3=nh3+k
                id1=2*nh1-in1
                id2=2*nh2-in2
                id3=2*nh3-in3
                ix1=iray(in1,in2)
                ix2=iray(id1,id2)
                IF (ix1.GT.0) THEN
                   iray(in1,in2)=ix1+1
                   IF (ix2.GT.0) THEN
                      IF ((in1.NE.id1).OR.(in2.NE.id2)) THEN
                         IF (paral%io_parent) WRITE(6,*) ' INCONSISTENT MESH?'
                         CALL stopgm('XFFT',' ',& 
                              __LINE__,__FILE__)
                      ENDIF
                   ENDIF
                ELSEIF (ix2.GT.0) THEN
                   iray(id1,id2)=ix2+1
                ELSE
                   iray(in1,in2)=1
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zfft
  ! ==================================================================
  SUBROUTINE czfft(iray,gvcut)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iray(fpar%kr2,*)
    REAL(real_8)                             :: gvcut

    INTEGER                                  :: i, icpu1, icpu2, id1, id2, &
                                                id3, in1, in2, in3, ir, j, &
                                                jmax, jmin, k, kmax, kmin, &
                                                nh1, nh2, nh3
    REAL(real_8)                             :: g2, t

! ==--------------------------------------------------------------==

    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    DO i=0,parm%nr3-1
       jmin=-parm%nr1+1
       jmax=parm%nr1-1
       IF (i.EQ.0) THEN
          jmin=0
       ENDIF
       DO j=jmin,jmax
          kmin=-parm%nr2+1
          kmax=parm%nr2-1
          IF (i.EQ.0.AND.j.EQ.0) THEN
             kmin=0
          ENDIF
          DO k=kmin,kmax
             g2=0._real_8
             DO ir=1,3
                t=REAL(i,kind=real_8)*gvec_com%b3(ir)+REAL(j,kind=real_8)*gvec_com%b1(ir)+REAL(k,kind=real_8)*gvec_com%b2(ir)
                g2=g2+t*t
             ENDDO
             IF (compare_lt(g2,gvcut)) THEN
!             IF (g2.LT.gvcut) THEN
                in1=nh1+i
                in2=nh2+j
                in3=nh3+k
                id1=2*nh1-in1
                id2=2*nh2-in2
                id3=2*nh3-in3
                icpu1=iray(in1,in2)
                icpu2=iray(id1,id2)
                IF (icpu2.EQ.0) THEN
                   iray(id1,id2)=icpu1
                ELSEIF (icpu1.EQ.0) THEN
                   iray(in1,in2)=icpu2
                ELSEIF (icpu1.NE.icpu2) THEN
                   IF (paral%io_parent)&
                        WRITE(6,*) ' INCONSISTENT XRAY FIELDS',icpu1,icpu2
                   CALL stopgm('CZFFT',' ',&
                        __LINE__,__FILE__)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE czfft
  ! ==================================================================
  SUBROUTINE gorder
    ! ==--------------------------------------------------------------==
    ! == CHECK THE G ORDERING USING MAPGP(1:NHG)                      ==
    ! == FOR A MONO-PROCESSOR, SUM_I^NHG = 1/2 * NHGS * (NHGS+1)      ==
    ! == IN PARALLEL VERSION, WE DO THE SAME THINGS:                  ==
    ! == FOR THE IP PROC, MAPGP(IG) = IG                              ==
    ! == FOR EACH JP /= IP, WE ARE LOOKING FOR THE INDEX J:           ==
    ! ==     HG_(in JP)(J) <= HG_(in IP)(IG) < HG_(in JP)(J)          ==
    ! ==     and MAPGP(IG) = MAPGP(IG) + J                            ==
    ! == IF WE SUM ALL MAPGP WE SHOULD HAVE 1/2 * NHGS* (NHGS+1)      ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'gorder'

    INTEGER                                  :: i, ierr, ip, ipp, isub, j, l, &
                                                mho, msg1, nhgmax, nho
    REAL(real_8)                             :: xm, xt
    REAL(real_8), ALLOCATABLE                :: ho(:), hx(:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    nhgmax=0
    DO ip=0,parai%nproc-1
       nhgmax=MAX(nhgmax,parap%sparm(1,ip))
    ENDDO
    ALLOCATE(ho(nhgmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hx(nhgmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ho)!,nhgmax)
    CALL zeroing(hx)!,nhgmax)
    CALL dcopy(ncpw%nhg,hg,1,hx,1)
    msg1 = nhgmax * 8
    CALL mp_sync(parai%allgrp)
    DO ip=1,parai%nproc-1
       ipp=MOD(parai%mepos+ip,parai%nproc)
       CALL my_shift(hx,ho,msg1,parai%mepos,-ip,parai%allgrp)
       nho=parap%sparm(1,ipp)
       l=1
       DO i=1,ncpw%nhg
          IF (hg(i).LT.ho(1)) THEN
             mho=0
          ELSE
             DO j=l,nho-1
                IF (compare_ge(hg(i),ho(j)).AND.compare_lt(hg(i),ho(j+1))) THEN
!                IF (hg(i).GE.ho(j).AND.hg(i).LT.ho(j+1)) THEN
                   l=j
                   mho=j
                   GOTO 100
                ENDIF
             ENDDO
             l=nho
             mho=nho
100          CONTINUE
          ENDIF
          mapgp(i)=mapgp(i)+mho
       ENDDO
    ENDDO
    ! ..check ordering
    xm=0._real_8
    DO i=1,ncpw%ngw
       xm=xm+mapgp(i)
    ENDDO
    CALL mp_sum(xm,parai%allgrp)
    xt=REAL(spar%ngws,kind=real_8)
    xt=0.5_real_8*xt*(xt+1._real_8)-xm
    IF(ABS(XT).GT.0.1_real_8.AND.paral%io_PARENT) THEN
       WRITE(6,*) 'PROGRAMING ERROR. INFORM THE PROGRAMMER'
       WRITE(6,*) procedureN,'ERROR IN G-VEC ORDERING (NGW)'
    ENDIF
    xm=0._real_8
    DO i=1,ncpw%nhg
       xm=xm+mapgp(i)
    ENDDO
    CALL mp_sum(xm,parai%allgrp)
    xt=REAL(spar%nhgs,kind=real_8)
    xt=0.5_real_8*xt*(xt+1._real_8)-xm
    IF(ABS(XT).GT.0.1_real_8.AND.paral%io_PARENT) THEN
       WRITE(6,*) 'PROGRAMING ERROR. INFORM THE PROGRAMMER'
       wRITE(6,*) procedureN,'ERROR IN G-VEC ORDERING (NHG)'
    ENDIF
    DEALLOCATE(ho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gorder
  ! ==================================================================

  SUBROUTINE gsort(hg,inyh,nhg)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: inyh(3,*), nhg
    REAL(real_8)                             :: hg(nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gsort'

    INTEGER                                  :: icurr, ierr, ig, it, j
    INTEGER, ALLOCATABLE                     :: INDEX(:)

! REORDER THE G S IN ORDER OF INCREASING MAGNITUDE

    ALLOCATE(INDEX(nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL sort2(hg,nhg,index)
    DO 25 ig=1,nhg-1
       icurr=ig
30     CONTINUE
       IF (INDEX(icurr).NE.ig) THEN
          DO j=1,3
             it=inyh(j,icurr)
             inyh(j,icurr)=inyh(j,INDEX(icurr))
             inyh(j,INDEX(icurr))=it
          END DO
          it=icurr
          icurr=INDEX(icurr)
          INDEX(it)=it
          IF (INDEX(icurr).EQ.ig) THEN
             INDEX(icurr)=icurr
             GOTO 25
          ENDIF
          GOTO 30
       ENDIF
25  CONTINUE
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gsort

  ! ******************************************************************************
  PURE FUNCTION compare_ge(a,b) RESULT(greater_equal)
    REAL(real_8), INTENT(IN)                 :: a,b
    REAL(real_8), PARAMETER                  :: eps=EPSILON(1.0_real_8)
    LOGICAL                                  :: greater_equal
    
    greater_equal=.FALSE.
    IF(ABS(a-b).LT.eps)THEN
       !a.eq.b
       greater_equal=.TRUE.
    ELSE
       IF(a.GT.b)THEN
          greater_equal=.TRUE.
       END IF
    END IF
    
  END FUNCTION compare_ge
  ! ******************************************************************************
  PURE FUNCTION compare_lt(a,b) RESULT(lower)
    REAL(real_8), INTENT(IN)                 :: a,b
    REAL(real_8), PARAMETER                  :: eps=EPSILON(1.0_real_8)
    LOGICAL                                  :: lower
    
    lower=.FALSE.
    IF(ABS(a-b).LT.eps)THEN
       !a.eq.b
    ELSE
       IF(a.LT.b)THEN
          lower=.TRUE.
       END IF
    END IF
    
  END FUNCTION compare_lt
  ! ******************************************************************************
  SUBROUTINE Distribute_Sticks( n1, n2, iray, iray2 )
    IMPLICIT NONE

    INTEGER, INTENT(IN)                         :: n1, n2
    INTEGER, INTENT(INOUT)                      :: iray( n1, n2 )
    INTEGER, OPTIONAL, INTENT(INOUT)            :: iray2( n1, n2 )

!    REAL(real_8), PARAMETER                     :: thresh
    INTEGER                                     :: m1, m2, ifa, jfa, i, j, highest, proc, c_proc, indx1, indx2, k
    LOGICAL                                     :: finished
    INTEGER                                     :: nst( parai%nproc ), ngv( parai%nproc )

    finished = .false.

    m1 = ( n1 / 2 ) + 1
    m2 = ( n2 / 2 ) + 1

    nst = 0
    ngv = 0

    IF( present( iray2 ) ) THEN

       DO i = 1, n1
          DO j = 1, n2
             IF( iray2( i, j ) .ne. 0 ) THEN
                nst( -iray2( i, j ) ) = nst( -iray2( i, j ) ) + 1
                ngv( -iray2( i, j ) ) = ngv( -iray2( i, j ) ) + iray( i, j )
                iray( i, j ) = iray2( i, j )
             END IF
          ENDDO
       ENDDO

       tfft%npst = tfft%nwst
       k = tfft%nwst
       DO ifa = 1, n2
          i = MOD( m2 + ifa - 1 - 1, n2 ) + 1
          DO jfa = 1, n1
             j = MOD( m1 + jfa - 1 - 1, n1 ) + 1
             IF( iray( j, i ) .gt. 0 ) THEN
                tfft%npst = tfft%npst + 1
                tfft%indx_map( j, i ) = tfft%npst
             END IF
          ENDDO
       ENDDO

    ELSE

       tfft%nwst = 0
       k = 0
       DO ifa = 1, n2
          i = MOD( m2 + ifa - 1 - 1, n2 ) + 1
          DO jfa = 1, n1
             j = MOD( m1 + jfa - 1 - 1, n1 ) + 1
             IF( iray( j, i ) .gt. 0 ) THEN
                tfft%nwst = tfft%nwst + 1
                tfft%indx_map( j, i ) = tfft%nwst
             END IF
          ENDDO
       ENDDO

    END IF


    DO while( .not. finished )

       highest = 0
       indx1 = 0
       indx2 = 0
       k = k + 1
       DO ifa = 1, n2
          i = MOD( m2 + ifa - 1 - 1, n2 ) + 1
          DO jfa = 1, n1
             j = MOD( m1 + jfa - 1 - 1, n1 ) + 1
             IF( iray( j, i ) .gt. highest ) THEN
                highest = iray( j, i )
                indx1 = j
                indx2 = i
             END IF
          ENDDO
       ENDDO

       IF( highest .ne. 0 ) THEN
          c_proc = 1
          DO proc = 1, parai%nproc
             IF( ngv(proc) .lt. ngv(c_proc) ) THEN
                c_proc = proc
             ELSEIF( ngv(proc) .eq. ngv(c_proc) .and. nst(proc) .lt. nst(c_proc) ) THEN
                c_proc = proc
             END IF
          ENDDO
          iray( indx1, indx2 ) = - c_proc
          ngv(c_proc) = ngv(c_proc) + highest
          nst(c_proc) = nst(c_proc) + 1
          tfft%indx( k ) = tfft%indx_map( indx1, indx2 )
          tfft%ind1( tfft%indx( k ) ) = indx1
          tfft%ind2( tfft%indx( k ) ) = indx2
       ELSE
          finished = .true.
       END IF

    END DO


  END SUBROUTINE
  ! ******************************************************************************
  SUBROUTINE hpsort_eps (n, ra, ind, eps)
    !---------------------------------------------------------------------
    !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
    !! and considering two elements being equal if their values differ
    !! for less than "eps".
    !! \(\text{n}\) is input, \(\text{ra}\) is replaced on output by its
    !! sorted rearrangement.
    !! Create an index table (ind) by making an exchange in the index array
    !! whenever an exchange is made on the sorted data array (\(\text{ra}\)).
    !! In case of equal values in the data array (\(\text{ra}\)) the values
    !! in the index array (ind) are used to order the entries.
    !! If on input ind(1) = 0 then indices are initialized in the routine,
    !! if on input ind(1) != 0 then indices are assumed to have been
    !! initialized before entering the routine and these indices are carried
    !! around during the sorting process.
    !
    ! no work space needed !
    ! free us from machine-dependent sorting-routines !
    !
    ! adapted from Numerical Recipes pg. 329 (new edition)
    !
    implicit none
    !-input/output variables
    integer, intent(in) :: n
    integer, intent(inout) :: ind (*)
    real(real_8), intent(inout) :: ra (*)
    real(real_8), intent(in) :: eps
    !-local variables
    integer :: i, ir, j, l, iind
    real(real_8) :: rra
    ! initialize index array
    if (ind (1) .eq.0) then
       do i = 1, n
          ind (i) = i
       enddo
    endif
    ! nothing to order
    if (n.lt.2) return
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1

    ir = n

    sorting: do

      ! still in hiring phase
      if ( l .gt. 1 ) then
         l    = l - 1
         rra  = ra (l)
         iind = ind (l)
         ! in retirement-promotion phase.
      else
         ! clear a space at the end of the array
         rra  = ra (ir)
         !
         iind = ind (ir)
         ! retire the top of the heap into it
         ra (ir) = ra (1)
         !
         ind (ir) = ind (1)
         ! decrease the size of the corporation
         ir = ir - 1
         ! done with the last promotion
         if ( ir .eq. 1 ) then
            ! the least competent worker at all !
            ra (1)  = rra
            !
            ind (1) = iind
            exit sorting
         endif
      endif
      ! wheter in hiring or promotion phase, we
      i = l
      ! set up to place rra in its proper level
      j = l + l
      !
      do while ( j .le. ir )
         if ( j .lt. ir ) then
            ! compare to better underling
            if ( abs(ra(j)-ra(j+1)).ge.eps ) then
               if (ra(j).lt.ra(j+1)) j = j + 1
            else
               ! this means ra(j) == ra(j+1) within tolerance
               if (ind (j) .lt.ind (j + 1) ) j = j + 1
            endif
         endif
         ! demote rra
         if ( abs(rra - ra(j)).ge.eps ) then
            if (rra.lt.ra(j)) then
               ra (i) = ra (j)
               ind (i) = ind (j)
               i = j
               j = j + j
            else
               ! set j to terminate do-while loop
               j = ir + 1
            end if
         else
            !this means rra == ra(j) within tolerance
            ! demote rra
            if (iind.lt.ind (j) ) then
               ra (i) = ra (j)
               ind (i) = ind (j)
               i = j
               j = j + j
            else
               ! set j to terminate do-while loop
               j = ir + 1
            endif
         end if
      enddo
      ra (i) = rra
      ind (i) = iind

    end do sorting
    !
  END SUBROUTINE hpsort_eps
  ! ******************************************************************************
  SUBROUTINE SetupArrays( ihray, ixray )
    IMPLICIT NONE

    INTEGER, ALLOCATABLE, INTENT(IN) :: ihray(:,:), ixray(:,:)
    CHARACTER(*), PARAMETER                     :: procedureN = 'SetupArrays'

    INTEGER                                     :: ierr, i, istick, i1, i2, ix, iy, ixM, i1M, i2M, k, f, j, ifa, jfa

    tfft%which = 1

    ALLOCATE(tfft%nr3p(parai%nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%nr3p_offset(parai%nproc),STAT=ierr)

    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    tfft%nr3p = 0
    tfft%nr3p_offset = 0
    j = 0
    DO i = 1, fpar%kr3s
       j = mod( j, parai%nproc ) + 1
       tfft%nr3p( j ) = tfft%nr3p( j ) + 1
    ENDDO
    DO i = 1, parai%nproc
       IF( i .eq. 1 ) THEN
          tfft%nr3p_offset(i) = 0
       ELSE
          tfft%nr3p_offset(i) = tfft%nr3p_offset(i-1) + tfft%nr3p(i-1)
       END IF
    ENDDO
    tfft%my_nr3p = tfft%nr3p( parai%me+1 )
    tfft%nr3px   = MAXVAL ( tfft%nr3p )

    ALLOCATE(tfft%ir1w(fpar%kr1s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%ir1p(fpar%kr1s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%indw(fpar%kr1s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%indp(fpar%kr1s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%nsw(parai%nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tfft%nsp(parai%nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    tfft%nsw = 0
    tfft%nsp = 0
    tfft%ir1w = 0
    tfft%ir1p = 0
    tfft%nr1w = 0
    tfft%nr1p = 0
    DO i = 1, tfft%npst
       istick = tfft%indx( i )
       i1 = tfft%ind1( istick )
       i2 = tfft%ind2( istick )

       ix = i1 - ( fpar%kr1s / 2 )

       IF( ihray( i1, i2 ) .ne. 0 ) THEN

          IF( ixray( i1, i2 ) .ne. 0 ) THEN

             IF( tfft%ir1w( ix ) .eq. 0 ) THEN
                tfft%ir1w(ix) = -1
                tfft%ir1p(ix) = -1
             END IF
             tfft%nsw( - ixray( i1, i2 ) ) = tfft%nsw( - ixray( i1, i2 ) ) + 1
             tfft%nsp( - ixray( i1, i2 ) ) = tfft%nsp( - ixray( i1, i2 ) ) + 1

          ELSE

             IF( tfft%ir1p( ix ) .eq. 0 ) THEN
                tfft%ir1p(ix) = -1
             END IF
             tfft%nsp( - ihray( i1, i2 ) ) = tfft%nsp( - ihray( i1, i2 ) ) + 1

          END IF


          i1M = 2 * ( fpar%kr1s / 2 ) + 2 - i1
          i2M = 2 * ( fpar%kr2s / 2 ) + 2 - i2
          IF( i1M .eq. i1 .and. i2M .eq. i2 ) CYCLE
          ixM = i1M - ( fpar%kr1s / 2 ) + fpar%kr1s
          IF( ixM .gt. fpar%kr1s ) ixM = ixM - fpar%kr1s
          IF( ixray( i1M , i2M ) .ne. 0 ) THEN

             IF( tfft%ir1w( ixM ) .eq. 0 ) THEN
                tfft%ir1w(ixM) = -1
                tfft%ir1p(ixM) = -1
             END IF
             tfft%nsw( - ixray( i1M, i2M ) ) = tfft%nsw( - ixray( i1M, i2M ) ) + 1
             tfft%nsp( - ixray( i1M, i2M ) ) = tfft%nsp( - ixray( i1M, i2M ) ) + 1

          ELSE

             IF( tfft%ir1p( ixM ) .eq. 0 ) THEN
                tfft%ir1p(ixM) = -1
             END IF
             tfft%nsp( - ihray( i1M, i2M ) ) = tfft%nsp( - ihray( i1M, i2M ) ) + 1

          END IF


       END IF

    ENDDO

    DO i = 1, fpar%kr1s

       IF( tfft%ir1w( i ) .lt. 0 ) THEN
          tfft%nr1w = tfft%nr1w + 1
          tfft%ir1w( i ) = tfft%nr1w
          tfft%ir1p( i ) = tfft%nr1w
          tfft%indw( tfft%nr1w ) = i
          tfft%indp( tfft%nr1w ) = i
       END IF

    ENDDO
    tfft%nr1p = tfft%nr1w
    DO i = 1, fpar%kr1s

       IF( tfft%ir1p( i ) .lt. 0 ) THEN
          tfft%nr1p = tfft%nr1p + 1
          tfft%ir1p( i ) = tfft%nr1p
          tfft%indp( tfft%nr1p ) = i
       END IF

    ENDDO


    ALLOCATE(tfft%iss(parai%nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    DO i = 1, parai%nproc
      IF( i .eq. 1 ) THEN
        tfft%iss( i ) = 0
      ELSE
        tfft%iss( i ) = tfft%iss( i - 1 ) + tfft%nsp( i - 1 )
      ENDIF
    ENDDO


    ALLOCATE(tfft%ismap(fpar%kr1s*fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    tfft%ismap = 0
    tfft%nsp = 0
    k = 0
    f = 0
    DO ifa = 1, fpar%kr2s
       i = MOD( ( fpar%kr2s / 2 ) + 1 + ifa - 1 - 1, fpar%kr2s ) + 1
       DO jfa = 1, fpar%kr1s
          j = MOD( ( fpar%kr1s / 2 ) + 1 + jfa - 1 - 1, fpar%kr1s ) + 1
          f = f + 1
          IF( ixray( j, i ) .lt. 0 ) THEN
             tfft%nsp( - ixray( j, i ) ) = tfft%nsp( - ixray( j, i ) ) + 1
             tfft%ismap( tfft%nsp( - ixray( j, i ) ) + tfft%iss( - ixray( j, i ) ) ) = f
          END IF
       ENDDO
    ENDDO

    k = 0
    f = 0
    DO ifa = 1, fpar%kr2s
       i = MOD( ( fpar%kr2s / 2 ) + 1 + ifa - 1 - 1, fpar%kr2s ) + 1
       DO jfa = 1, fpar%kr1s
          j = MOD( ( fpar%kr1s / 2 ) + 1 + jfa - 1 - 1, fpar%kr1s ) + 1
          f = f + 1
          IF( ixray( j, i ) .eq. 0 .and. ihray( j, i ) .lt. 0 ) THEN
             tfft%nsp( - ihray( j, i ) ) = tfft%nsp( - ihray( j, i ) ) + 1
             tfft%ismap( tfft%nsp( - ihray( j, i ) ) + tfft%iss( - ihray( j, i ) ) ) = f
          END IF
       ENDDO
    ENDDO

    tfft%small_chunks(1) = tfft%nr3px * MAXVAL( tfft%nsw )
    tfft%small_chunks(2) = tfft%nr3px * MAXVAL( tfft%nsp )
    tfft%big_chunks(1)   = tfft%small_chunks(1) * parai%max_node_nproc * parai%max_node_nproc
    tfft%big_chunks(2)   = tfft%small_chunks(2) * parai%max_node_nproc * parai%max_node_nproc
    tfft%tscale = 1.0d0 / dble( fpar%kr1s * fpar%kr2s * fpar%kr3s )

  END SUBROUTINE SetupArrays

END MODULE loadpa_utils
