MODULE density_utils
  USE kinds,                           ONLY: real_8
  USE fft,                             ONLY: FFT_TYPE_DESCRIPTOR
  USE system,                          ONLY: fpar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: build_density_real
  PUBLIC :: build_density_imag
  PUBLIC :: build_density_sum
  PUBLIC :: build_density_sum_batch
  PUBLIC :: build_density_sum_manualOMP

CONTAINS

  ! ==================================================================
  SUBROUTINE build_density_real(alpha,psi,rho,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: alpha
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rho(:)
    INTEGER                                  :: n

    INTEGER                                  :: l

    !$omp parallel do private(L)
    DO l=1,n
       rho(l)=rho(l)+alpha*REAL(psi(l))**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_real
  ! ==================================================================
  SUBROUTINE build_density_imag(alpha,psi,rho,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: alpha
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rho(:)
    INTEGER                                  :: n

    INTEGER                                  :: l

    !$omp parallel do private(L)
    DO l=1,n
       rho(l)=rho(l)+alpha*AIMAG(psi(l))**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_imag
  ! ==================================================================
  SUBROUTINE build_density_sum(alpha_real,alpha_imag,psi,rho,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: alpha_real, alpha_imag
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rho(:)
    INTEGER                                  :: n

    INTEGER                                  :: l

    !$omp parallel do private(L)
    DO l=1,n
       rho(l)=rho(l)+alpha_real*REAL(psi(l))**2&
            +alpha_imag*AIMAG(psi(l))**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_sum
  ! ==================================================================
  SUBROUTINE build_density_sum_manualOMP( tfft, alpha_real, alpha_imag, psi, rho, mythread, spins, nspin)
    ! ==--------------------------------------------------------------==
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN)    :: tfft
    COMPLEX(real_8), INTENT(IN)              :: psi( tfft%my_nr3p * fpar%kr2s * fpar%kr1s, * )
    REAL(real_8), INTENT(INOUT)              :: rho( tfft%my_nr3p * fpar%kr2s * fpar%kr1s, * )
    INTEGER, INTENT(IN)                      :: mythread, nspin
    REAL(real_8), INTENT(IN)                 :: alpha_real, alpha_imag
    INTEGER, INTENT(IN)                      :: spins( 2 )

    INTEGER                                  :: i, j

    IF( nspin .eq. 1 ) THEN

          DO j = tfft%thread_rspace_start( mythread+1 ), tfft%thread_rspace_end( mythread+1 )
             rho( j, 1 ) = rho( j, 1 ) + alpha_real *  REAL( psi( j, 1 ) )**2 &
                                       + alpha_imag * AIMAG( psi( j, 1 ) )**2
          ENDDO

    ELSE

          DO j = tfft%thread_rspace_start( mythread+1 ), tfft%thread_rspace_end( mythread+1 )
             rho( j, spins( 1 ) ) = rho( j, spins( 1 ) ) + alpha_real *  REAL( psi( j, 1 ) )**2
             rho( j, spins( 2 ) ) = rho( j, spins( 2 ) ) + alpha_imag * AIMAG( psi( j, 1 ) )**2
          ENDDO

    END IF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_sum_manualOMP
  ! ==================================================================
  SUBROUTINE build_density_sum_batch(alpha_real,alpha_imag,psi,rho,n1,n2,n3,spins,nspin)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n1, n2, n3, nspin, spins(2,n2)
    REAL(real_8), INTENT(IN)                 :: alpha_real(n2), alpha_imag(n2)
    COMPLEX(real_8), INTENT(IN)              :: psi(n1,n2,n3)
    REAL(real_8), INTENT(INOUT)              :: rho(n1,n3,nspin)

    INTEGER                                  :: l1,l2,l3

    IF(nspin.EQ.1)THEN
       !$omp parallel do private(l1,l2,l3)
       DO l3=1,n3
          DO l2=1,n2
             DO l1=1,n1
                rho(l1,l3,1)=rho(l1,l3,1)&
                     +alpha_real(l2)*REAL(psi(l1,l2,l3),KIND=real_8)**2
                rho(l1,l3,1)=rho(l1,l3,1)&
                     +alpha_imag(l2)*AIMAG(psi(l1,l2,l3))**2
             END DO
          END DO
       END DO
    ELSE
       !$omp parallel do private(l1,l2,l3)
       DO l3=1,n3
          DO l2=1,n2
             DO l1=1,n1
                rho(l1,l3,spins(1,l2))=rho(l1,l3,spins(1,l2))&
                     +alpha_real(l2)*REAL(psi(l1,l2,l3),KIND=real_8)**2
                rho(l1,l3,spins(2,l2))=rho(l1,l3,spins(2,l2))&
                     +alpha_imag(l2)*AIMAG(psi(l1,l2,l3))**2
             END DO
          END DO
       END DO
    END IF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_sum_batch

END MODULE density_utils
