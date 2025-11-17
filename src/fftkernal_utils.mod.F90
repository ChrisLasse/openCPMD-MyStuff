#include "cpmd_global.h"


MODULE fftkernal_utils

  USE fft,                             ONLY: FFT_TYPE_DESCRIPTOR
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntl,&
                                             fpar

  IMPLICIT NONE

  PRIVATE

  ABSTRACT INTERFACE
    SUBROUTINE z2y_fillcom_t( tfft, aux, comm_mem, batch_size, remswitch, mythread, nss )
      IMPORT FFT_TYPE_DESCRIPTOR, real_8, fpar
      TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
      INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
      INTEGER, INTENT(IN) :: nss( * )
      COMPLEX(real_8), INTENT(IN)  :: aux ( fpar%kr3s , * )
      COMPLEX(real_8), INTENT(OUT) :: comm_mem( * )
    END SUBROUTINE z2y_fillcom_t
  END INTERFACE
  ABSTRACT INTERFACE
    SUBROUTINE y2z_fillcom_t( tfft, aux, comm_mem, batch_size, mythread, ispec, map_y2z )
      IMPORT FFT_TYPE_DESCRIPTOR, real_8
      TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
      INTEGER, INTENT(IN) :: batch_size, mythread, ispec
      INTEGER, INTENT(IN) :: map_y2z( * )
      COMPLEX(real_8), INTENT(IN)  :: aux ( * )
      COMPLEX(real_8), INTENT(OUT) :: comm_mem( * )
    END SUBROUTINE y2z_fillcom_t
  END INTERFACE
  ABSTRACT INTERFACE
    SUBROUTINE y2z_voidrecv_t( tfft, comm_mem_recv, aux, batch_size, remswitch, mythread, factor, nss )
      IMPORT FFT_TYPE_DESCRIPTOR, real_8, fpar
      TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
      INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
      INTEGER, INTENT(IN) :: nss( * )
      DOUBLE PRECISION, INTENT(IN) :: factor
      COMPLEX(real_8), INTENT(IN)  :: comm_mem_recv( * )
      COMPLEX(real_8), INTENT(OUT) :: aux ( fpar%kr3s , * )
    END SUBROUTINE y2z_voidrecv_t
  END INTERFACE

  PROCEDURE(z2y_fillcom_t),  POINTER :: z2y_fillsend
  PROCEDURE(z2y_fillcom_t),  POINTER :: z2y_fillrecv
  PROCEDURE(y2z_fillcom_t),  POINTER :: y2z_fillsend
  PROCEDURE(y2z_fillcom_t),  POINTER :: y2z_fillrecv
  PROCEDURE(y2z_voidrecv_t), POINTER :: y2z_voidrecv

  PUBLIC :: select_kernals
  PUBLIC :: z2y_fillsend
  PUBLIC :: z2y_fillrecv
  PUBLIC :: y2z_fillsend
  PUBLIC :: y2z_fillrecv
  PUBLIC :: y2z_voidrecv

CONTAINS

  SUBROUTINE select_kernals()
    IMPLICIT NONE

    IF( cntl%fft_distmem ) THEN
       z2y_fillsend => z2y_fillsend_distmem
       z2y_fillrecv => z2y_fillrecv_distmem
       y2z_fillsend => y2z_fillsend_distmem
       y2z_fillrecv => y2z_fillrecv_distmem
       y2z_voidrecv => y2z_voidrecv_distmem
    ELSE
       z2y_fillsend => z2y_fillsend_sharedmem
       z2y_fillrecv => z2y_fillrecv_sharedmem
       y2z_fillsend => y2z_fillsend_sharedmem
       y2z_fillrecv => y2z_fillrecv_sharedmem
       y2z_voidrecv => y2z_voidrecv_sharedmem
    END IF

  END SUBROUTINE select_kernals

  SUBROUTINE z2y_fillsend_sharedmem( tfft, aux, comm_mem_send, batch_size, remswitch, mythread, nss )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
    INTEGER, INTENT(IN) :: nss( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( fpar%kr3s , * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_send( * )

    INTEGER :: i, l, j, m, offset, k, kdest

    IF( parai%nnode .ne. 1 ) THEN

       j = 0
       DO l = 1, parai%nnode
          IF( l .eq. parai%my_node+1 ) THEN
             j = j + parai%node_nproc_overview( l )
             CYCLE
          END IF
          DO m = 1, parai%node_nproc_overview( l )
             j = j + 1
             !     ( Where am I on the node + to which proc does it go ) * Package size
             offset = ( parai%node_me + ((l-1)*parai%max_node_nproc+(m-1))*parai%max_node_nproc ) * tfft%small_chunks(tfft%which) + (l-1)*(batch_size-1) * tfft%big_chunks(tfft%which)
             DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
                kdest = offset + tfft%nr3px * mod( (k-1), nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%big_chunks(tfft%which)
                DO i = 1, tfft%nr3p( j )
                   comm_mem_send( kdest + i ) = aux( i + tfft%nr3p_offset( j ), k )
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    END IF

  END SUBROUTINE z2y_fillsend_sharedmem

  SUBROUTINE z2y_fillsend_distmem( tfft, aux, comm_mem_send, batch_size, remswitch, mythread, nss )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
    INTEGER, INTENT(IN) :: nss( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( fpar%kr3s , * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_send( * )

    INTEGER :: i, l, j, m, offset, k, kdest

    IF( parai%cp_nproc .ne. 1 ) THEN

       j = 0
       DO l = 1, parai%nnode
          DO m = 1, parai%node_nproc_overview( l )
             j = j + 1
!             IF( parai%cp_me+1 .eq. j ) CYCLE
             !     ( Where am I on the node + to which proc does it go ) * Package size
             offset = (j-1) * tfft%small_chunks(tfft%which) * batch_size
             DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
                kdest = offset + tfft%nr3px * mod( (k-1), nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%small_chunks(tfft%which)
                DO i = 1, tfft%nr3p( j )
                   comm_mem_send( kdest + i ) = aux( i + tfft%nr3p_offset( j ), k )
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    END IF

  END SUBROUTINE z2y_fillsend_distmem

  SUBROUTINE z2y_fillrecv_sharedmem( tfft, aux, comm_mem_recv, batch_size, remswitch, mythread, nss )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
    INTEGER, INTENT(IN) :: nss( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( fpar%kr3s , * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_recv( * )

    INTEGER :: i, m, offset, k, kdest

    DO m = 1, parai%node_nproc
       offset = ( parai%node_me + (parai%my_node*parai%max_node_nproc+(m-1))*parai%max_node_nproc) * tfft%small_chunks(tfft%which) + parai%my_node*(batch_size-1) * tfft%big_chunks(tfft%which)
       DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
          kdest = offset + tfft%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%big_chunks(tfft%which)
          DO i = 1, tfft%nr3p( parai%node_grpindx(m)+1 )
             comm_mem_recv( kdest + i ) = aux( i + tfft%nr3p_offset( parai%node_grpindx(m)+1 ), k )
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE z2y_fillrecv_sharedmem

  SUBROUTINE z2y_fillrecv_distmem( tfft, aux, comm_mem_recv, batch_size, remswitch, mythread, nss )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
    INTEGER, INTENT(IN) :: nss( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( fpar%kr3s , * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_recv( * )

!    INTEGER :: i, k, kdest
    INTEGER :: j, l, m, offset, i, kdest, k

!    DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
!       kdest = parai%cp_me * tfft%small_chunks(tfft%which) * batch_size + tfft%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%small_chunks(tfft%which)
!       DO i = 1, tfft%nr3p( parai%cp_me+1 )
!          comm_mem_recv( kdest + i ) = aux( i + tfft%nr3p_offset( parai%cp_me+1 ), k )
!       ENDDO
!    ENDDO

    j = 0
    DO l = 1, parai%nnode
       DO m = 1, parai%node_nproc_overview( l )
          j = j + 1
          !     ( Where am I on the node + to which proc does it go ) * Package size
          offset = (j-1) * tfft%small_chunks(tfft%which) * batch_size
          DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
             kdest = offset + tfft%nr3px * mod( (k-1), nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%small_chunks(tfft%which)
             DO i = 1, tfft%nr3p( j )
                comm_mem_recv( kdest + i ) = aux( i + tfft%nr3p_offset( j ), k )
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE z2y_fillrecv_distmem

  SUBROUTINE y2z_fillsend_sharedmem( tfft, aux, comm_mem_send, batch_size, mythread, ispec, map_y2z )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, ispec
    INTEGER, INTENT(IN) :: map_y2z( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_send( * )

    INTEGER :: i, l, j, m, offset, k, kdest, offset2

    offset2 =  (ispec-1) * tfft%big_chunks(tfft%which)

    IF( parai%nnode .ne. 1 ) THEN

       i = 0
       DO l = 1, parai%nnode
          IF( l .eq. parai%my_node+1 ) THEN
             i = i + parai%node_nproc_overview( l )
             CYCLE
          END IF
          DO m = 1, parai%node_nproc_overview( l )
             i = i + 1
             offset = ( parai%node_me + ((l-1)*parai%max_node_nproc+(m-1))*parai%max_node_nproc ) * tfft%small_chunks(tfft%which) + (l-1)*(batch_size-1) * tfft%big_chunks(tfft%which)
             DO j = tfft%thread_z_start( mythread+1, 3, i, tfft%which ), tfft%thread_z_end( mythread+1, 3, i, tfft%which )
                DO k = 1, tfft%my_nr3p
                   comm_mem_send( offset + offset2 + (j-1)*tfft%nr3px + k ) = &
                   aux( map_y2z( (i-1)*tfft%small_chunks(tfft%which) + (j-1)*tfft%nr3px + k ) )
                END DO
             END DO
          END DO
       END DO

    END IF

  END SUBROUTINE y2z_fillsend_sharedmem

  SUBROUTINE y2z_fillsend_distmem( tfft, aux, comm_mem_send, batch_size, mythread, ispec, map_y2z )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, ispec
    INTEGER, INTENT(IN) :: map_y2z( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_send( * )

    INTEGER :: i, l, j, m, offset, k, kdest, offset2

    offset2 =  (ispec-1) * tfft%small_chunks(tfft%which)

    IF( parai%cp_nproc .ne. 1 ) THEN

       i = 0
       DO l = 1, parai%nnode
          DO m = 1, parai%node_nproc_overview( l )
             i = i + 1
             IF( parai%cp_me+1 .eq. i ) CYCLE
             offset = (i-1) * tfft%small_chunks(tfft%which) * batch_size
             DO j = tfft%thread_z_start( mythread+1, 3, i, tfft%which ), tfft%thread_z_end( mythread+1, 3, i, tfft%which )
                DO k = 1, tfft%my_nr3p
                   comm_mem_send( offset + offset2 + (j-1)*tfft%nr3px + k ) = &
                   aux( map_y2z( (i-1)*tfft%small_chunks(tfft%which) + (j-1)*tfft%nr3px + k ) )
                END DO
             END DO
          END DO
       END DO

    END IF

  END SUBROUTINE y2z_fillsend_distmem

  SUBROUTINE y2z_fillrecv_sharedmem( tfft, aux, comm_mem_recv, batch_size, mythread, ispec, map_y2z )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, ispec
    INTEGER, INTENT(IN) :: map_y2z( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_recv( * )

    INTEGER :: j, m, offset, k, offset2

    offset2 =  (ispec-1) * tfft%small_chunks(tfft%which)

    DO m = 1, parai%node_nproc
       offset = ( parai%node_me + (parai%my_node*parai%max_node_nproc+(m-1))*parai%max_node_nproc ) * tfft%small_chunks(tfft%which) + parai%my_node*(batch_size-1) * tfft%big_chunks(tfft%which)
       DO j = tfft%thread_z_start( mythread+1, 3, parai%node_grpindx(m)+1, tfft%which ), tfft%thread_z_end( mythread+1, 3, parai%node_grpindx(m)+1, tfft%which )
          DO k = 1, tfft%my_nr3p
             comm_mem_recv( offset + offset2 + (j-1)*tfft%nr3px + k ) = &
             aux( map_y2z( parai%node_grpindx(m)*tfft%small_chunks(tfft%which) + (j-1)*tfft%nr3px + k ) )
          END DO
       END DO
    END DO

  END SUBROUTINE y2z_fillrecv_sharedmem

  SUBROUTINE y2z_fillrecv_distmem( tfft, aux, comm_mem_recv, batch_size, mythread, ispec, map_y2z )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, ispec
    INTEGER, INTENT(IN) :: map_y2z( * )
    COMPLEX(real_8), INTENT(IN)  :: aux ( * )
    COMPLEX(real_8), INTENT(OUT) :: comm_mem_recv( * )

    INTEGER :: j, k, offset2

    offset2 =  (ispec-1) * tfft%small_chunks(tfft%which)

    DO j = tfft%thread_z_start( mythread+1, 3, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, 3, parai%me+1, tfft%which )
       DO k = 1, tfft%my_nr3p
          comm_mem_recv( parai%cp_me * tfft%small_chunks(tfft%which) * batch_size + offset2 + (j-1)*tfft%nr3px + k ) = &
          aux( map_y2z( parai%me*tfft%small_chunks(tfft%which) + (j-1)*tfft%nr3px + k ) )
       END DO
    END DO

  END SUBROUTINE y2z_fillrecv_distmem

  SUBROUTINE y2z_voidrecv_sharedmem( tfft, comm_mem_recv, aux, batch_size, remswitch, mythread, factor, nss )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
    INTEGER, INTENT(IN) :: nss( * )
    DOUBLE PRECISION, INTENT(IN) :: factor
    COMPLEX(real_8), INTENT(IN)  :: comm_mem_recv( * )
    COMPLEX(real_8), INTENT(OUT) :: aux ( fpar%kr3s , * )

    INTEGER :: j, m, l, i, offset, k, kfrom

    m = 0
    DO j = 1, parai%nnode
       DO l = 1, parai%node_nproc_overview( j )
          m = m + 1
          offset = ( parai%node_me*parai%max_node_nproc + (l-1) ) * tfft%small_chunks(tfft%which) + (j-1)*batch_size * tfft%big_chunks(tfft%which)
          DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
             kfrom = offset + tfft%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%big_chunks(tfft%which)
             DO i = 1, tfft%nr3p( m )
                aux( tfft%nr3p_offset( m ) + i, k ) = comm_mem_recv( kfrom + i ) * factor
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE y2z_voidrecv_sharedmem

  SUBROUTINE y2z_voidrecv_distmem( tfft, comm_mem_recv, aux, batch_size, remswitch, mythread, factor, nss )
    IMPLICIT NONE

    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(IN) :: tfft
    INTEGER, INTENT(IN) :: batch_size, mythread, remswitch
    INTEGER, INTENT(IN) :: nss( * )
    DOUBLE PRECISION, INTENT(IN) :: factor
    COMPLEX(real_8), INTENT(IN)  :: comm_mem_recv( * )
    COMPLEX(real_8), INTENT(OUT) :: aux ( fpar%kr3s , * )

    INTEGER :: j, m, l, i, offset, k, kfrom

    m = 0
    DO j = 1, parai%nnode
       DO l = 1, parai%node_nproc_overview( j )
          m = m + 1
          offset = (m-1) * tfft%small_chunks(tfft%which) * batch_size
          DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
             kfrom = offset + tfft%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%small_chunks(tfft%which)
             DO i = 1, tfft%nr3p( m )
                aux( tfft%nr3p_offset( m ) + i, k ) = comm_mem_recv( kfrom + i ) * factor
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE y2z_voidrecv_distmem

END MODULE fftkernal_utils
