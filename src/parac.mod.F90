MODULE parac
#ifdef __PARALLEL
    USE mpi_f08
#endif

  IMPLICIT NONE

  ! ==--------------------------------------------------------------==
  ! == PARENT  : .TRUE if IP = SOURCE                               ==
  ! == IO_PARENT: .TRUE if IP = IO_SOURCE (cp_grp communicator)     ==
  ! ==--------------------------------------------------------------==
  TYPE :: paral_t
     LOGICAL :: parent = .FALSE.
     LOGICAL :: qmmmparent = .FALSE.
     LOGICAL :: qmnode = .FALSE.
     LOGICAL :: io_parent = .FALSE.
     LOGICAL :: cp_inter_io_parent = .FALSE.
  END TYPE paral_t
  TYPE(paral_t), SAVE :: paral
  ! ==--------------------------------------------------------------==
  ! == NCPUS   : value of NCPUS environment variable                ==
  ! == NPROC   : TOTAL NUMBER OF PROCESSORS                         ==
  ! == ME      : CPU index                                          ==
  ! == MEPOS   : equal to ME                                        ==
  ! == SOURCE  : MASTER CPU index                                   ==
  ! == IO_SOURCE: MASTER CPU index for IO (cp_grp communicator)     ==
  ! == IGEQ0   : PROCESSOR INDEX OF G=0 COMPONENT                   ==
  ! == ALLGRP  : GIVEN VALUE FOR BROADCAST COMMUNICATION            ==
  ! == CP_GRP  : CPMD communicator (should be use instead of        ==
  ! MPI_COMM_WORLD)                                    ==
  ! == CP_NPROC: Nbr processes in the cpmd communicator             ==
  ! == CP_ME   : id of the processes in the cpmd communicator       ==
  ! == NHRAYS  : number of rays for the processor for the density   ==
  ! == NGRAYS  : number of rays for the processor for the wavefunc. ==
  ! == CP_INTER_GRP: CPMD group communicator                        ==
  ! == CP_INTER_ME: CPMD group index                                ==
  ! == CP_NOGRP: number of CPMD group                               ==
  ! == CP_INTER_IO_SOURCE: MASTER intergroup index for IO           ==
  ! == LOC_GRP
  ! == LOC_ME
  ! == LOC_NPROC
  ! == LOC_INTER_GRP
  ! == node_grp : subgroup of allgrp
  ! == node_nproc : nbr of procs
  ! == node_me : index
  ! == cp_inter_node_grp : sub group of cp_inter_grp
  ! == cp_inter_node_nproc : nbr of procs
  ! == cp_inter_node_me : index
  ! ==================================================================
  ! == NEEDED FOR NEW GDISTRIBUTION FFT
  ! ==================================================================
  ! == send_handle : send handle for shared memory communication
  ! == recv_handle : recv handle for shared memory communication
  ! == c2_send_handle : send handle for c2 communication
  ! == c2_recv_handle : recv handle for c2 communication
  ! ==--------------------------------------------------------------==
  TYPE :: parai_t
     INTEGER :: ncpus = HUGE(0)
     INTEGER :: ncpus_FFT = HUGE(0)
     INTEGER :: nproc = HUGE(0)
     INTEGER :: me = HUGE(0)
     INTEGER :: mepos = HUGE(0)
     INTEGER :: source = HUGE(0)
     INTEGER :: igeq0 = HUGE(0)
#ifdef __PARALLEL
     type(MPI_COMM) :: allgrp
#else
     INTEGER :: allgrp = HUGE(0)
#endif
     INTEGER :: nhrays = HUGE(0)
     INTEGER :: ngrays = HUGE(0)
     INTEGER :: qmmmnproc = HUGE(0)
#ifdef __PARALLEL
     type(MPI_COMM) :: qmmmgrp
#else
     INTEGER :: qmmmgrp = HUGE(0)
#endif

     INTEGER :: qmmmsource = HUGE(0)
#ifdef __PARALLEL
     type(MPI_COMM) :: cp_grp
#else
     INTEGER :: cp_grp = HUGE(0)
#endif
     INTEGER :: cp_nproc = HUGE(0)
     INTEGER :: cp_me = HUGE(0)
#ifdef __PARALLEL
     type(MPI_COMM) :: cp_inter_grp
#else
     INTEGER :: cp_inter_grp = HUGE(0)
#endif
     INTEGER :: cp_inter_me = HUGE(0)
     INTEGER :: cp_nogrp = HUGE(0)
     INTEGER :: io_source = HUGE(0)
     INTEGER :: cp_inter_io_source = HUGE(0)
#ifdef __PARALLEL
     type(MPI_COMM) :: loc_grp
#else
     INTEGER :: loc_grp = HUGE(0)
#endif
     INTEGER :: loc_me = HUGE(0)
     INTEGER :: loc_nproc = HUGE(0)
#ifdef __PARALLEL
     type(MPI_COMM) :: loc_inter_grp
#else
     INTEGER :: loc_inter_grp = HUGE(0)
#endif
#ifdef __PARALLEL
     type(MPI_COMM) :: node_grp
#else
     INTEGER :: node_grp = HUGE(0)
#endif
     INTEGER :: node_nproc = HUGE(0)
     INTEGER :: node_me = HUGE(0)
#ifdef __PARALLEL
     type(MPI_COMM) :: cp_inter_node_grp
#else
     INTEGER :: cp_inter_node_grp = HUGE(0)
#endif
     INTEGER :: cp_inter_node_nproc = HUGE(0)
     INTEGER :: cp_inter_node_me = HUGE(0)
#ifdef __PARALLEL
     type(MPI_REQUEST), ALLOCATABLE :: sendrecv_handle(:,:,:,:)
     type(MPI_REQUEST), ALLOCATABLE :: c2_send_handle(:,:)
     type(MPI_REQUEST), ALLOCATABLE :: c2_recv_handle(:,:)
     type(MPI_REQUEST), ALLOCATABLE :: c2_comb_handle(:,:)
     type(MPI_REQUEST), ALLOCATABLE :: c2_send_handle_counter(:)
     type(MPI_REQUEST), ALLOCATABLE :: c2_recv_handle_counter(:)
     type(MPI_REQUEST), ALLOCATABLE :: c2_comb_handle_counter(:)
     type(MPI_REQUEST), ALLOCATABLE :: sendrecv_handle2(:,:,:,:)
#endif
     INTEGER :: nnode = HUGE(0)
     INTEGER :: my_node = HUGE(0)
     INTEGER, ALLOCATABLE :: cp_overview(:,:)
     INTEGER :: max_node_nproc = HUGE(0)
     INTEGER, ALLOCATABLE :: node_nproc_overview(:)
     INTEGER, ALLOCATABLE :: node_grpindx(:)
  END TYPE parai_t
  TYPE(parai_t), SAVE :: parai

END MODULE parac
