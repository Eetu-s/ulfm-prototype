! $Id$
module Ft_mpi
!
  implicit none
!
  include 'mpif.h'
  include 'mpif-ext.h'
!
  private
!
  integer, public :: pc_comm_world  = MPI_COMM_NULL
  integer, public :: pc_comm_active = MPI_COMM_NULL
  integer, public :: pc_comm_spare  = MPI_COMM_NULL
  integer, public :: pc_world_rank  = -1
  integer, public :: pc_world_size  = 0
  integer, public :: pc_active_rank = -1
  integer, public :: pc_active_size = 0
  integer, public :: pc_n_spares    = 0
  logical, public :: pc_is_spare          = .false.
  logical, public :: pc_was_promoted      = .false.
  logical, public :: pc_need_restart_step = .false.
!
  public :: ft_init, ft_check, ft_repair, ft_log
  public :: ft_set_neighbour_comm_setup, ft_set_topology_verify
!
  abstract interface
    subroutine neighbour_comm_setup_iface(parent_comm)
      integer, intent(in) :: parent_comm
    end subroutine neighbour_comm_setup_iface
    subroutine topology_verify_iface()
    end subroutine topology_verify_iface
  end interface
  procedure(neighbour_comm_setup_iface), pointer :: pc_neighbour_setup => null()
  procedure(topology_verify_iface),      pointer :: pc_topology_verify => null()
!
contains
!***********************************************************************
    subroutine ft_set_neighbour_comm_setup(p)
!
      procedure(neighbour_comm_setup_iface) :: p
!
      pc_neighbour_setup => p
!
    endsubroutine ft_set_neighbour_comm_setup
!***********************************************************************
    subroutine ft_set_topology_verify(p)
!
      procedure(topology_verify_iface) :: p
!
      pc_topology_verify => p
!
    endsubroutine ft_set_topology_verify
!***********************************************************************
    subroutine ft_init(n_spares_default)
!
      integer, intent(in) :: n_spares_default
      integer :: ierr, color, n_env
      character(len=32) :: envstr
!
      pc_comm_world = MPI_COMM_WORLD
!
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
!
      pc_n_spares = n_spares_default
      call get_environment_variable('PENCIL_FT_SPARES', envstr)
      if (len_trim(envstr) > 0) then
        read(envstr,*,iostat=ierr) n_env
        if (ierr == 0) pc_n_spares = n_env
      endif
!
      call MPI_Comm_size(MPI_COMM_WORLD, pc_world_size, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, pc_world_rank, ierr)
!
      if (pc_n_spares < 0)              pc_n_spares = 0
      if (pc_n_spares >= pc_world_size) pc_n_spares = 0
!
      if (pc_world_rank < pc_world_size - pc_n_spares) then
        color       = 0
        pc_is_spare = .false.
      else
        color       = 1
        pc_is_spare = .true.
      endif
!
      if (.not. pc_is_spare) then
        call MPI_Comm_split(MPI_COMM_WORLD, color, pc_world_rank, pc_comm_active, ierr)
        call MPI_Comm_set_errhandler(pc_comm_active, MPI_ERRORS_RETURN, ierr)
        call MPI_Comm_size(pc_comm_active, pc_active_size, ierr)
        call MPI_Comm_rank(pc_comm_active, pc_active_rank, ierr)
        if (pc_world_rank == 0) then
          write(*,'(A,I0,A,I0)') ' [ft_init] active=', pc_active_size, &
                                  ' spare=', pc_n_spares
          flush(6)
        endif
      else
        call MPI_Comm_split(MPI_COMM_WORLD, color, pc_world_rank, pc_comm_spare, ierr)
        call MPI_Comm_set_errhandler(pc_comm_spare, MPI_ERRORS_RETURN, ierr)
        call spare_wait_loop()
      endif
!
    endsubroutine ft_init
!***********************************************************************
    subroutine spare_wait_loop()
!
      integer :: ierr, shrunk_world
      integer, allocatable :: dead_ids(:)
!
      call MPIX_Comm_shrink(MPI_COMM_WORLD, shrunk_world, ierr)
      call MPI_Comm_set_errhandler(shrunk_world, MPI_ERRORS_RETURN, ierr)
!
      allocate(dead_ids(0))
      call ft_join_after_shrink(shrunk_world, .true., 0, dead_ids)
      deallocate(dead_ids)
!
    endsubroutine spare_wait_loop
!***********************************************************************
    subroutine ft_check(ierr, comm_to_revoke)
!
!  comm_to_revoke: if present, MPIX_Comm_revoke is called on it before
!  ft_repair. The first rank to detect a failure (e.g. on a ghost-cell
!  exchange) passes the comm in. Revoke propagates the failure to every
!  other rank's next op on that comm, so non-neighbour ranks discover
!  the failure without a per-timestep consensus allreduce.
!
      integer, intent(in) :: ierr
      integer, intent(in), optional :: comm_to_revoke
      integer :: rev_ierr
!
      if (ierr == MPI_SUCCESS) return
      if (ierr == MPIX_ERR_PROC_FAILED .or. &
          ierr == MPIX_ERR_PROC_FAILED_PENDING .or. &
          ierr == MPIX_ERR_REVOKED) then
        if (present(comm_to_revoke)) &
          call MPIX_Comm_revoke(comm_to_revoke, rev_ierr)
        call ft_log('ft_check: detected ULFM error, invoking ft_repair')
        call ft_repair()
        pc_need_restart_step = .true.
      endif
!
    endsubroutine ft_check
!***********************************************************************
    subroutine ft_repair()
!
      integer :: ierr, shrunk_world
      integer :: failed_group, active_group, n_dead, j
      integer, allocatable :: dead_ids(:), src_ranks(:)
!
      call MPIX_Comm_get_failed(pc_comm_active, failed_group, ierr)
      call MPI_Group_size(failed_group, n_dead, ierr)
!
      allocate(dead_ids(max(1, n_dead)))
      if (n_dead > 0) then
        call MPI_Comm_group(pc_comm_active, active_group, ierr)
        allocate(src_ranks(n_dead))
        do j = 1, n_dead
          src_ranks(j) = j - 1
        enddo
        call MPI_Group_translate_ranks(failed_group, n_dead, src_ranks, &
                                       active_group, dead_ids, ierr)
        deallocate(src_ranks)
        call MPI_Group_free(active_group, ierr)
      endif
      call MPI_Group_free(failed_group, ierr)
!
      call MPIX_Comm_revoke(pc_comm_active, ierr)
!
      call MPIX_Comm_shrink(MPI_COMM_WORLD, shrunk_world, ierr)
      call MPI_Comm_set_errhandler(shrunk_world, MPI_ERRORS_RETURN, ierr)
!
      call ft_join_after_shrink(shrunk_world, .false., n_dead, dead_ids)
!
      deallocate(dead_ids)
!
    endsubroutine ft_repair
!***********************************************************************
    subroutine ft_join_after_shrink(shrunk_world, came_from_spare, n_dead_in, dead_ids_in)
!
      integer, intent(in) :: shrunk_world
      logical, intent(in) :: came_from_spare
      integer, intent(in) :: n_dead_in
      integer, intent(in) :: dead_ids_in(:)
      integer :: ierr, sw_size, sw_rank
      integer :: my_key, agree_val, i, n_dead, ncpus, n_survivors, spare_idx
      integer, allocatable :: dead_ids(:)
      integer :: new_pc_active
!
      call MPI_Comm_size(shrunk_world, sw_size, ierr)
      call MPI_Comm_rank(shrunk_world, sw_rank, ierr)
!
      if (came_from_spare) then
        n_dead = 0
      else
        n_dead = n_dead_in
      endif
      call MPI_Bcast(n_dead, 1, MPI_INTEGER, 0, shrunk_world, ierr)
!
      allocate(dead_ids(max(1, n_dead)))
      if (.not. came_from_spare .and. n_dead > 0) dead_ids(1:n_dead) = dead_ids_in(1:n_dead)
      if (n_dead > 0) call MPI_Bcast(dead_ids, n_dead, MPI_INTEGER, 0, shrunk_world, ierr)
!
      ncpus       = pc_world_size - pc_n_spares
      n_survivors = ncpus - n_dead
!
      if (.not. came_from_spare) then
        my_key = pc_active_rank
      else
        spare_idx = sw_rank - n_survivors
        if (spare_idx >= 0 .and. spare_idx < n_dead) then
          my_key = dead_ids(spare_idx + 1)
          pc_was_promoted = .true.
        else
          my_key = pc_world_size
        endif
      endif
!
      agree_val = n_dead
      call MPIX_Comm_agree(shrunk_world, agree_val, ierr)
!
      call MPI_Comm_split(shrunk_world, 0, my_key, new_pc_active, ierr)
      call MPI_Comm_set_errhandler(new_pc_active, MPI_ERRORS_RETURN, ierr)
!
      if (pc_comm_active /= MPI_COMM_NULL .and. .not. came_from_spare) &
        call MPI_Comm_free(pc_comm_active, ierr)
      if (came_from_spare .and. pc_comm_spare /= MPI_COMM_NULL) then
        call MPI_Comm_free(pc_comm_spare, ierr)
        pc_comm_spare = MPI_COMM_NULL
      endif
      call MPI_Comm_free(shrunk_world, ierr)
!
      pc_comm_active = new_pc_active
      call MPI_Comm_size(pc_comm_active, pc_active_size, ierr)
      call MPI_Comm_rank(pc_comm_active, pc_active_rank, ierr)
!
      if (came_from_spare) then
        pc_is_spare          = .false.
        pc_was_promoted      = .true.
        pc_need_restart_step = .true.
        call ft_log_pair('ft_repair: promoted spare ', pc_world_rank, &
                         ' -> rank ', my_key)
      else
        if (sw_rank == 0) then
          do i = 1, n_dead
            call ft_log_pair('ft_repair: dead cart id ', dead_ids(i), &
                             '; n_dead=', n_dead)
          enddo
        endif
        pc_need_restart_step = .true.
      endif
!
      if (associated(pc_neighbour_setup)) call pc_neighbour_setup(pc_comm_active)
      if (associated(pc_topology_verify)) call pc_topology_verify()
!
      deallocate(dead_ids)
!
    endsubroutine ft_join_after_shrink
!***********************************************************************
    subroutine ft_log(msg)
!
      character(len=*), intent(in) :: msg
!
      write(*,'(A,I0,A,A)') ' [ft] rank=', pc_world_rank, ' ', trim(msg)
      flush(6)
!
    endsubroutine ft_log
!***********************************************************************
    subroutine ft_log_pair(prefix, a, mid, b)
!
      character(len=*), intent(in) :: prefix, mid
      integer, intent(in) :: a, b
!
      write(*,'(A,I0,A,A,I0,A,I0)') ' [ft] rank=', pc_world_rank, ' ', &
                                     trim(prefix), a, trim(mid), b
      flush(6)
!
    endsubroutine ft_log_pair
!***********************************************************************
endmodule Ft_mpi
