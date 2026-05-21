! $Id$
!
!  Failure injector for ULFM testing.
!  Reads PENCIL_FT_KILL_RANK and PENCIL_FT_KILL_STEP env vars ONCE on the
!  first call and caches them; if the current rank matches kill_rank and
!  the current iteration matches kill_step, raises SIGKILL on itself via
!  the libc kill() syscall. Unset (or PENCIL_FT_KILL_STEP=-1) disables.
!
!  Caching matters for failure-free benchmarks: inject_failure(it) is
!  called every timestep from run.f90, so two getenv() + parse calls per
!  step would add unrelated overhead to the ULFM performance measurement.
!  After this change the per-step cost is two integer comparisons.
!
module Ft_inject
!
  use iso_c_binding, only: c_int
  use Ft_mpi, only: pc_world_rank
!
  implicit none
  private
!
  public :: inject_failure
!
  interface
    function c_kill(pid, sig) bind(c, name='kill') result(ret)
      import :: c_int
      integer(c_int), value :: pid, sig
      integer(c_int)        :: ret
    end function c_kill
    function c_getpid() bind(c, name='getpid') result(pid)
      import :: c_int
      integer(c_int) :: pid
    end function c_getpid
  end interface
!
  logical, save :: cached      = .false.
  integer, save :: kill_rank_c = -1
  integer, save :: kill_step_c = -1
!
contains
!***********************************************************************
    subroutine read_env_once()
!
      integer :: ios
      character(len=32) :: env
!
      if (cached) return
      cached = .true.
!
      call get_environment_variable('PENCIL_FT_KILL_RANK', env)
      if (len_trim(env) > 0) then
        read(env,*,iostat=ios) kill_rank_c
        if (ios /= 0) kill_rank_c = -1
      endif
      call get_environment_variable('PENCIL_FT_KILL_STEP', env)
      if (len_trim(env) > 0) then
        read(env,*,iostat=ios) kill_step_c
        if (ios /= 0) kill_step_c = -1
      endif
!
    endsubroutine read_env_once
!***********************************************************************
    subroutine inject_failure(it)
!
      integer, intent(in) :: it
      integer :: ret
!
      if (.not. cached) call read_env_once()
!
      if (kill_step_c < 0)             return     ! disabled
      if (pc_world_rank /= kill_rank_c) return
      if (it /= kill_step_c)           return
!
      write(*,'(A,I0,A,I0)') ' [ft_inject] rank=', pc_world_rank, ' SIGKILL at it=', it
      flush(6)
      ret = c_kill(c_getpid(), 9_c_int)
!
    endsubroutine inject_failure
!***********************************************************************
endmodule Ft_inject
