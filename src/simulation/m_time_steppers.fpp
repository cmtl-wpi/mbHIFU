!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

#:include 'macros.fpp'

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_ibm

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_boundary_conditions

    use m_helper

    use m_sim_helpers

    use m_fftw

    use m_nvtx

    use m_thermochem

    use m_body_forces

    use m_particles             !< Lagrangian solver

    use m_kernel_functions

    use m_hifu                 !< HIFU

    ! ==========================================================================

    implicit none

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), q_cons_ts)
    !! Cell-average conservative variables at each time-stage (TS)

    @:CRAY_DECLARE_GLOBAL(type(scalar_field), dimension(:), q_prim_vf)
    !! Cell-average primitive variables at the current time-stage

    @:CRAY_DECLARE_GLOBAL(type(scalar_field), dimension(:), rhs_vf)
    !! Cell-average RHS variables at the current time-stage

    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), q_prim_ts)
    !! Cell-average primitive variables at consecutive TIMESTEPS

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :, :), rhs_pb)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :, :), rhs_mv)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension( :, :, :), max_dt)

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

    !$acc declare link(q_cons_ts,q_prim_vf,rhs_vf,q_prim_ts, rhs_mv, rhs_pb, max_dt)
#else
    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_pb

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_mv

    real(kind(0d0)), allocatable, dimension(:, :, :) :: max_dt

    integer, private :: num_ts, num_ts_hifu !<
    !! Number of time stages in the time-stepping scheme

    !$acc declare create(q_cons_ts,q_prim_vf,rhs_vf,q_prim_ts, rhs_mv, rhs_pb, max_dt)
#endif

    type(vector_field), allocatable, dimension(:) :: rhs_vp_adapt
    !! Adaptive 4th and 5th order Runge-Kutta-Cash-Karp time stepper

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module

        type(int_bounds_info) :: ix_t, iy_t, iz_t !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (coupledflag) then !Euler-Lagrangian solver
            num_ts = 2
        else
            if (time_stepper == 1) then
                num_ts = 1
            elseif (any(time_stepper == (/2, 3/))) then
                num_ts = 2
            end if
        end if

        if (hifu_intensityFlag .or. hifu_heateqnFlag) num_ts_hifu = 3

        ! Setting the indical bounds in the x-, y- and z-directions
        ix_t%beg = -buff_size; ix_t%end = m + buff_size

        if (n > 0) then
            iy_t%beg = -buff_size; iy_t%end = n + buff_size

            if (p > 0) then
                iz_t%beg = -buff_size; iz_t%end = p + buff_size
            else
                iz_t%beg = 0; iz_t%end = 0
            end if
        else
            iy_t%beg = 0; iy_t%end = 0
            iz_t%beg = 0; iz_t%end = 0
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE_GLOBAL(q_cons_ts(1:max(num_ts,num_ts_hifu)))

        do i = 1, num_ts
            @:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                @:ALLOCATE(q_cons_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
            end do
            @:ACC_SETUP_VFs(q_cons_ts(i))
        end do

        if (hifu_intensityFlag .or. hifu_heateqnFlag) then
            @:ALLOCATE(q_cons_ts(num_ts_hifu)%vf(1:max(sys_size,sys_size_hifu)))
            do j = 1, max(sys_size,sys_size_hifu)
                @:ALLOCATE(q_cons_ts(num_ts_hifu)%vf(j)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
            end do
        end if
        do i=1, num_ts_hifu
           if (proc_rank==0) print*,'DiegoV: q',i, 'size: ', size(q_cons_ts(i)%vf), sys_size
        end do        

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE_GLOBAL(q_prim_ts(0:3))

            do i = 0, 3
                @:ALLOCATE(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                        iy_t%beg:iy_t%end, &
                        iz_t%beg:iz_t%end))
                end do
            end do

            do i = 0, 3
                @:ACC_SETUP_VFs(q_prim_ts(i))
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE_GLOBAL(q_prim_vf(1:sys_size))

        do i = 1, adv_idx%end
            @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end))
            @:ACC_SETUP_SFs(q_prim_vf(i))
        end do

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
            if (adv_n) then
                @:ALLOCATE(q_prim_vf(n_idx)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(n_idx))
            end if
        end if

        if (hypoelasticity) then

            do i = stress_idx%beg, stress_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
        end if

        if (sigma /= dflt_real) then
            @:ALLOCATE(q_prim_vf(c_idx)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end))
            @:ACC_SETUP_SFs(q_prim_vf(c_idx))
        end if

        if (chemistry) then
            do i = chemxb, chemxe
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do

            @:ALLOCATE(q_prim_vf(tempxb)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end))
            @:ACC_SETUP_SFs(q_prim_vf(tempxb))
        end if

        @:ALLOCATE_GLOBAL(pb_ts(1:2))
        !Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE_GLOBAL(rhs_pb(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE_GLOBAL(rhs_pb(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
        end if

        @:ALLOCATE_GLOBAL(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE_GLOBAL(rhs_mv(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))

        else if (qbmm .and. polytropic) then
            @:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE_GLOBAL(rhs_mv(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
        end if

        ! Allocating the cell-average RHS variables
        @:ALLOCATE_GLOBAL(rhs_vf(1:sys_size))

        do i = 1, sys_size
            @:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(rhs_vf(i))
        end do

        ! Allocating the cell-average RHS variable for adaptive method, Lagrangian solver
        if (coupledflag .or. (solverapproach == 2)) then
            @:ALLOCATE_GLOBAL(rhs_vp_adapt(1:6))
            do i = 1, 6
                @:ALLOCATE(rhs_vp_adapt(i)%vf(1:sys_size))
            end do
            do i = 1, 6
                do j = 1, sys_size
                    @:ALLOCATE(rhs_vp_adapt(i)%vf(j)%sf(0:m,0:n,0:p))
                end do
                @:ACC_SETUP_SFs(rhs_vp_adapt(i))
            end do
        end if

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

        if (cfl_dt) then
            @:ALLOCATE_GLOBAL(max_dt(0:m, 0:n, 0:p))
        end if

    end subroutine s_initialize_time_steppers_module

    !Forward finite difference approximation for dT/dt
    subroutine s_time_stepper_heatEqn(t_step)

        integer, intent(IN) :: t_step

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish
        integer :: unitFile
        character(LEN=path_len + 3*name_len) :: file_path !<
        logical :: axialCondition, radialCondition_extra, radialCondition, condition
        real(kind(0d0)) :: xloc_focal

        ! Stage 1 of 1 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        if (t_step == t_step_start .and. proc_rank==0) print*, 'HIFU simulation >>>> Stage 3: Finding the final temperature distribution'

        call s_rhs_heatEqn(q_cons_ts(3)%vf, q_cons_ts(1)%vf, t_step)

        if (t_step == t_step_stop) return

        l = 0
        do  j = 0, m
            do k = 0, n
  
                !Forward euler time scheme, explicit
                q_cons_ts(3)%vf(T_hifu_idx)%sf(j, k, l) = q_cons_ts(3)%vf(T_hifu_idx)%sf(j, k, l) + q_cons_ts(3)%vf(T_hifu_idx+1)%sf(j, k, l)*dt

                if (ieee_is_nan(q_cons_ts(3)%vf(T_hifu_idx)%sf(j, k, l))) then
                    call s_mpi_abort('Temperature value is NaN!!')
                end if

            end do
        end do

        call nvtxEndRange

        call cpu_time(finish)

        ! ============================================================================

    end subroutine s_time_stepper_heatEqn

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step, time_avg)

        integer, intent(in) :: t_step
        real(kind(0d0)), intent(inout) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: nR3bar
        real(kind(0d0)) :: e_mix

        real(kind(0d0)) :: T
        real(kind(0d0)), dimension(num_species) :: Ys

        ! Stage 1 of 1 =====================================================

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg)

        if (ib .and. t_step == 1) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        call nvtxStartRange("body_forces")
        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, dt)
        call nvtxEndRange

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        ! ==================================================================

    end subroutine s_1st_order_tvd_rk

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step, time_avg)

        integer, intent(in) :: t_step
        real(kind(0d0)), intent(inout) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: nR3bar

        ! Stage 1 of 2 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg)

        if (ib .and. t_step == 1) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        call nvtxStartRange("body_forces")
        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, dt)
        call nvtxEndRange

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if
        ! ==================================================================

        ! Stage 2 of 2 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        call nvtxStartRange("body_forces")
        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, 2d0*dt/3d0)
        call nvtxEndRange

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)
        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q !< Generic loop iterator
        real(kind(0d0)) :: ts_error, denom, error_fraction, time_step_factor !< Generic loop iterator
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: nR3bar

        ! Stage 1 of 3 =====================================================

        if (.not. adap_dt) then
            call cpu_time(start)
            call nvtxStartRange("Time_Step")
        end if

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (hifu_intensityFlag) then !HIFU obtain intentities
            call s_update_HIFU_vars_stg2(q_cons_ts(1)%vf, q_prim_vf, q_cons_ts(3)%vf, t_step, dt)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        call nvtxStartRange("body_forces")
        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, dt)
        call nvtxEndRange

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if
        ! ==================================================================

        ! Stage 2 of 3 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            (3d0*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/4d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        call nvtxStartRange("body_forces")
        if (bodyForces) call s_apply_bodyforces(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, dt/4d0)
        call nvtxEndRange

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if
        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + 2d0*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + 2d0*dt*rhs_vf(i)%sf(j, k, l))/3d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*pb_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt*rhs_pb(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*mv_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt*rhs_mv(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        call nvtxStartRange("body_forces")
        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, 2d0*dt/3d0)
        call nvtxEndRange

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        if (.not. adap_dt) then
            call nvtxEndRange
            call cpu_time(finish)

            time = time + (finish - start)
        end if
        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk

    !> Strang splitting scheme with 3rd order TVD RK time-stepping algorithm for
        !!      the flux term and adaptive time stepping algorithm for
        !!      the source term
        !! @param t_step Current time-step
    subroutine s_strang_splitting(t_step, time_avg)

        integer, intent(in) :: t_step
        real(kind(0d0)), intent(inout) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        ! Stage 1 of 3 =====================================================
        call s_adaptive_dt_bubble(t_step)

        ! Stage 2 of 3 =====================================================
        call s_3rd_order_tvd_rk(t_step, time_avg)

        ! Stage 3 of 3 =====================================================
        call s_adaptive_dt_bubble(t_step)

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

        ! ==================================================================

    end subroutine s_strang_splitting

    !> Bubble source part in Strang operator splitting scheme
        !! @param t_step Current time-step
    subroutine s_adaptive_dt_bubble(t_step)

        integer, intent(in) :: t_step

        type(int_bounds_info) :: ix, iy, iz
        type(vector_field) :: gm_alpha_qp

        integer :: i, j, k, l, q !< Generic loop iterator

        ix%beg = 0; iy%beg = 0; iz%beg = 0
        ix%end = m; iy%end = n; iz%end = p
        call s_convert_conservative_to_primitive_variables( &
            q_cons_ts(1)%vf, &
            q_prim_vf, &
            gm_alpha_qp%vf, &
            ix, iy, iz)

        call s_compute_bubble_source(q_cons_ts(1)%vf, q_prim_vf, t_step, rhs_vf)

    end subroutine s_adaptive_dt_bubble

    subroutine s_compute_dt()

        real(kind(0d0)) :: rho        !< Cell-avg. density
        real(kind(0d0)), dimension(num_dims) :: vel        !< Cell-avg. velocity
        real(kind(0d0)) :: vel_sum    !< Cell-avg. velocity sum
        real(kind(0d0)) :: pres       !< Cell-avg. pressure
        real(kind(0d0)), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(kind(0d0)) :: gamma      !< Cell-avg. sp. heat ratio
        real(kind(0d0)) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(kind(0d0)) :: c          !< Cell-avg. sound speed
        real(kind(0d0)) :: H          !< Cell-avg. enthalpy
        real(kind(0d0)), dimension(2) :: Re         !< Cell-avg. Reynolds numbers
        type(vector_field) :: gm_alpha_qp
        real(kind(0d0)) :: dt_local
        type(int_bounds_info) :: ix, iy, iz
        integer :: i, j, k, l, q !< Generic loop iterators

        ix%beg = 0; iy%beg = 0; iz%beg = 0
        ix%end = m; iy%end = n; iz%end = p

        call s_convert_conservative_to_primitive_variables( &
            q_cons_ts(1)%vf, &
            q_prim_vf, &
            gm_alpha_qp%vf, &
            ix, iy, iz)

        !$acc parallel loop collapse(3) gang vector default(present) private(vel, alpha, Re)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)

                    ! Compute mixture sound speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, c)

                    call s_compute_dt_from_cfl(vel, c, max_dt, rho, Re, j, k, l)
                end do
            end do
        end do

        !$acc kernels
        dt_local = minval(max_dt)
        !$acc end kernels

        if (num_procs == 1) then
            dt = dt_local
        else
            call s_mpi_allreduce_min(dt_local, dt)
        end if

        !$acc update device(dt)

    end subroutine s_compute_dt

    !> This subroutine applies the body forces source term at each
        !! Runge-Kutta stage
    subroutine s_apply_bodyforces(q_cons_vf, q_prim_vf, rhs_vf, ldt)

        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(1:sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)), intent(in) :: ldt !< local dt

        integer :: i, j, k, l

        call s_compute_body_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) + &
                                                   ldt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

    end subroutine s_apply_bodyforces

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step)

        integer, intent(in) :: t_step

        integer :: i !< Generic loop iterator

        do i = 1, sys_size
            !$acc update host(q_prim_vf(i)%sf)
        end do

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling

    !> Cash-Karp Runge-Kutta 4th/5th order time-stepping algorithm
        !! @param realtime
        !! @param hnext
        !! @param hdid 
        !! @param t_step Current time-step
    subroutine rkqs(realtime, hnext, hdid, t_step)

        logical :: largestep
        real(kind(0.d0)) :: newtime, errmax, errmax_glb, qtime, hdid, hnext, dttarget
        real(kind(0.d0)) :: RKh, RKh_glb, htemp, SAFETY = 0.9d0, PGROW = -0.2d0, &
                            PSHRNK = -0.25d0, ERRCON = 1.89d-4
        integer :: i, j, k
        real(kind(0.d0)), intent(in) :: realtime
        integer, intent(in) :: t_step

        qtime = realtime
        dttarget = dt

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        !> Starting adaptive Runge-Kutta
        RKh = min(hnext, dttarget)
        RKh = max(Rkh, 1.0d-12)
        if (num_procs > 1) then
            call s_mpi_allreduce_min(RKh, RKh_glb)
            RKh = RKh_glb
        end if

        largestep = .false.
        if (coupledFlag .or. bubblesources) then
            call s_RK_particle_dynamics(qtime, 1, q_cons_ts(1)%vf, t_step, q_prim_vf, rhs_vp_adapt(1)%vf)
        else
            call s_RK_particle_dynamics(qtime, 1, q_cons_ts(1)%vf, t_step, q_prim_vf)
        end if

        !> Take a step
502     errmax = 0.0d0
        call rkck(qtime, RKh, errmax, largestep, t_step)

        if (largestep) then ! Negative radius, need to reduce time step
            if (cfl_dt) then
                if (RKh .gt. 1.0d-14) then
                    RKh = RKh/2.0d0
                    if (proc_rank==0) print*, '>>>>> WARNING: Reducing dt and restarting time step, now dt: ', RKh
                    largestep = .false.
                    goto 502
                else
                    call s_mpi_abort('Time step smaller than 1e-14')
                end if
            else
                call s_mpi_abort('Time step too large, please reduce dt or enable cfl_adapt_dt')
            end if
        end if

        if (cfl_dt) then !Check truncation error
            errmax = min(errmax,1.0d0)
            if (num_procs > 1) then
                call s_mpi_allreduce_max(errmax, errmax_glb)
                errmax=errmax_glb
            end if
            errmax=errmax/RKeps !Scale relative to USER required tolerance.
            if ((errmax .gt. 1.0d0)) then !Truncation error too large, reduce stepsize.
                htemp=SAFETY*RKh*(errmax**PSHRNK)
                RKh=sign(max(abs(htemp),0.1d0*abs(RKh)),RKh)  ! No more than a factor of 10.
                if (proc_rank==0) print*, '>>>>> WARNING: Truncation error found. Reducing dt and restaring time step, now dt: ', RKh
                goto 502         
            else ! Step succeeded. Compute size of next step.
                if (errmax .gt. ERRCON) then
                    hnext=SAFETY*RKh*(errmax**PGROW) ! No more than a factor of 5 increase. 
                else    
                    hnext=2.0d0*RKh !Truncation error too small (< 1.89e-4), increase time step
                end if    
            end if 
            hnext = min(hnext, dt0)
            
        else
            hnext = RKh
        end if
        
        !if (proc_rank==0) print*, hnext, RKh, errmax, PGROW, SAFETY
        hdid = RKh

        !> Update values
        qtime = qtime + hdid

        if (hifu_intensityFlag) then !HIFU obtain intentities
            call s_update_RK(q_cons_ts, .true., q_prim_vf, q_cons_ts(3)%vf, hdid)
            call s_update_HIFU_vars_stg2(q_cons_ts(1)%vf, q_prim_vf, q_cons_ts(3)%vf, t_step, hdid)
        else
            call s_update_RK(q_cons_ts, .true., q_prim_vf)
        end if

        if (avgdensflag) call s_write_void_evol(qtime)
        if (particlestatFlag) call s_calculate_particle_stats()

        return

    end subroutine rkqs

    !> Cash-Karp Runge-Kutta step
    subroutine rkck(qtime, RKh, errmax, largestep, t_step)
        !> USES derivs
        !> Given values for n variables y and their derivatives dydx known at x, use the .fth-order
        !> Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        !> the incremented variables as yout. Also return an estimate of the local truncation error
        !> in yout using the embedded fourth-order method. The user supplies the subroutine
        !> derivs(x,y,dydx), which returns derivatives dydx at x.

        logical :: largestep
        real(kind(0.d0)) :: RKh, qtime, errmax
        integer, intent(in) :: t_step
        integer :: i, j, k, l
        real(kind(0.d0)) :: A2 = 0.2d0, A3 = 0.3d0, A4 = 0.6d0, A5 = 1.0d0, A6 = 0.875d0
        real(kind(0.d0)), dimension(6) :: &
            RKcoef1 = (/0.2d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/), &
            RKcoef2 = (/3.0d0/40.0d0, 9.0d0/40.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/), &
            RKcoef3 = (/0.3d0, -0.9d0, 1.2d0, 0.0d0, 0.0d0, 0.0d0/), &
            RKcoef4 = (/-11.0d0/54.0d0, 2.5d0, -70.0d0/27.0d0, 35.d0/27.d0, 0.0d0, 0.0d0/), &
            RKcoef5 = (/1631.0d0/55296.0d0, 175.0d0/512.0d0, 575.d0/13824.d0, 44275.d0/110592.d0, 253.d0/4096.d0, 0.0d0/), &
            RKcoef6 = (/37.d0/378.d0, 0.0d0, 250.d0/621.d0, 125.0d0/594.0d0, 0.0d0, 512.0d0/1771.0d0/), &
            RKcoefE = (/37.d0/378.d0 - 2825.0d0/27648.0d0, 0.0d0, 250.d0/621.d0 - 18575.0d0/48384.0d0, &
                        125.0d0/594.0d0 - 13525.0d0/55296.0d0, -277.0d0/14336.0d0, 512.0d0/1771.0d0 - 0.25d0/)

        if (coupledFlag .or. bubblesources) then

            !> First step
            !if (proc_rank == 0) print *, 'rkqs 1st step at', qtime
            call s_update_particle(RKh, 1, RKcoef1, largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf, .true.)
            if (largestep) return

            !> Second step
            !if (proc_rank == 0) print *, 'rkqs 2nd step at', qtime + A2*RKh
            call s_RK_particle_dynamics(qtime + A2*RKh, 2, q_cons_ts(2)%vf, t_step, q_prim_vf, rhs_vp_adapt(2)%vf)
            call s_update_particle(RKh, 2, RKcoef2, largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)
            if (largestep) return

            !> Third step
            !if (proc_rank == 0) print *, 'rkqs 3rd step at', qtime + A3*RKh
            call s_RK_particle_dynamics(qtime + A3*RKh, 3, q_cons_ts(2)%vf, t_step, q_prim_vf, rhs_vp_adapt(3)%vf)
            call s_update_particle(RKh, 3, RKcoef3, largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)
            if (largestep) return

            !> Fourth step
            !if (proc_rank == 0) print *, 'rkqs 4th step at', qtime + A4*RKh
            call s_RK_particle_dynamics(qtime + A4*RKh, 4, q_cons_ts(2)%vf, t_step, q_prim_vf, rhs_vp_adapt(4)%vf)
            call s_update_particle(RKh, 4, RKcoef4, largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)
            if (largestep) return

            !> Fifth step
            !if (proc_rank == 0) print *, 'rkqs 5th step at', qtime + A5*RKh
            call s_RK_particle_dynamics(qtime + A5*RKh, 5, q_cons_ts(2)%vf, t_step, q_prim_vf, rhs_vp_adapt(5)%vf)
            call s_update_particle(RKh, 5, RKcoef5, largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)
            if (largestep) return

            !> Sixth step
            !if (proc_rank == 0) print *, 'rkqs 6th step at', qtime + A6*RKh
            call s_RK_particle_dynamics(qtime + A6*RKh, 6, q_cons_ts(2)%vf, t_step, q_prim_vf, rhs_vp_adapt(6)%vf)
            call s_update_particle(RKh, 6, RKcoef6, largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)
            if (largestep) return

            ! Configuring Coordinate Direction indexes =========================
            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

            if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            ! ==================================================================

            do i = 1, cont_idx%end
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            q_prim_vf(i)%sf(j, k, l) = q_cons_ts(1)%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            do i = adv_idx%beg, sys_size
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            q_prim_vf(i)%sf(j, k, l) = q_cons_ts(1)%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

            call s_calculate_RKerror(qtime + RKh, RKh, RKcoefE, errmax, t_step, q_cons_ts, q_prim_vf, rhs_vp_adapt)

        else

            !> First step
            call s_update_particle(RKh, 1, RKcoef1, largestep)
            if (largestep) return

            !> Second step
            call s_RK_particle_dynamics(qtime + A2*RKh, 2, q_cons_ts(1)%vf, t_step, q_prim_vf)
            call s_update_particle(RKh, 2, RKcoef2, largestep)
            if (largestep) return

            !> Third step
            call s_RK_particle_dynamics(qtime + A3*RKh, 3, q_cons_ts(1)%vf, t_step, q_prim_vf)
            call s_update_particle(RKh, 3, RKcoef3, largestep)
            if (largestep) return

            !> Fourth step
            call s_RK_particle_dynamics(qtime + A4*RKh, 4, q_cons_ts(1)%vf, t_step, q_prim_vf)
            call s_update_particle(RKh, 4, RKcoef4, largestep)
            if (largestep) return

            !> Fifth step
            call s_RK_particle_dynamics(qtime + A5*RKh, 5, q_cons_ts(1)%vf, t_step, q_prim_vf)
            call s_update_particle(RKh, 5, RKcoef5, largestep)
            if (largestep) return

            !> Sixth step
            call s_RK_particle_dynamics(qtime + A6*RKh, 6, q_cons_ts(1)%vf, t_step, q_prim_vf)
            call s_update_particle(RKh, 6, RKcoef6, largestep)
            if (largestep) return

            call s_calculate_RKerror(qtime + RKh, RKh, RKcoefE, errmax, t_step)

        end if

    end subroutine rkck

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, min(num_ts,num_ts_hifu)

            do j = 1, sys_size
                @:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
            end do

        end do

        if (hifu_intensityFlag .or. hifu_heateqnFlag) then
            do j = 1, max(sys_size,sys_size_hifu)
                @:DEALLOCATE(q_cons_ts(num_ts_hifu)%vf(j)%sf)
            end do
        end if

        do i = 1, max(num_ts,num_ts_hifu)
            @:DEALLOCATE(q_cons_ts(i)%vf)
        end do

        @:DEALLOCATE_GLOBAL(q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
                    @:DEALLOCATE(q_prim_ts(i)%vf(j)%sf)
                end do
                @:DEALLOCATE(q_prim_ts(i)%vf)
            end do
            @:DEALLOCATE_GLOBAL(q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = 1, adv_idx%end
            @:DEALLOCATE(q_prim_vf(i)%sf)
        end do

        if (hypoelasticity) then
            do i = stress_idx%beg, stress_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        @:DEALLOCATE_GLOBAL(q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            @:DEALLOCATE(rhs_vf(i)%sf)
        end do

        @:DEALLOCATE_GLOBAL(rhs_vf)

        ! Deallocating the cell-average RHS variable for adaptive method, Lagrangian solver
        if (coupledflag .or. (solverapproach == 2)) then
            do i = 1, 6
                do j = 1, adv_idx%end
                    deallocate (rhs_vp_adapt(i)%vf(j)%sf)
                end do
                deallocate (rhs_vp_adapt(i)%vf)
            end do
            deallocate (rhs_vp_adapt)
        end if

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

        if (hifu_intensityFlag) call s_close_run_time_information_samplingHIFU()

    end subroutine s_finalize_time_steppers_module

end module m_time_steppers
