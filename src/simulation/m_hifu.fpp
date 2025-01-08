!>
!! @file m_hifu.f90
!! @brief Contains module m_hifu

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to study HIFU
module m_hifu

    ! Dependencies ============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_bubbles_EL

    use m_bubbles_EL_kernels
    ! ==========================================================================

    implicit none

    type(scalar_field), allocatable, dimension(:), public :: q_hifu !< HIFU vector field
    type(vector_field), public :: q_hifu_3d
    !$acc declare create(q_hifu, q_hifu_3d)

contains

    !> Initializes the hifu model
    subroutine s_initialize_HIFU_module()

        integer :: i

        ! Define hifu indexes
        hifu_params%T_idx = 1
        hifu_params%tsamp_idx = 3
        hifu_params%qus_idx = 4
        hifu_params%qvis_idx = 5
        hifu_params%qth_idx = 7
        hifu_params%qus_prms_idx = 9
        hifu_params%P_idx = 10
        hifu_params%u_idx = 12
        hifu_params%v_idx = 14

        !$acc update device(hifu, hifu_params, sys_size_hifu)

        ! Allocating the cell-average RHS variables
        @:ALLOCATE(q_hifu(1:sys_size_hifu))

        do i = 1, sys_size_hifu
            @:ALLOCATE(q_hifu(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_hifu(i))
        end do

        !$acc update device(q_hifu)

    end subroutine s_initialize_HIFU_module

    !> Populate HIFU vars with user inputs and zeroing the time-averaged vars.
    subroutine s_start_HIFU_vars()

        integer :: i, j, k, l

        !Zeroing all the hifu variables

        !$acc parallel loop collapse(4) gang vector default(present)
        do l = 1, sys_size_hifu
            do k = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do i = idwbuff(1)%beg, idwbuff(1)%end
                        q_hifu(l)%sf(i, j, k) = 0._wp
                    end do
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = idwbuff(3)%beg, idwbuff(3)%end
            do j = idwbuff(2)%beg, idwbuff(2)%end
                do i = idwbuff(1)%beg, idwbuff(1)%end
                    !Initial Temperature
                    q_hifu(hifu_params%T_idx)%sf(i, j, k) = hifu_params%Tref
                    !Initialize Pmax
                    q_hifu(hifu_params%P_idx)%sf(i, j, k) = min(dflt_real, -dflt_real)
                    !Initialize Pmin
                    q_hifu(hifu_params%P_idx + 1)%sf(i, j, k) = max(dflt_real, -dflt_real)
                end do
            end do
        end do

        call s_open_run_time_information_samplingHIFU()

    end subroutine s_start_HIFU_vars

    subroutine s_restart_hifu_stages()

        ! Starting fresh
        if (cfl_dt) then
            if (n_start == 0) then
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 1: developing hydrodynamic field'
                return
            end if
        else
            if (t_step_start == 0) then
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 1: developing hydrodynamic field'
                return
            end if
        end if

        ! Restart during stg1
        if (hifu_params%stg1) then
            if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 1 -> restarting'
            return
        end if

        ! Restart during stg2
        if (hifu_params%stg2) then

            if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 2 -> restarting'
            hifu_params%sampling = .true.
            hifu_params%heatSolver = .false.

            if (cfl_dt) then
                t_stop = hifu_params%t_stop_stg2
            else
                t_step_stop = hifu_params%t_step_stop_stg2
            end if

            !$acc device update(hifu_params)

            return
        end if

        ! Restart during stg3
        if (hifu_params%stg3) then

            hifu_params%sampling = .false.
            hifu_params%heatSolver = .true.
            dt = hifu_params%dt_stg3

            cfl_dt = .false.
            t_step_save = hifu_params%t_step_save_stg3
            t_step_stop = hifu_params%t_step_stop_stg3
            if (hifu_params%stg3_3d) then
                p = p_hifu
                p_glb = p
                num_dims = 3
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3 -> restarting (2D -> 3D)'
            else
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3 -> restarting'
            end if

            !$acc device update(hifu_params, dt, p, num_dims)

            return
        end if

    end subroutine s_restart_hifu_stages

    !> The idea is to jump form one stage into the other with this procedure
    subroutine s_HIFU_stages(t_step, hifu_write_output)

        integer, intent(inout) :: t_step
        logical, intent(out) :: hifu_write_output

        integer :: save_count

        hifu_write_output = .false.

        ! 1st to 2nd stage
        if (cfl_dt) then
            if (mytime >= hifu_params%t_stop_stg1) then
                ! Define params to start stage 2 (Apadt dt)
                ! Stg 2 uses the same dt as in stg 1
                if (.not. hifu_params%stg2) return

                hifu_params%sampling = .true.
                hifu_params%heatSolver = .false.
                dt = hifu_params%dt_stg2
                t_stop = hifu_params%t_stop_stg2

                call s_start_HIFU_vars()
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 2: sampling heat sources'
                hifu_write_output = .true.
                return
            end if
        else
            if (t_step == hifu_params%t_step_stop_stg1) then
                ! Define params to start stage 2 (Constant dt)
                ! Stg 2 uses the same dt as in stg 1
                if (.not. hifu_params%stg2) return

                hifu_params%sampling = .true.
                hifu_params%heatSolver = .false.
                dt = hifu_params%dt_stg2
                t_step_stop = hifu_params%t_step_stop_stg2
                finaltime = t_step_stop*dt

                call s_start_HIFU_vars()
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 2: sampling heat sources'
                hifu_write_output = .true.
                return
            end if
        end if

        ! 2nd to 3rd stage
        if (cfl_dt) then
            if (mytime >= hifu_params%t_stop_stg2) then
                ! Define params to start stage 3 (Move to constant dt)
                call s_close_run_time_information_samplingHIFU()
                if (.not. hifu_params%stg3) return

                hifu_params%sampling = .false.
                hifu_params%heatSolver = .true.
                cfl_dt = .false.
                mytime = 0._wp
                t_step = 1
                dt = hifu_params%dt_stg3
                t_step_save = hifu_params%t_step_save_stg3
                t_step_stop = hifu_params%t_step_stop_stg3
                finaltime = t_step_stop*dt
                ! if (bc_x%beg == -20) bc_x%beg = -6 !In case acoustic BC is active
                if (hifu_params%stg3_3d) then
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D -> 3D)'
                else
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D)'
                end if

            end if
        else
            if (t_step == hifu_params%t_step_stop_stg2) then
                ! Define params to start stage 3 (Constant dt)
                call s_close_run_time_information_samplingHIFU()
                if (.not. hifu_params%stg3) return

                hifu_params%sampling = .false.
                hifu_params%heatSolver = .true.
                mytime = 0._wp
                t_step = 0
                dt = hifu_params%dt_stg3
                t_step_save = hifu_params%t_step_save_stg3
                t_step_stop = hifu_params%t_step_stop_stg3
                finaltime = t_step_stop*dt
                ! if (bc_x%beg == -20) bc_x%beg = -6 !In case acoustic BC is active
                if (hifu_params%stg3_3d) then
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D -> 3D)'
                else
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D)'
                end if

            end if
        end if

        !$acc device update(hifu_params, dt)

    end subroutine s_HIFU_stages

    !> The purpose of this procedure is to take samples needed to calculate the time averaged heat sources. It calculates
        !!      the generated heat source "q_us_ac", from the primary ultrasound source.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !!  @param t_step Current time step
        !!  @param hdid Physical period advanced in the last time step.
    subroutine s_update_HIFU_vars_sampling(q_cons_vf, q_prim_vf, t_step, hdid)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: t_step
        real(wp), intent(in) :: hdid

        real(wp) :: rho_h, pres_h, gamma_h, pi_inf_h, G_h, T_h, c_c_h
        real(wp), dimension(num_dims) :: vel_h
        real(wp) :: qv_h, cson_h
        real(wp), dimension(2) :: Re_h
        real(wp), dimension(num_fluids) :: alpha_h
        real(wp) :: rhoYks_h(1:num_species)

        logical :: axialCondition, radialCondition, condition
        real(wp) :: shearVisc, bulkVisc, absCoef
        real(wp) :: varA, varB
        real(wp) :: duxdx, duxdr, durdx, durdr, ep11, ep22, ep33, ep13
        real(wp) :: intensity_ac, sumIntensity_ac, tmp, focalIntensity_ac, intensity_ac_prms
        real(wp) :: focalIntensity_th, sumIntensity_th
        real(wp) :: sumIntensity_vis, focalIntensity_vis, focalIntensity_ac_prms
        real(wp) :: focal_u, focal_v

        integer :: i, j, k, l, s

        ! Zeroing out flow variables for all processors
        rho_h = 0._wp
        do s = 1, num_dims
            vel_h(s) = 0._wp
        end do
        pres_h = 0._wp
        cson_h = 0._wp
        gamma_h = 0._wp
        pi_inf_h = 0._wp
        sumIntensity_ac = 0._wp
        sumIntensity_vis = 0._wp
        sumIntensity_th = 0._wp

        if (cyl_coord .and. p == 0) then  !Axysimetric

            !$acc parallel loop collapse(4) gang vector default(present)
            !$acc reduction(+: sumIntensity_ac, sumIntensity_vis, sumIntensity_th)
            do l = 0, p
                do j = 0, m
                    do k = 0, n
                        !Get viscosities (and absorption coeff.) which are user inputs
                        shearVisc = 0._wp
                        bulkVisc = 0._wp
                        absCoef = 0._wp

                        !$acc loop seq
                        do i = 1, num_fluids
                            shearVisc = shearVisc + q_prim_vf(E_idx + i)%sf(j, k, l)*fluid_pp(i)%Re(1)
                            bulkVisc = bulkVisc + q_prim_vf(E_idx + i)%sf(j, k, l)*fluid_pp(i)%Re(2)
                            absCoef = absCoef + q_prim_vf(E_idx + i)%sf(j, k, l)*fluid_pp(i)%absCoef
                            alpha_h(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                        end do
                        shearVisc = 1/shearVisc
                        bulkVisc = 1/bulkVisc

                        if (f_is_default(absCoef)) call s_mpi_abort('Check absCoef values!')

                        !>> Get the strain rate tensor (using central finite difference)
                        varA = 0._wp
                        varB = 0._wp

                        ! Only for axysimmetric assumption
                        duxdx = (q_prim_vf(mom_idx%beg)%sf(j + 1, k, 0) - q_prim_vf(mom_idx%beg)%sf(j - 1, k, 0))/(x_cc(j + 1) - x_cc(j - 1))
                        duxdr = (q_prim_vf(mom_idx%beg)%sf(j, k + 1, 0) - q_prim_vf(mom_idx%beg)%sf(j, k - 1, 0))/(y_cc(k + 1) - y_cc(k - 1))

                        durdx = (q_prim_vf(mom_idx%beg + 1)%sf(j + 1, k, 0) - q_prim_vf(mom_idx%beg + 1)%sf(j - 1, k, 0))/(x_cc(j + 1) - x_cc(j - 1))
                        durdr = (q_prim_vf(mom_idx%beg + 1)%sf(j, k + 1, 0) - q_prim_vf(mom_idx%beg + 1)%sf(j, k - 1, 0))/(y_cc(k + 1) - y_cc(k - 1))

                        !>> Get pressure, density and speed of sound
                        call s_convert_to_mixture_variables(q_cons_vf, j, k, l, rho_h, gamma_h, pi_inf_h, qv_h, &
                                                            Re_h, G_h, fluid_pp(:)%G)
                        !$acc loop seq
                        do s = 1, num_dims
                            vel_h(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho_h
                        end do

                        call s_compute_pressure(q_cons_vf(E_idx)%sf(j, k, l), &
                                                0._wp, 0.5_wp*rho_h*dot_product(vel_h, vel_h), pi_inf_h, gamma_h, rho_h, qv_h, rhoYks_h, pres_h, T_h)

                        call s_compute_speed_of_sound(pres_h, rho_h, gamma_h, pi_inf_h, &
                                                      ((gamma_h + 1._wp)*pres_h + pi_inf_h)/rho_h, alpha_h, 0._wp, c_c_h, cson_h)

                        !>> Compute intensity form acoustic damping

                        intensity_ac_prms = 0._wp
                        intensity_ac_prms = intensity_ac_prms + absCoef*(q_hifu(hifu_params%P_idx)%sf(j, k, l) - &
                                                                         hifu_params%atmPres)**2/(rho_h*cson_h)

                        intensity_ac = 0._wp
                        ep11 = durdr
                        ep22 = vel_h(2)/y_cc(k)
                        ep33 = duxdx
                        ep13 = 0.5_wp*(durdx + duxdr)
                        varA = ep11**2.0 + ep22**2.0 + ep33**2.0
                        varB = (8._wp/3._wp)*varA - (4._wp/3._wp)*(ep11*ep22 + ep11*ep33 + ep22*ep33) + 6._wp*(ep13**2._wp)
                        intensity_ac = intensity_ac + bulkVisc*varA + 2._wp*shearVisc*varB !intensity is "q_us_ac"

                        q_hifu(hifu_params%tsamp_idx)%sf(j, k, l) = q_hifu(hifu_params%tsamp_idx)%sf(j, k, l) + &
                                                                    +hdid   ! Update total sampling time
                        q_hifu(hifu_params%qus_idx)%sf(j, k, l) = q_hifu(hifu_params%qus_idx)%sf(j, k, l) + &
                                                                  intensity_ac*hdid   ! Sampling acoustic intensity
                        q_hifu(hifu_params%qus_prms_idx)%sf(j, k, l) = q_hifu(hifu_params%qus_prms_idx)%sf(j, k, l) + &
                                                                       intensity_ac_prms*hdid   ! Sampling acoustic intensity (prms)

                        !Update average velocities for streaming
                        q_hifu(hifu_params%u_idx)%sf(j, k, l) = q_hifu(hifu_params%u_idx)%sf(j, k, l) + vel_h(1)*hdid ! Sampling x-vel
                        q_hifu(hifu_params%v_idx)%sf(j, k, l) = q_hifu(hifu_params%v_idx)%sf(j, k, l) + vel_h(2)*hdid ! Sampling y-vel

                        !Get focal intensity and velocities
                        axialCondition = (dy(k) > y_cc(k) .and. y_cc(k) > 0.0)
                        radialCondition = (x_cb(j - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < x_cb(j))
                        condition = (axialCondition .and. radialCondition)
                        if (condition) then
                            focalIntensity_ac = q_hifu(hifu_params%qus_idx)%sf(j, k, l)
                            focalIntensity_ac_prms = q_hifu(hifu_params%qus_prms_idx)%sf(j, k, l)
                            focalIntensity_vis = q_hifu(hifu_params%qvis_idx)%sf(j, k, l)
                            focalIntensity_th = q_hifu(hifu_params%qth_idx)%sf(j, k, l)
                            focal_u = q_hifu(hifu_params%u_idx)%sf(j, k, l)
                            focal_v = q_hifu(hifu_params%v_idx)%sf(j, k, l)
                        end if

                        !Intensity summation through the domain, avoid acoustic source influence (0.8*Focal length)
                        sumIntensity_ac = sumIntensity_ac + q_hifu(hifu_params%qus_idx)%sf(j, k, l)
                        sumIntensity_vis = sumIntensity_vis + q_hifu(hifu_params%qvis_idx)%sf(j, k, l)
                        sumIntensity_th = sumIntensity_th + q_hifu(hifu_params%qth_idx)%sf(j, k, l)

                        !Obtaining Pmax and Pmin fields
                        q_hifu(hifu_params%P_idx)%sf(j, k, l) = max(q_hifu(hifu_params%P_idx)%sf(j, k, l), pres_h)
                        q_hifu(hifu_params%P_idx + 1)%sf(j, k, l) = min(q_hifu(hifu_params%P_idx + 1)%sf(j, k, l), pres_h)

                    end do
                end do
            end do

            tmp = sumIntensity_ac
            call s_mpi_allreduce_sum(tmp, sumIntensity_ac)
            tmp = sumIntensity_vis
            call s_mpi_allreduce_sum(tmp, sumIntensity_vis)
            tmp = sumIntensity_th
            call s_mpi_allreduce_sum(tmp, sumIntensity_th)

            tmp = focalIntensity_vis
            call s_mpi_allreduce_max(tmp, focalIntensity_vis)
            tmp = focalIntensity_ac
            call s_mpi_allreduce_max(tmp, focalIntensity_ac)
            tmp = focalIntensity_ac_prms
            call s_mpi_allreduce_max(tmp, focalIntensity_ac_prms)

            tmp = focal_u
            call s_mpi_allreduce_sum(tmp, focal_u)
            tmp = focal_v
            call s_mpi_allreduce_sum(tmp, focal_v)

            if (proc_rank == 0) write (99, '(6x,I24.8,f24.8,9e24.8)') &
                t_step, q_hifu(hifu_params%tsamp_idx)%sf(0, 0, 0), focalIntensity_ac, &
                focalIntensity_ac_prms, focalIntensity_vis, focalIntensity_th, &
                sumIntensity_ac, sumIntensity_vis, sumIntensity_th, focal_u, focal_v

        end if

    end subroutine s_update_HIFU_vars_sampling

    !> The purpose of this procedure is to write the maximum and minimum pressure through time
        !!      along the axisymmetric and radial axes
        !! @param save_count File identifier
    subroutine s_write_Pmax(save_count)

        integer, intent(in) :: save_count

        integer :: j, k, l
        logical :: axialCondition, radialCondition, condition

        if (cyl_coord .and. p == 0) then
            l = 0
            do j = 0, m
                do k = 0, n
                    ! Specify enough conditions for axial and radial probe lines
                    axialCondition = (dy(k) > y_cc(k) .and. y_cc(k) > 0)
                    radialCondition = (x_cb(j - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < x_cb(j))
                    condition = (axialCondition .or. radialCondition)
                    if (condition) then
                        write (100, '(6x,I24,4E24.8)') &
                            save_count, &
                            x_cc(j), &
                            y_cc(k), &
                            q_hifu(hifu_params%P_idx)%sf(j, k, l), &
                            q_hifu(hifu_params%P_idx + 1)%sf(j, k, l)
                    end if
                end do
            end do
        end if

    end subroutine s_write_Pmax

    !Initilazile 3d cylindrical domain to solve heat equation.
    subroutine s_initialize_from_2d_to_3d()

        integer :: i, j, k, l
        real(wp) :: dz_val

        ! Deallocate some old vars
        do i = 1, sys_size_HIFU
            MPI_IO_HIFU_DATA%var(i)%sf => null()
        end do
        call s_finalize_mpi_proxy_module()

        p = p_hifu
        p_glb = p
        num_dims = 3

        if (bc_x%beg == -20) bc_x%beg = -6    ! from -20: acoustic bc
        bc_z%beg = -1; bc_z%end = -1        ! Assume entire cylindrical ring is taking care by one processor
        if (bc_y%beg == -2) bc_y%beg = -14 !   from -2: reflective boundary
        grid_geometry = 3

        idwbuff(3)%beg = -buff_size
        idwbuff(3)%end = p - idwbuff(3)%beg

        !$acc update device(p, num_dims, bc_x, bc_y, bc_z, grid_geometry, idwbuff)

        !> Allocate variables
        !Theta axis
        @:ALLOCATE(z_cb(-1 - buff_size:p + buff_size))
        @:ALLOCATE(z_cc(-buff_size:p + buff_size))
        @:ALLOCATE(dz(-buff_size:p + buff_size))

        !New mpi vars
        do i = 1, sys_size_HIFU
            allocate (MPI_IO_HIFU_DATA%var(i)%sf(0:m, 0:n, 0:p))
            MPI_IO_HIFU_DATA%var(i)%sf => null()
        end do
        call s_initialize_mpi_proxy_module()

        !Allocate 3d q_hifu
        @:ALLOCATE(q_hifu_3d%vf(1:sys_size_hifu))

        do i = 1, sys_size_hifu
            @:ALLOCATE(q_hifu_3d%vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
        end do
        @:ACC_SETUP_VFs(q_hifu_3d)

        !> Write grid z dir
        call s_write_parallel_grid_zdir()

        !> Populate vars
        !Theta axis keeping same processor distribution
        dz_val = (2._wp*pi - 0._wp)/real(p + 1, wp)
        do i = 0, p
            z_cb(i - 1) = 0._wp + dz_val*real(i, wp)
        end do
        z_cb(p) = 2._wp*pi
        ! Computing the cell width distribution
        dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
        ! Computing the cell center locations
        z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2._wp

        ! Population of Buffers in z-direction =============================

        ! Populating cell-width distribution buffer, at the beginning of the
        ! coordinate direction
        if (p == 0) then
            return
        elseif (bc_z%beg <= -3) then
            do i = 1, buff_size
                dz(-i) = dz(0)
            end do
        elseif (bc_z%beg == -2) then
            do i = 1, buff_size
                dz(-i) = dz(i - 1)
            end do
        elseif (bc_z%beg == -1) then
            do i = 1, buff_size
                dz(-i) = dz(p - (i - 1))
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(3, -1)
        end if

        ! Computing the cell-boundary locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            z_cb(-1 - i) = z_cb(-i) - dz(-i)
        end do
        ! Computing the cell-center locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            z_cc(-i) = z_cc(1 - i) - (dz(1 - i) + dz(-i))/2._wp
        end do

        ! Populating the cell-width distribution buffer, at the end of the
        ! coordinate direction
        if (bc_z%end <= -3) then
            do i = 1, buff_size
                dz(p + i) = dz(p)
            end do
        elseif (bc_z%end == -2) then
            do i = 1, buff_size
                dz(p + i) = dz(p - (i - 1))
            end do
        elseif (bc_z%end == -1) then
            do i = 1, buff_size
                dz(p + i) = dz(i - 1)
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(3, 1)
        end if

        ! Populating the cell-boundary locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            z_cb(p + i) = z_cb(p + (i - 1)) + dz(p + i)
        end do
        ! Populating the cell-center locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            z_cc(p + i) = z_cc(p + (i - 1)) + (dz(p + (i - 1)) + dz(p + i))/2._wp
        end do

        ! END: Population of Buffers in z-direction ========================

        !$acc update device(dz, z_cb, z_cc)

        !3d q_hifu
        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size_hifu
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_hifu_3d%vf(i)%sf(j, k, l) = q_hifu(i)%sf(j, k, 0)
                    end do
                end do
            end do
        end do

        !Smear qvis and qth from the bubbles
        if (bubbles_lagrange) call s_smoothfunction(nBubs, bub_qvis, bub_qth, &
                                                    mtn_s, mtn_pos, q_hifu_3d)

        ! Add 3rd component of probe points
        do i = 1, num_probes
            probe(i)%z = 0._wp
        end do

    end subroutine s_initialize_from_2d_to_3d

    subroutine s_write_parallel_grid_zdir()

        ! Locations of cell boundaries
        real(wp), allocatable, dimension(:) :: z_cb_glb !<
            !! Locations of cell boundaries

        character(LEN=path_len + name_len) :: file_loc !<
            !! Generic string used to store the address of a file

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        real(wp) :: dz_val
        integer :: i, j !< Generic loop integers

        allocate (z_cb_glb(-1:p_glb))

        ! Grid generation in the z-direction
        if (p_glb > 0) then

            !Theta axis keeping same processor distribution
            dz_val = (2._wp*pi - 0._wp)/real(p_glb + 1, wp)
            do i = 0, p_glb
                z_cb_glb(i - 1) = 0._wp + dz_val*real(i, wp)
            end do
            z_cb_glb(p_glb) = 2._wp*pi

            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
            data_size = p_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)
            call MPI_FILE_WRITE(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)

        end if

        deallocate (z_cb_glb)

    end subroutine s_write_parallel_grid_zdir

    !Finalize 3d cylindrical domain to solve heat equation.
    subroutine s_finalize_from_2d_to_3d()

        integer :: i

        if ((.not. hifu_params%stg1) .and. (.not. hifu_params%stg2)) return

        !Deallocate vars
        @:DEALLOCATE(z_cb, z_cc, dz)

        do i = 1, sys_size_hifu
            @:DEALLOCATE(q_hifu_3d%vf(i)%sf)
        end do
        @:DEALLOCATE(q_hifu_3d%vf)

        ! Restore 2D params
        p = 0
        p_glb = p
        num_dims = 2
        if (bc_x%beg == -6) bc_x%beg = -20            ! from -20: acoustic bc
        bc_z%beg = dflt_int; bc_z%end = dflt_int    ! Assume entire cylindrical ring is taking care by one processor
        if (bc_y%beg == -14) bc_y%beg = -2            ! from -2: reflective boundary
        grid_geometry = 2

        !$acc update device(p, num_dims, bc_x, bc_y, bc_z, grid_geometry)

    end subroutine s_finalize_from_2d_to_3d

    !Calculate the rhs value from heat transfer eqn discretized with finite volumes.
    subroutine s_rhs_heatEqn(q_cons_vf, pb, mv, t_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb
        real(wp), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: mv
        integer, intent(in) :: t_step

        real(wp) :: dTdx_L, dTdx_R, dTdr_L, dTdr_R
        real(wp) :: Tx_L, Tx_R, Tr_L, Tr_R
        real(wp) :: dTdz_L, dTdz_R, Tz_L, Tz_R
        real(wp) :: Ux_L, Ux_R, Ur_L, Ur_R
        real(wp) :: rho_cp, tdiff, alpha
        integer :: i, j, k, l, qus_hifu_idx_ht

        if (hifu_params%intPrms) then
            qus_hifu_idx_ht = hifu_params%qus_prms_idx
        else
            qus_hifu_idx_ht = hifu_params%qus_idx
        end if

        !< Axisymmetric rhs
        if (p == 0) then

            ! call s_populate_HIFU_variables_buffers(q_hifu)
            call s_populate_variables_buffers(q_hifu, pb, mv)

            !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht)
            do l = 0, p
                do j = 0, m
                    do k = 0, n

                        !<  Zeroing RHS_heat
                        q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) = 0._wp

                        !> Find temperature derivatives at the faces of the cell
                        dTdx_L = (q_hifu(hifu_params%T_idx)%sf(j, k, l) - q_hifu(hifu_params%T_idx)%sf(j - 1, k, l))/(x_cc(j) - x_cc(j - 1))
                        dTdx_R = (q_hifu(hifu_params%T_idx)%sf(j + 1, k, l) - q_hifu(hifu_params%T_idx)%sf(j, k, l))/(x_cc(j + 1) - x_cc(j))
                        dTdr_L = (q_hifu(hifu_params%T_idx)%sf(j, k, l) - q_hifu(hifu_params%T_idx)%sf(j, k - 1, l))/(y_cc(k) - y_cc(k - 1))
                        dTdr_R = (q_hifu(hifu_params%T_idx)%sf(j, k + 1, l) - q_hifu(hifu_params%T_idx)%sf(j, k, l))/(y_cc(k + 1) - y_cc(k))

                        !> Find temperature and streaming velocities at the faces of the cell
                        Tx_L = (q_hifu(hifu_params%T_idx)%sf(j, k, l) + q_hifu(hifu_params%T_idx)%sf(j - 1, k, l))/2._wp
                        Ux_L = (q_hifu(hifu_params%u_idx)%sf(j, k, l) + q_hifu(hifu_params%u_idx)%sf(j - 1, k, l))/2._wp
                        Tx_R = (q_hifu(hifu_params%T_idx)%sf(j, k, l) + q_hifu(hifu_params%T_idx)%sf(j + 1, k, l))/2._wp
                        Ux_R = (q_hifu(hifu_params%u_idx)%sf(j, k, l) + q_hifu(hifu_params%u_idx)%sf(j + 1, k, l))/2._wp
                        Tr_L = (q_hifu(hifu_params%T_idx)%sf(j, k, l) + q_hifu(hifu_params%T_idx)%sf(j, k - 1, l))/2._wp
                        Ur_L = (q_hifu(hifu_params%v_idx)%sf(j, k, l) + q_hifu(hifu_params%v_idx)%sf(j, k - 1, l))/2._wp
                        Tr_R = (q_hifu(hifu_params%T_idx)%sf(j, k, l) + q_hifu(hifu_params%T_idx)%sf(j, k + 1, l))/2._wp
                        Ur_R = (q_hifu(hifu_params%v_idx)%sf(j, k, l) + q_hifu(hifu_params%v_idx)%sf(j, k + 1, l))/2._wp

                        !> Get thermal properties
                        alpha = 0._wp
                        rho_cp = 0._wp
                        tdiff = 0._wp

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha = q_cons_vf(advxb + i - 1)%sf(j, k, 0)
                            rho_cp = rho_cp + alpha*fluid_pp(i)%rho_cp
                            tdiff = tdiff + alpha*fluid_pp(i)%tdiff
                        end do

                        if (f_is_default(rho_cp) .or. f_is_default(tdiff)) then
                            print *, 'alpha, rho_cp, tdiff', alpha, rho_cp, tdiff
                            call s_mpi_abort('Check thermal properties HIFU!')
                        end if

                        !> Obtain rhs (see notes)
                        q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) + &
                                                                    tdiff*(1._wp/dx(j))*(dTdx_R - dTdx_L) + &
                                                                    tdiff*(1._wp/(2._wp*y_cc(k)*dy(k)))*((2._wp*y_cc(k) + dy(k))*dTdr_R - (2._wp*y_cc(k) - dy(k))*dTdr_L)

                        if ((q_hifu(hifu_params%tsamp_idx)%sf(j, k, l) > 0._wp) .and. (t_step < hifu_params%stepStopSource)) then
                            !> Adding the heat source terms
                            q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) + &
                                                                        (1._wp/(rho_cp))*(1._wp/q_hifu(hifu_params%tsamp_idx)%sf(j, k, l))*q_hifu(qus_hifu_idx_ht)%sf(j, k, l) + &  !Acoustic intensity
                                                                        (1._wp/(rho_cp))*(1._wp/q_hifu(hifu_params%tsamp_idx)%sf(j, k, l))*q_hifu(hifu_params%qvis_idx)%sf(j, k, l) + &    !Viscous intensity
                                                                        (1._wp/(rho_cp))*(1._wp/q_hifu(hifu_params%tsamp_idx)%sf(j, k, l))*q_hifu(hifu_params%qth_idx)%sf(j, k, l)         !Thermal intensity

                            if (hifu_params%streaming) then
                                !> Convected heat flux
                                q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) - &
                                                                            (1._wp/q_hifu(hifu_params%tsamp_idx)%sf(j, k, l))*( & !Double check this
                                                                            (1._wp/dx(j))*(Ux_R*Tx_R - Ux_L*Tx_L) + &
                                                                            (1._wp/(2._wp*y_cc(k)*dy(k)))*((2._wp*y_cc(k) + dy(k))*Ur_R*Tr_R - (2._wp*y_cc(k) - dy(k))*Ur_L*Tr_L))
                            end if
                        end if
                    end do
                end do
            end do

            !< 3D Cylindrical rhs
        else

            call s_populate_variables_buffers(q_hifu_3d%vf, pb, mv)

            !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht)
            do l = 0, p
                do j = 0, m
                    do k = 0, n

                        !<  Zeroing RHS_heat
                        q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) = 0._wp

                        !> Find temperature derivatives at the faces of the cell
                        dTdx_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j - 1, k, l))/(x_cc(j) - x_cc(j - 1))
                        dTdx_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j + 1, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(x_cc(j + 1) - x_cc(j))
                        dTdr_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k - 1, l))/(y_cc(k) - y_cc(k - 1))
                        dTdr_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k + 1, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(y_cc(k + 1) - y_cc(k))
                        dTdz_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l - 1))/(z_cc(l) - z_cc(l - 1))
                        dTdz_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + 1) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(z_cc(l + 1) - z_cc(l))

                        !> Find temperature and streaming velocities at the faces of the cell
                        Tx_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu_3d%vf(hifu_params%T_idx)%sf(j - 1, k, l))/2._wp
                        Tx_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu_3d%vf(hifu_params%T_idx)%sf(j + 1, k, l))/2._wp
                        Tr_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k - 1, l))/2._wp
                        Tr_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k + 1, l))/2._wp
                        Tz_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l - 1))/2._wp
                        Tz_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + 1))/2._wp

                        !> Get thermal properties
                        alpha = 0._wp
                        rho_cp = 0._wp
                        tdiff = 0._wp

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha = q_cons_vf(advxb + i - 1)%sf(j, k, 0) !From 2D solution
                            rho_cp = rho_cp + alpha*fluid_pp(i)%rho_cp
                            tdiff = tdiff + alpha*fluid_pp(i)%tdiff
                        end do

                        if (f_is_default(rho_cp) .or. f_is_default(tdiff)) then
                            print *, 'alpha, rho_cp, tdiff', alpha, rho_cp, tdiff
                            call s_mpi_abort('Check thermal properties HIFU!')
                        end if

                        !> Obtain rhs (see notes)
                        q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) + tdiff*( &
                                                                          (1._wp/dx(j))*(dTdx_R - dTdx_L) + &
                                                                          (1._wp/(y_cc(k)*dy(k)))*((y_cc(k) + 0.5_wp*dy(k))*dTdr_R - &
                                                                                                   (y_cc(k) - 0.5_wp*dy(k))*dTdr_L) + &
                                                                          (1._wp/(dz(l)*y_cc(k)**2._wp))*(dTdz_R - dTdz_L))

                        if ((q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l) > 0._wp) .and. (t_step < hifu_params%stepStopSource)) then
                            !> Adding the heat source terms
                            q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) + &
                                                                              (1._wp/(rho_cp*q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l)))*( &
                                                                              q_hifu_3d%vf(qus_hifu_idx_ht)%sf(j, k, l) + &     !Acoustic intensity
                                                                              q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l) + &     !Viscous intensity
                                                                              q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l))     !Thermal intensity

                            if (hifu_params%streaming) then
                                !> Convected heat flux
                                call s_mpi_abort('No streaming valid for HIFU 3D heat eqn')
                            end if
                        end if
                    end do
                end do
            end do

        end if

        if (proc_rank == 0 .and. t_step == hifu_params%stepStopSource) print *, 'WARNING :: Turn off HIFU source'

    end subroutine s_rhs_heatEqn ! =============================================

    subroutine s_open_run_time_information_samplingHIFU()

        character(LEN=path_len + 3*name_len) :: file_path

        !Open files to save Pmax data at the axial and radial axes
        write (file_path, '(A,I0,A)') '/D/Pmax_', proc_rank, '.dat'
        file_path = trim(case_dir)//trim(file_path)
        open (100, FILE=trim(file_path), FORM='formatted', STATUS='unknown')
        write (100, *) 'timeStep, x_cc, y_cc, Pmax, Pmin, j, k, l'

        if (proc_rank == 0) then

            !Open files to save intensity sampling information at focus
            write (file_path, '(A)') '/D/sumIntensity-HIFU.dat'
            file_path = trim(case_dir)//trim(file_path)
            open (99, FILE=trim(file_path), FORM='formatted', POSITION='append', STATUS='unknown')
            write (99, *) 'timeStep, totalSamplingTime, acousticFocalIntensity, acousticFocalIntensityPRMS, ', &
                'viscousFocalIntensity, thermalFocalIntensity, sumAcousticIntensity, ', &
                'sumViscousIntensity, sumThermalIntensity, focalxVel, focalyVel'

            !Open files to save viscous and thermal intensity sampling information for a single bubble
            write (file_path, '(A,I0,A)') '/D/viscous_thermal_kernel-HIFU_', proc_rank, '.dat'
            file_path = trim(case_dir)//trim(file_path)
            open (98, FILE=trim(file_path), FORM='formatted', POSITION='append', STATUS='unknown')
            write (98, *) 'Recommended to use only with one particle to test and compare the performance of the smootheing function'
            write (98, *) 'Requires to uncomment some command lines in s_update_RK (m_particle.fpp)'
            write (98, *) 'dt_did, totalSamplingTime, viscousIntensity_beforeKernel, viscousIntensity_afterKernel, ', &
                'thermalIntensity_beforeKernel, thermalIntensity_afterKernel, radius, velocity'

        end if

    end subroutine s_open_run_time_information_samplingHIFU

    subroutine s_close_run_time_information_samplingHIFU()

        !Close files to save Pmax data at the axial and radial axes
        close (100)

        !Close file to save intensity sampling information at focus
        if (proc_rank == 0) close (99)

        !Close file to save viscous and thermal intensity sampling information for a single bubble
        close (98)

    end subroutine s_close_run_time_information_samplingHIFU

    subroutine s_finalize_HIFU_module()

        integer :: i

        do i = 1, sys_size_hifu
            @:DEALLOCATE(q_hifu(i)%sf)
        end do

    end subroutine s_finalize_HIFU_module

end module m_hifu
