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

    !use m_data_output

    use m_bubbles_EL

    use m_bubbles_EL_kernels
    ! ==========================================================================

    implicit none

    type(scalar_field), allocatable, dimension(:) :: q_hifu !< HIFU vector field
    type(vector_field) :: q_hifu_3d
    !$acc declare create(q_hifu, q_hifu_3d)

    real(wp), allocatable, dimension(:) :: shear_viscous_fluids, bulk_viscous_fluids, abs_coef_fluids, rho_cp_fluids, tdiff_fluids
    !$acc declare create(shear_viscous_fluids, bulk_viscous_fluids, abs_coef_fluids, rho_cp_fluids, tdiff_fluids)

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

        ! Allocating the cell-average RHS variables
        @:ALLOCATE(q_hifu(1:sys_size_hifu))

        do i = 1, sys_size_hifu
            @:ALLOCATE(q_hifu(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_hifu(i))
        end do

        ! Fluids' properties needed (GPU)
        @:ALLOCATE(shear_viscous_fluids(1: num_fluids))
        @:ALLOCATE(bulk_viscous_fluids(1: num_fluids))
        @:ALLOCATE(abs_coef_fluids(1: num_fluids))
        @:ALLOCATE(rho_cp_fluids(1: num_fluids))
        @:ALLOCATE(tdiff_fluids(1: num_fluids))

        do i = 1, num_fluids
            shear_viscous_fluids(i) = fluid_pp(i)%Re(1)
            bulk_viscous_fluids(i) = fluid_pp(i)%Re(2)
            abs_coef_fluids(i) = fluid_pp(i)%absCoef
            rho_cp_fluids(i) = fluid_pp(i)%rho_cp
            tdiff_fluids(i) = fluid_pp(i)%tdiff
        end do
        !$acc update device(shear_viscous_fluids, bulk_viscous_fluids, abs_coef_fluids, rho_cp_fluids, tdiff_fluids)

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
        
        !!!!>>> NO GPU NEEDED <<<!!!!!!!

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
                p = hifu_params%p
                p_glb = p
                num_dims = 3
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3 -> restarting (2D -> 3D)'
            else
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3 -> restarting'
            end if

            return
        end if

    end subroutine s_restart_hifu_stages

    !> The idea is to jump form one stage into the other with this procedure
    subroutine s_HIFU_stages(t_step, hifu_write_output, exitFlag)

        integer, intent(inout) :: t_step
        logical, intent(out) :: hifu_write_output
        logical, intent(inout) :: exitFlag

        hifu_write_output = .false.

        ! 1st to 2nd stage
        if (cfl_dt) then
            if (mytime >= hifu_params%t_stop_stg1 .and. .not. hifu_params%sampling) then
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
                !$acc update device(hifu_params, dt)

                exitFlag = .false.
                return
            end if
        else
            if (t_step == hifu_params%t_step_stop_stg1 .and. .not. hifu_params%sampling) then
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
                !$acc update device(hifu_params, dt)

                exitFlag = .false.
                return
            end if
        end if

        ! 2nd to 3rd stage
        if (cfl_dt) then
            if (mytime >= hifu_params%t_stop_stg2 .and. .not. hifu_params%heatSolver) then
                ! Define params to start stage 3 (Move to constant dt)
                call s_close_run_time_information_samplingHIFU()
                if (.not. hifu_params%stg3) return

                hifu_params%sampling = .false.
                hifu_params%heatSolver = .true.
                cfl_dt = .false.
                mytime = 0._wp
                t_step_start = 0
                t_step = 0
                dt = hifu_params%dt_stg3
                t_step_save = hifu_params%t_step_save_stg3
                t_step_stop = hifu_params%t_step_stop_stg3
                finaltime = t_step_stop*dt

                if (hifu_params%stg3_3d) then
                    if (proc_rank == 0) then
                        print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D -> 3D)'
                        print '(" Simulating a ", A, " ", I0, "x", I0, "x", I0, " case on ", I0, " rank(s) ", A, ".")', &
#:if not MFC_CASE_OPTIMIZATION
                        "regular", &
#:else
                        "case-optimized", &
#:endif
                        m_glb, n_glb, hifu_params%p, num_procs, &
#ifdef MFC_OpenACC
!&<
                        "with OpenACC offloading"
!&>
#else
                        "on CPUs"
#endif
                    end if
                else
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D)'
                end if

                !$acc update device(hifu_params, dt)

                exitFlag = .false.

            end if
        else
            if (t_step == hifu_params%t_step_stop_stg2  .and. .not. hifu_params%heatSolver) then
                ! Define params to start stage 3 (Constant dt)
                call s_close_run_time_information_samplingHIFU()
                if (.not. hifu_params%stg3) return

                hifu_params%sampling = .false.
                hifu_params%heatSolver = .true.
                mytime = 0._wp
                t_step_start = 0
                t_step = 0
                dt = hifu_params%dt_stg3
                t_step_save = hifu_params%t_step_save_stg3
                t_step_stop = hifu_params%t_step_stop_stg3
                finaltime = t_step_stop*dt

                if (hifu_params%stg3_3d) then
                    if (proc_rank == 0) then
                        print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D -> 3D)'
                        print '(" Simulating a ", A, " ", I0, "x", I0, "x", I0, " case on ", I0, " rank(s) ", A, ".")', &
#:if not MFC_CASE_OPTIMIZATION
                        "regular", &
#:else
                        "case-optimized", &
#:endif
                        m_glb, n_glb, hifu_params%p, num_procs, &
#ifdef MFC_OpenACC
!&<
                        "with OpenACC offloading"
!&>
#else
                        "on CPUs"
#endif
                    end if
                else
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D)'
                end if

                !$acc update device(hifu_params, dt)

                exitFlag = .false.

            end if
        end if

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

        real(wp) :: rho_h, pres_h, gamma_h, pi_inf_h, T_h, c_c_h
        real(wp), dimension(num_dims) :: vel_h
        real(wp) :: qv_h, cson_h
        real(wp), dimension(2) :: Re_h
        real(wp), dimension(contxe) :: myalpha_rho, myalpha
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
        logical :: abortFlag

        focalIntensity_ac = 0._wp
        focalIntensity_ac_prms = 0._wp
        sumIntensity_ac = 0._wp

        if (cyl_coord .and. p == 0) then  !Axysimetric

            !$acc parallel loop collapse(3) gang vector default(present) reduction(+: sumIntensity_ac) &
            !$acc reduction(MAX: focalIntensity_ac, focalIntensity_ac_prms) private(myalpha_rho, myalpha, vel_h, Re_h, rhoYks_h) &
            !$acc  copy(sumIntensity_ac, focalIntensity_ac, focalIntensity_ac_prms)
            do l = 0, p
                do j = 0, m
                    do k = 0, n
                        abortFlag = .false.

                        !Get viscosities (and absorption coeff.) which are user inputs
                        shearVisc = 0._wp
                        bulkVisc = 0._wp
                        absCoef = 0._wp

                        !$acc loop seq
                        do i = 1, num_fluids
                            shearVisc = shearVisc + q_prim_vf(E_idx + i)%sf(j, k, l) * shear_viscous_fluids(i)
                            bulkVisc = bulkVisc + q_prim_vf(E_idx + i)%sf(j, k, l) * bulk_viscous_fluids(i)
                            absCoef = absCoef + q_prim_vf(E_idx + i)%sf(j, k, l) * abs_coef_fluids(i)
                        end do
                        shearVisc = 1._wp/shearVisc
                        bulkVisc = 1._wp/bulkVisc

                        if (f_is_default(absCoef)) stop "HIFU: Check absCoef values!"

                        !>> Get the strain rate tensor (using central finite difference)
                        varA = 0._wp
                        varB = 0._wp

                        ! Only for axysimmetric assumption
                        duxdx = (q_prim_vf(contxe + 1)%sf(j + 1, k, 0) - q_prim_vf(contxe + 1)%sf(j - 1, k, 0))/(x_cc(j + 1) - x_cc(j - 1))
                        duxdr = (q_prim_vf(contxe + 1)%sf(j, k + 1, 0) - q_prim_vf(contxe + 1)%sf(j, k - 1, 0))/(y_cc(k + 1) - y_cc(k - 1))

                        durdx = (q_prim_vf(contxe + 2)%sf(j + 1, k, 0) - q_prim_vf(contxe + 2)%sf(j - 1, k, 0))/(x_cc(j + 1) - x_cc(j - 1))
                        durdr = (q_prim_vf(contxe + 2)%sf(j, k + 1, 0) - q_prim_vf(contxe + 2)%sf(j, k - 1, 0))/(y_cc(k + 1) - y_cc(k - 1))

                        !>> Get pressure, density and speed of sound
                        !$acc loop seq
                        do i = 1, contxe
                            myalpha_rho(i) = q_prim_vf(advxb + i - 1)%sf(j, k, l)* &
                                            q_prim_vf(i)%sf(j, k, l)
                            myalpha(i) = q_prim_vf(advxb + i - 1)%sf(j, k, l)
                        end do

                        call s_convert_species_to_mixture_variables_acc(rho_h, gamma_h, pi_inf_h, qv_h, myalpha, &
                                                                myalpha_rho, Re_h, j, k, l)

                        !$acc loop seq
                        do s = 1, num_dims
                            vel_h(s) = q_cons_vf(s + contxe)%sf(j, k, l)/rho_h
                        end do

                        call s_compute_pressure(q_cons_vf(E_idx)%sf(j, k, l), 0._wp, 0.5_wp*rho_h*dot_product(vel_h, vel_h), &
                                                                        pi_inf_h, gamma_h, rho_h, qv_h, rhoYks_h, pres_h, T_h)

                        call s_compute_speed_of_sound(pres_h, rho_h, gamma_h, pi_inf_h, &
                                                      ((gamma_h + 1._wp)*pres_h + pi_inf_h)/rho_h, myalpha, 0._wp, c_c_h, cson_h)
                        
                        !Obtaining Pmax and Pmin fields
                        q_hifu(hifu_params%P_idx)%sf(j, k, l) = max(q_hifu(hifu_params%P_idx)%sf(j, k, l), pres_h)
                        q_hifu(hifu_params%P_idx + 1)%sf(j, k, l) = min(q_hifu(hifu_params%P_idx + 1)%sf(j, k, l), pres_h)

                        !>> Compute intensity form acoustic damping

                        ! PRMS method (calculate only during the last time step in stage2 -> need developed Pmax field)
                        intensity_ac_prms = 0._wp
                        if (cfl_dt) then
                            if (mytime >= t_stop) then
                                intensity_ac_prms = absCoef*(q_hifu(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h)
                            end if
                        else
                            if (t_step == t_step_stop - 1) then
                                intensity_ac_prms = absCoef*(q_hifu(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h)
                            end if
                        end if

                        ! Shear stress method
                        intensity_ac = 0._wp
                        ep11 = durdr
                        ep22 = vel_h(2)/y_cc(k)
                        ep33 = duxdx
                        ep13 = 0.5_wp*(durdx + duxdr)
                        varA = ep11**2._wp + ep22**2._wp + ep33**2._wp
                        varB = (8._wp/3._wp)*varA - (4._wp/3._wp)*(ep11*ep22 + ep11*ep33 + ep22*ep33) + 6._wp*(ep13**2._wp)
                        intensity_ac = intensity_ac + bulkVisc*varA + 2._wp*shearVisc*varB !intensity is "q_us_ac"

                        q_hifu(hifu_params%tsamp_idx)%sf(j, k, l) = q_hifu(hifu_params%tsamp_idx)%sf(j, k, l) &
                                                                                                        + hdid      ! Update total sampling time
                        q_hifu(hifu_params%qus_idx)%sf(j, k, l) = q_hifu(hifu_params%qus_idx)%sf(j, k, l) &
                                                                                            + intensity_ac*hdid     ! Sampling acoustic intensity
                        q_hifu(hifu_params%qus_prms_idx)%sf(j, k, l) = intensity_ac_prms * &
                                                                       q_hifu(hifu_params%tsamp_idx)%sf(j, k, l)    ! Sampling acoustic intensity (prms)
                        
                        ! Checking for NaNs
                        if (q_hifu(hifu_params%qus_idx)%sf(j, k, l) /= q_hifu(hifu_params%qus_idx)%sf(j, k, l)) then
                            print*, 'Acoustic intensity is NaN', j, k, l, hdid, intensity_ac
                            print*, 'viscosities (bulk & shear)', bulkVisc, shearVisc
                            print*, 'var: A, B', varA, varB, ep11, ep22, ep33, ep13
                            print*, 'ep22:', vel_h(2), y_cc(k), rho_h
                            abortFlag = .true.
                        end if

                        if (q_hifu(hifu_params%qus_prms_idx)%sf(j, k, l) /= q_hifu(hifu_params%qus_prms_idx)%sf(j, k, l)) then
                            print*, 'Acoustic intensity PRMS is NaN', j, k, l, hdid, intensity_ac_prms
                            print*, 'absCoef*(Pres - atmPres)**2/(rho_h*cson_h)', absCoef,q_hifu(hifu_params%P_idx)%sf(j, k, l), &
                                                                                                    hifu_params%atmPres, rho_h, cson_h
                            abortFlag = .true.
                        end if

                        if (abortFlag) stop "NaNs in Acoustic intensity (prms)"

                        !Update average velocities for streaming
                        q_hifu(hifu_params%u_idx)%sf(j, k, l) = q_hifu(hifu_params%u_idx)%sf(j, k, l) + vel_h(1)*hdid ! Sampling x-vel
                        q_hifu(hifu_params%v_idx)%sf(j, k, l) = q_hifu(hifu_params%v_idx)%sf(j, k, l) + vel_h(2)*hdid ! Sampling y-vel

                        !Get focal intensity and velocities
                        axialCondition = (dy(k) > y_cc(k) .and. y_cc(k) > 0._wp)
                        radialCondition = (x_cb(j - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < x_cb(j))
                        condition = (axialCondition .and. radialCondition)
                        if (condition) then
                            focalIntensity_ac = max(focalIntensity_ac, q_hifu(hifu_params%qus_idx)%sf(j, k, l))
                            focalIntensity_ac_prms = max(focalIntensity_ac_prms, q_hifu(hifu_params%qus_prms_idx)%sf(j, k, l))
                        end if

                        !Intensity summation through the domain
                        sumIntensity_ac = sumIntensity_ac + q_hifu(hifu_params%qus_idx)%sf(j, k, l)

                    end do
                end do
            end do

            if (num_procs > 1) then
                tmp = sumIntensity_ac
                call s_mpi_allreduce_sum(tmp, sumIntensity_ac)

                tmp = focalIntensity_ac
                call s_mpi_allreduce_max(tmp, focalIntensity_ac)

                tmp = focalIntensity_ac_prms
                call s_mpi_allreduce_max(tmp, focalIntensity_ac_prms)
            end if

            !$acc update host(q_hifu(hifu_params%tsamp_idx)%sf)

            if (proc_rank == 0) write (99, '(6x,5E24.8)') &
                                        mytime, &
                                        q_hifu(hifu_params%tsamp_idx)%sf(0, 0, 0), &
                                        focalIntensity_ac, &
                                        focalIntensity_ac_prms, &
                                        sumIntensity_ac
        else
            call s_mpi_abort('Getting HIFU samples (stage 2) works only with axisymmetric assumption so far!')
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
                    axialCondition = (dy(k) > y_cc(k) .and. y_cc(k) > 0._wp)
                    radialCondition = (x_cb(j - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < x_cb(j))
                    condition = (axialCondition .or. radialCondition)
                    if (condition) then
                        write (100, '(6x,5E24.8)') &
                            mytime, &
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
        real(wp) :: dz_val, z_max

        ! Azimutal direction range
        !z_max = min(2._wp*pi, abs(hifu_params%z_max))
        z_max = 2._wp*pi

        ! Deallocate some old vars
        do i = 1, sys_size_HIFU
            MPI_IO_HIFU_DATA%var(i)%sf => null()
        end do
        call s_finalize_mpi_proxy_module()

        p = hifu_params%p
        p_glb = p
        num_dims = 3

        if (bc_x%beg == -20) bc_x%beg = -6    ! from -20: acoustic bc
        bc_z%beg = -1; bc_z%end = -1          ! Assume entire cylindrical ring is taking care by one processor
        if (bc_y%beg == -2) bc_y%beg = -21    ! from -2: reflective boundary -> -21: cilyndrical sector (or)
                                              !                                 -14: full cylinder
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
        sys_size = 1
        do i = 1, sys_size
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
        call s_write_parallel_grid_zdir(z_max)

        !> Populate vars
        !Theta axis keeping same processor distribution
        dz_val = (z_max - 0._wp)/real(p + 1, wp)
        do i = 0, p
            z_cb(i - 1) = 0._wp + dz_val*real(i, wp)
        end do
        z_cb(p) = z_max
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

        ! 3d q_hifu
        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size_hifu
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_hifu_3d%vf(i)%sf(j, k, l) = q_hifu(i)%sf(j, k, 0)
                        if (i==1 .and. j==0 .and. k==0 .and. l==0) print*, q_hifu_3d%vf(i)%sf(j, k, l), q_hifu(i)%sf(j, k, 0), hifu_params%Tref
                    end do
                end do
            end do
        end do

        !Smear qvis and qth from the bubbles
        if (bubbles_lagrange) then 
            ! call s_acc_identify_bubble_angle_rotation()
            call s_smoothfunction(nBubs, intfc_rad, intfc_vel, &
                                                    mtn_s, mtn_pos, q_hifu_3d, bub_qvis, bub_qth)
        end if

        ! Add 3rd component of probe points (if any)
        if (probe_wrt) then
            do i = 1, num_probes
                probe(i)%z = 0._wp
            end do
        end if

    end subroutine s_initialize_from_2d_to_3d

    !> Identifies the correct angle to rotate the bubble cloud to obtain the maximum, minimum or average temperature profiles.
    !>      Every hifu_params%z_max rads, it calculates the summation of the viscous dissipation of each bubble weighted with the distance
    !>      bewtwen the transducer's focal point and the bubble's location
    ! subroutine s_acc_identify_bubble_angle_rotation()

    !     integer :: i, l, n_sector, numB_in
    !     real(wp), dimension(360) :: sumvals
    !     real(wp) :: distance, thetaPos, bound_down, bound_up, gpu_sum_vis, gpu_sum_th

    !     if (num_procs > 1) call s_mpi_abort('s_acc_identify_bubble_angle_rotation works with one processor only!')
    !     n_sector = 2*int(2._wp*pi/hifu_params%z_max) - 1

    !     !$acc parallel loop gang vector default(present) copyin(n_sector) private(sumvals)
    !     do i = 1, n_sector

    !         bound_down = real(i-1)*hifu_params%z_max*0.5_wp
    !         bound_up = bound_down + hifu_params%z_max 
    !         gpu_sum_vis = 0._wp
    !         gpu_sum_th = 0._wp
    !         numB_in = 0

    !         !$acc loop seq
    !         do l = 1, nBubs

    !             thetaPos = mtn_pos(l, 3, 1)
    !             if (thetaPos<0._wp) thetaPos = 2._wp*pi + mtn_pos(l, 3, 1)

    !             if (thetaPos >= bound_down .and. thetaPos<= bound_up) then
    !                 distance = sqrt((mtn_pos(l, 1, 1) - acoustic_bc_params%focLen)**2._wp + mtn_pos(l, 2, 1)**2._wp)
    !                 gpu_sum_vis = gpu_sum_vis + bub_qvis(l)*(1._wp/distance)
    !                 gpu_sum_th = gpu_sum_th + bub_qth(l)*(1._wp/distance)
    !                 numB_in = numB_in + 1
    !                 !if (i==1) print*, l, bub_qvis(l), distance, gpu_sum_vis
    !             end if
                
    !         end do

    !         sumvals(i) = gpu_sum_vis + gpu_sum_th

    !         print*, i, numB_in, gpu_sum_vis, sumvals(i)

    !     end do 

    ! end subroutine s_acc_identify_bubble_angle_rotation

    subroutine s_write_parallel_grid_zdir(z_max)

        real(wp), intent(in) :: z_max

#ifdef MFC_MPI

        ! Locations of cell boundaries
        real(wp), allocatable, dimension(:) :: z_cb_glb !<
            !! Locations of cell boundaries

        character(LEN=path_len + name_len) :: file_loc !<
            !! Generic string used to store the address of a file

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        real(wp) :: dz_val
        integer :: i !< Generic loop integers

        allocate (z_cb_glb(-1:p_glb))

        ! Grid generation in the z-direction
        if (p_glb > 0) then

            !Theta axis keeping same processor distribution
            dz_val = (z_max - 0._wp)/real(p_glb + 1, wp)
            do i = 0, p_glb
                z_cb_glb(i - 1) = 0._wp + dz_val*real(i, wp)
            end do
            z_cb_glb(p_glb) = z_max

            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
            data_size = p_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)
            call MPI_FILE_WRITE(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)

        end if

        deallocate (z_cb_glb)

#endif

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
        if (bc_x%beg == -6) bc_x%beg = -20                      ! from -20: acoustic bc
        bc_z%beg = dflt_int; bc_z%end = dflt_int                ! Assume entire cylindrical ring is taking care by one processor
        if (bc_y%beg == -14 .or. bc_y%beg == -21) bc_y%beg = -2 ! from -2: reflective boundary
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
            call s_mpi_barrier()

            !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht, t_step)
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
                            rho_cp = rho_cp + alpha * rho_cp_fluids(i)
                            tdiff = tdiff + alpha * tdiff_fluids(i)
                        end do

                        if (f_is_default(rho_cp) .or. f_is_default(tdiff)) then
                            print *, 'alpha, rho_cp, tdiff', alpha, rho_cp, tdiff
                            stop "HIFU: Check thermal properties!"
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

                        !Checking NaNs
                        if (q_hifu(hifu_params%T_idx + 1)%sf(j, k, l) /= q_hifu(hifu_params%T_idx + 1)%sf(j, k, l)) then
                            print*, 'NaNs in q hifu rhs', q_hifu(hifu_params%T_idx + 1)%sf(j, k, l), j, k, l
                            print*, 'Therm. properties', tdiff, rho_cp
                            stop "NaNs in q hifu rhs"
                        end if

                    end do
                end do
            end do

            !< 3D Cylindrical rhs
        else

            call s_populate_variables_buffers(q_hifu_3d%vf, pb, mv)

            !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht, t_step)
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
                            rho_cp = rho_cp + alpha * rho_cp_fluids(i)
                            tdiff = tdiff + alpha * tdiff_fluids(i)
                        end do

                        if (f_is_default(rho_cp) .or. f_is_default(tdiff)) then
                            print *, 'alpha, rho_cp, tdiff', alpha, rho_cp, tdiff
                            stop "HIFU: Check thermal properties!"
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
                                                                            q_hifu_3d%vf(qus_hifu_idx_ht)%sf(j, k, l) + &         !Acoustic intensity
                                                                            q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l) + &    !Viscous intensity
                                                                            q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l))        !Thermal intensity

                            if (hifu_params%streaming) then
                                !> Convected heat flux
                                stop "HIFU: No streaming valid for 3D heat solver!"
                            end if
                        end if

                        !Checking NaNs
                        if (q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) /= q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l)) then
                            print*, 'NaNs in q hifu rhs', q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l), j, k, l
                            print*, 'Therm. properties', tdiff, rho_cp
                            stop "NaNs in q hifu rhs"
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
        write (100, *) 'mytime, x_cc, y_cc, Pmax, Pmin'

        if (proc_rank == 0) then

            !Open files to save intensity sampling information at focus
            write (file_path, '(A)') '/D/sumIntensity-HIFU.dat'
            file_path = trim(case_dir)//trim(file_path)
            open (99, FILE=trim(file_path), FORM='formatted', POSITION='append', STATUS='unknown')
            write (99, *) 'mytime, numSamples, acousticFocalIntensity, acousticFocalIntensityPRMS, sumAcousticIntensity'

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

        @:DEALLOCATE(shear_viscous_fluids)
        @:DEALLOCATE(bulk_viscous_fluids)
        @:DEALLOCATE(abs_coef_fluids)
        @:DEALLOCATE(rho_cp_fluids)
        @:DEALLOCATE(tdiff_fluids)

    end subroutine s_finalize_HIFU_module

end module m_hifu
