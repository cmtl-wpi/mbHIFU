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

    type(vector_field) :: q_hifu, q_hifu_3d !< HIFU vector fields
    !$acc declare create(q_hifu, q_hifu_3d)

    real(wp), allocatable, dimension(:) :: shear_viscous_fluids, bulk_viscous_fluids, abs_coef_fluids, rho_cp_fluids, tdiff_fluids
    !$acc declare create(shear_viscous_fluids, bulk_viscous_fluids, abs_coef_fluids, rho_cp_fluids, tdiff_fluids)

    integer :: bc_pole
    !$acc declare create(bc_pole)

    integer :: sys_size_hyd

contains

    !> Initializes the hifu model
    subroutine s_initialize_HIFU_module()

        integer :: i

        sys_size_hyd = sys_size

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
        @:ALLOCATE(q_hifu%vf(1:sys_size_hifu))

        do i = 1, sys_size_hifu
            @:ALLOCATE(q_hifu%vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
        end do
        @:ACC_SETUP_VFs(q_hifu)

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
                        q_hifu%vf(l)%sf(i, j, k) = 0._wp
                    end do
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = idwbuff(3)%beg, idwbuff(3)%end
            do j = idwbuff(2)%beg, idwbuff(2)%end
                do i = idwbuff(1)%beg, idwbuff(1)%end
                    !Initial Temperature
                    q_hifu%vf(hifu_params%T_idx)%sf(i, j, k) = hifu_params%Tref
                    !Initialize Pmax
                    q_hifu%vf(hifu_params%P_idx)%sf(i, j, k) = min(dflt_real, -dflt_real)
                    !Initialize Pmin
                    q_hifu%vf(hifu_params%P_idx + 1)%sf(i, j, k) = max(dflt_real, -dflt_real)
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
            if (.not. f_is_default(hifu_params%dt_stg3)) then
                dt = hifu_params%dt_stg3
            else if (.not. f_is_default(hifu_params%cfl_stg3)) then
                dt = f_dt_from_CFL(hifu_params%cfl_stg3)
            else
                stop 'Define dt_stg3 or cfl_stg3'
            end if

            cfl_dt = .false.
            t_step_save = hifu_params%t_step_save_stg3
            t_step_stop = hifu_params%t_step_stop_stg3 - 1

            if (hifu_params%stg3_3d) then
                p = hifu_params%p
                p_glb = p
                num_dims = 3
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3 -> restarting (2D -> 3D)'
            else
                if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3 -> restarting'
            end if

            if (.not. hifu_params%cartesian .and. hifu_params%stg3_3d) call s_reduce_heat_domain()

            return
        end if

    end subroutine s_restart_hifu_stages

    function f_dt_from_CFL(cfl_heat)

        real(wp), intent(in) :: cfl_heat
        real(wp) :: f_dt_from_CFL, ds, tdiff, val_tmp

        integer :: i, j, k, l

        ds = abs(dflt_real)
        tdiff = abs(dflt_real)

        if (hifu_params%stg3_3d) then
            stop 'Need to be implemented f_dt_from_CFL'
        else

            do j = 0, m
                ds = min(ds, dx(j)) 
            end do
            do k = 0, n
                ds = min(ds, dy(k)) 
            end do
            if (p>0) then
                do l = 0, p
                    ds = min(ds, dz(l)) 
                end do
            end if

            do i = 1, num_fluids
                tdiff = min(tdiff, tdiff_fluids(i))
            end do

            f_dt_from_CFL = cfl_heat*(ds**2._wp)/tdiff
            if (p > 0) f_dt_from_CFL = f_dt_from_CFL

        end if

        if (num_procs > 1) then
            val_tmp = f_dt_from_CFL
            call s_mpi_allreduce_sum(val_tmp, f_dt_from_CFL)            
        end if

    end function f_dt_from_CFL

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
                if (.not. f_is_default(hifu_params%dt_stg3)) then
                    dt = hifu_params%dt_stg3
                else if (.not. f_is_default(hifu_params%cfl_stg3)) then
                    dt = f_dt_from_CFL(hifu_params%cfl_stg3)
                else
                    stop 'Define dt_stg3 or cfl_stg3'
                end if
                t_step_save = hifu_params%t_step_save_stg3
                t_step_stop = hifu_params%t_step_stop_stg3 - 1
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
                else if (p == 0) then
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D)'
                else if (p > 0) then
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (3D)'
                    if (bc_x%beg == -20) bc_x%beg = -6
                        !$acc update device(bc_x)
                end if

                if (.not. hifu_params%cartesian .and. hifu_params%stg3_3d) call s_reduce_heat_domain()

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

                if (.not. f_is_default(hifu_params%dt_stg3)) then
                    dt = hifu_params%dt_stg3
                else if (.not. f_is_default(hifu_params%cfl_stg3)) then
                    dt = f_dt_from_CFL(hifu_params%cfl_stg3)
                else
                    stop 'Define dt_stg3 or cfl_stg3'
                end if

                t_step_save = hifu_params%t_step_save_stg3
                t_step_stop = hifu_params%t_step_stop_stg3 - 1
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
                else if (p == 0) then
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (2D)'
                else if (p > 0) then
                    if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3: solving heat equation (3D)'
                    if (bc_x%beg == -20) bc_x%beg = -6
                        !$acc update device(bc_x)
                end if

                if (.not. hifu_params%cartesian .and. hifu_params%stg3_3d) call s_reduce_heat_domain()

                !$acc update device(hifu_params, dt)

                exitFlag = .false.

            end if
        end if

    end subroutine s_HIFU_stages


    subroutine s_reduce_heat_domain()

        !no GPU needed!!!!
        integer :: j, k, l
        real(wp) :: fun_tmp_old, fun_tmp, val_tmp, n_cell_xy


        ! Identify m_beg from x_beg
        hifu_params%mb = 0
        if (.not. f_is_default(hifu_params%xb)) then
            fun_tmp_old = x_cc(0) - hifu_params%xb
            do j = 0, m
                fun_tmp = x_cc(j) - hifu_params%xb
                if (fun_tmp * fun_tmp_old <= 0._wp) then
                    hifu_params%mb = j
                    exit
                else
                    fun_tmp_old = fun_tmp
                end if
            end do
        end if

        ! Identify m_end from x_end
        hifu_params%me = m
        if (.not. f_is_default(hifu_params%xe)) then
            fun_tmp_old = x_cc(0) - hifu_params%xe
            do j = 0, m
                fun_tmp = x_cc(j) - hifu_params%xe
                if (fun_tmp * fun_tmp_old <= 0._wp) then
                    hifu_params%me = j
                    exit
                else
                    fun_tmp_old = fun_tmp
                end if
            end do
        end if

        ! Identify n_end from y_end
        hifu_params%ne = n
        if (.not. f_is_default(hifu_params%ye)) then
            fun_tmp_old = y_cc(0) - hifu_params%ye
            do j = 0, n
                fun_tmp = y_cc(j) - hifu_params%ye
                if (fun_tmp * fun_tmp_old <= 0._wp) then
                    hifu_params%ne = j
                    exit
                else
                    fun_tmp_old = fun_tmp
                end if
            end do
        end if

        ! Find processors completely out the heat domain
        if (.not. f_is_default(hifu_params%xb)) then
            if (x_cc(m) < hifu_params%xb) hifu_params%me = 0
        end if
        if (.not. f_is_default(hifu_params%xe)) then
            if (x_cc(0) > hifu_params%xe) hifu_params%me = 0
        end if
        if (.not. f_is_default(hifu_params%ye)) then
            if (y_cc(0) > hifu_params%ye) hifu_params%ne = 0
        end if

        ! Actual number of cells in xy plane 
        n_cell_xy = real((hifu_params%me - hifu_params%mb+1)*(hifu_params%ne+1))

        if (num_procs>1) then
            val_tmp = n_cell_xy
            call s_mpi_allreduce_sum(val_tmp, n_cell_xy)
        end if

        if (proc_rank == 0) print '(" Domain reduction from ", I0, "x", I0, " to ", I0, "x", I0, ".")', &
                                    m_glb*n_glb, hifu_params%p, int(n_cell_xy), hifu_params%p

    end subroutine s_reduce_heat_domain

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
        real(wp) :: varA, varB, varC
        real(wp) :: duxdx, duxdr, durdx, durdr, ep11, ep22, ep33, ep12, ep13, ep23
        real(wp), dimension(3) :: duxdn, duydn, duzdn
        real(wp) :: intensity_ac, sumIntensity_ac, tmp, focalIntensity_ac, intensity_ac_prms
        real(wp) :: focalIntensity_th, sumIntensity_th
        real(wp) :: sumIntensity_vis, focalIntensity_vis, focalIntensity_ac_prms
        real(wp) :: focal_u, focal_v

        integer :: i, j, k, l, s, mtd_idx
        integer :: abortFlag, abortFlag_max

        focalIntensity_ac = 0._wp
        focalIntensity_ac_prms = 0._wp
        sumIntensity_ac = 0._wp

        if (bubbles_lagrange .and. .not. adap_dt) call s_compute_bubble_heat_sources_HIFU(hdid)

        abortFlag_max = 0

        if (cyl_coord .and. p == 0) then  !Axysimetric

#ifdef MFC_DEBUG
            if (proc_rank==0) print*, 'Computing axysimetric acoustic damping', mytime, hdid
#endif

            !$acc parallel loop collapse(3) gang vector default(present) reduction(+: sumIntensity_ac) &
            !$acc reduction(MAX: focalIntensity_ac, focalIntensity_ac_prms, abortFlag_max) &
            !$acc private(myalpha_rho, myalpha, vel_h, Re_h, rhoYks_h) &
            !$acc copy(sumIntensity_ac, focalIntensity_ac, focalIntensity_ac_prms, abortFlag_max)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        abortFlag = 0

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

                        if (f_is_default(absCoef)) then
                            print*, "HIFU: Check absCoef values!"
                            abortFlag = 1
                        end if

                        !>> Get the strain rate tensor (using central finite difference)
                        varA = 0._wp
                        varB = 0._wp

                        ! Only for axysimmetric assumption
                        duxdx = (q_prim_vf(contxe + 1)%sf(j + 1, k, 0) - q_prim_vf(contxe + 1)%sf(j - 1, k, 0))/(x_cc(j + 1) - x_cc(j - 1))
                        duxdr = (q_prim_vf(contxe + 1)%sf(j, k + 1, 0) - q_prim_vf(contxe + 1)%sf(j, k - 1, 0))/(y_cc(k + 1) - y_cc(k - 1))

                        durdx = (q_prim_vf(contxe + 2)%sf(j + 1, k, 0) - q_prim_vf(contxe + 2)%sf(j - 1, k, 0))/(x_cc(j + 1) - x_cc(j - 1))
                        durdr = (q_prim_vf(contxe + 2)%sf(j, k + 1, 0) - q_prim_vf(contxe + 2)%sf(j, k - 1, 0))/(y_cc(k + 1) - y_cc(k - 1))

                        !>> Get pressure, density and speed of sound
                        do i = 1, contxe
                            myalpha_rho(i) = q_prim_vf(i)%sf(j, k, l)
                            myalpha(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
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
                        q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) = max(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l), pres_h)
                        q_hifu%vf(hifu_params%P_idx + 1)%sf(j, k, l) = min(q_hifu%vf(hifu_params%P_idx + 1)%sf(j, k, l), pres_h)

                        !>> Compute intensity form acoustic damping

                        ! PRMS method (calculate only during the last time step in stage2 -> need developed Pmax field)
                        intensity_ac_prms = 0._wp
                        if (cfl_dt) then
                            if (mytime >= t_stop) then
                                intensity_ac_prms = absCoef*(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h)
                            end if
                        else
                            if (t_step == t_step_stop - 1) then
                                intensity_ac_prms = absCoef*(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h)
                            end if
                        end if

                        ! Shear stress method
                        intensity_ac = 0._wp
                        ep11 = durdr
                        ep22 = vel_h(2)/y_cc(k)
                        ep33 = duxdx
                        ep13 = 0.5_wp*(durdx + duxdr)
                        varA = ep11**2._wp + ep22**2._wp + ep33**2._wp
                        !varB = (8._wp/3._wp)*varA - (4._wp/3._wp)*(ep11*ep22 + ep11*ep33 + ep22*ep33) + 6._wp*(ep13**2._wp)
                        varB = (2._wp/3._wp)*(ep11**2._wp + ep22**2._wp + ep33**2._wp - ep11*ep22 - ep11*ep33 - ep22*ep33) + &
                                                                                                            2._wp*(ep13**2._wp)
                        intensity_ac = intensity_ac + bulkVisc*varA + 2._wp*shearVisc*varB !intensity is "q_us_ac"
                        !intensity_ac = intensity_ac + bulkVisc*(ep11 + ep22 + ep33)**2._wp + 2._wp*shearVisc*varB !intensity is "q_us_ac"

                        q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l) &
                                                                                                        + hdid      ! Update total sampling time
                        q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l) &
                                                                                            + intensity_ac*hdid     ! Sampling acoustic intensity
                        q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l) = intensity_ac_prms * &
                                                                       q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l)    ! Sampling acoustic intensity (prms)
                        
                        ! Checking for NaNs
                        if (q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l) /= q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l)) then
                            print*, 'Acoustic intensity is NaN', j, k, l, hdid, intensity_ac
                            print*, 'viscosities (bulk & shear)', bulkVisc, shearVisc
                            print*, 'var: A, B', varA, varB, ep11, ep22, ep33, ep13
                            print*, 'ep22:', vel_h(2), y_cc(k), rho_h
                            abortFlag = 1
                        end if

                        if (q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l) /= q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l)) then
                            print*, 'Acoustic intensity PRMS is NaN', j, k, l, hdid, intensity_ac_prms, &
                                                        absCoef,q_hifu%vf(hifu_params%P_idx)%sf(j, k, l), &
                                                        hifu_params%atmPres, rho_h, cson_h
                            abortFlag = 1
                        end if

                        abortFlag_max = max(abortFlag_max, abortFlag)

                        !Update average velocities for streaming
                        q_hifu%vf(hifu_params%u_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%u_idx)%sf(j, k, l) + vel_h(1)*hdid ! Sampling x-vel
                        q_hifu%vf(hifu_params%v_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%v_idx)%sf(j, k, l) + vel_h(2)*hdid ! Sampling y-vel

                        !Get focal intensity and velocities
                        axialCondition = (dy(k) > y_cc(k) .and. y_cc(k) > 0._wp)
                        radialCondition = (x_cb(j - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < x_cb(j))
                        condition = (axialCondition .and. radialCondition)
                        if (condition) then
                            print*, cson_h, pres_h, rho_h
                            print*, absCoef, shearVisc, bulkVisc
                            print*, 'q_ac', (absCoef*(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h))*q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l), q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l)
                            print*, 'visc terms focus', varA, varB
                            focalIntensity_ac = max(focalIntensity_ac, q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l))
                            focalIntensity_ac_prms = max(focalIntensity_ac_prms, q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l))
                        end if

                        !Intensity summation through the domain
                        sumIntensity_ac = sumIntensity_ac + q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l)

                    end do
                end do
            end do

            if (abortFlag_max > 0) stop "NaNs in Acoustic intensity (prms)"

            if (num_procs > 1) then
                tmp = sumIntensity_ac
                call s_mpi_allreduce_sum(tmp, sumIntensity_ac)

                tmp = focalIntensity_ac
                call s_mpi_allreduce_max(tmp, focalIntensity_ac)

                tmp = focalIntensity_ac_prms
                call s_mpi_allreduce_max(tmp, focalIntensity_ac_prms)
            end if

            !$acc update host(q_hifu%vf(hifu_params%tsamp_idx)%sf)

            if (proc_rank == 0) write (99, '(6x,5E24.8)') &
                                        mytime, &
                                        q_hifu%vf(hifu_params%tsamp_idx)%sf(0, 0, 0), &
                                        focalIntensity_ac, &
                                        focalIntensity_ac_prms, &
                                        sumIntensity_ac

        else if (.not. cyl_coord .and. p > 0) then !Cartesian 3D

#ifdef MFC_DEBUG
            if (proc_rank==0) print*, 'Computing cartesian 3D acoustic damping', mytime, hdid
#endif
            !$acc parallel loop collapse(3) gang vector default(present) reduction(+: sumIntensity_ac) &
            !$acc reduction(MAX: focalIntensity_ac, focalIntensity_ac_prms) &
            !$acc private(myalpha_rho, myalpha, vel_h, Re_h, rhoYks_h, duxdn, duydn, duzdn) &
            !$acc copy(sumIntensity_ac, focalIntensity_ac, focalIntensity_ac_prms)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        abortFlag = 0

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

                        if (f_is_default(absCoef)) then
                            abortFlag = 1
                            print*, "HIFU: Check absCoef values!"
                        end if

                        !>> Get the strain rate tensor (using central finite difference)
                        varA = 0._wp
                        varB = 0._wp
                        varC = 0._wp

                        mtd_idx = 1
                        call s_space_derivative(q_prim_vf(contxe + 1), j, k, l, duxdn, mtd_idx)
                        call s_space_derivative(q_prim_vf(contxe + 2), j, k, l, duydn, mtd_idx)
                        call s_space_derivative(q_prim_vf(contxe + 3), j, k, l, duzdn, mtd_idx)

                        !> First order centered difference approximation
                        ! duxdn(1) = (q_prim_vf(contxe + 1)%sf(j + 1, k, l) - q_prim_vf(contxe + 1)%sf(j - 1, k, l))/ (x_cc(j + 1) - x_cc(j - 1))
                        ! duxdn(2) = (q_prim_vf(contxe + 1)%sf(j, k + 1, l) - q_prim_vf(contxe + 1)%sf(j, k - 1, l))/ (y_cc(k + 1) - y_cc(k - 1))
                        ! duxdn(3) = (q_prim_vf(contxe + 1)%sf(j, k, l + 1) - q_prim_vf(contxe + 1)%sf(j, k, l - 1))/ (z_cc(l + 1) - z_cc(l - 1))

                        ! duydn(1) = (q_prim_vf(contxe + 2)%sf(j + 1, k, l) - q_prim_vf(contxe + 2)%sf(j - 1, k, l))/ (x_cc(j + 1) - x_cc(j - 1))
                        ! duydn(2) = (q_prim_vf(contxe + 2)%sf(j, k + 1, l) - q_prim_vf(contxe + 2)%sf(j, k - 1, l))/ (y_cc(k + 1) - y_cc(k - 1))
                        ! duydn(3) = (q_prim_vf(contxe + 2)%sf(j, k, l + 1) - q_prim_vf(contxe + 2)%sf(j, k, l - 1))/ (z_cc(l + 1) - z_cc(l - 1))

                        ! duzdn(1) = (q_prim_vf(contxe + 3)%sf(j + 1, k, l) - q_prim_vf(contxe + 3)%sf(j - 1, k, l))/ (x_cc(j + 1) - x_cc(j - 1))
                        ! duzdn(2) = (q_prim_vf(contxe + 3)%sf(j, k + 1, l) - q_prim_vf(contxe + 3)%sf(j, k - 1, l))/ (y_cc(k + 1) - y_cc(k - 1))
                        ! duzdn(3) = (q_prim_vf(contxe + 3)%sf(j, k, l + 1) - q_prim_vf(contxe + 3)%sf(j, k, l - 1))/ (z_cc(l + 1) - z_cc(l - 1))

                        !>> Get pressure, density and speed of sound
                        do i = 1, contxe
                            myalpha_rho(i) = q_prim_vf(i)%sf(j, k, l)
                            myalpha(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
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
                        q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) = max(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l), pres_h)
                        q_hifu%vf(hifu_params%P_idx + 1)%sf(j, k, l) = min(q_hifu%vf(hifu_params%P_idx + 1)%sf(j, k, l), pres_h)

                        !>> Compute intensity form acoustic damping

                        ! PRMS method (calculate only during the last time step in stage2 -> need developed Pmax field)
                        intensity_ac_prms = 0._wp
                        if (cfl_dt) then
                            if (mytime >= t_stop) then
                                intensity_ac_prms = absCoef*(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h)
                            end if
                        else
                            if (t_step == t_step_stop - 1) then
                                intensity_ac_prms = absCoef*(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h)
                            end if
                        end if

                        ! Shear stress method
                        ! Intensity is "q_us_ac"
                        intensity_ac = 0._wp
                        ep11 = duxdn(1)
                        ep22 = duydn(2)
                        ep33 = duzdn(3)
                        ep12 = 0.5_wp*(duxdn(2) + duydn(1))
                        ep13 = 0.5_wp*(duxdn(3) + duzdn(1))
                        ep23 = 0.5_wp*(duydn(3) + duzdn(2))
                        !varA = ep11**2._wp + ep22**2._wp + ep33**2._wp
                        !varB = (ep11 - varA/3._wp)**2._wp + (ep22 - varA/3._wp)**2._wp + (ep33 - varA/3._wp)**2._wp
                        !varB = varB + 2._wp*(ep12**2._wp + ep13**2._wp + ep23**2._wp)
                        !intensity_ac = intensity_ac + bulkVisc*varA + 2._wp*shearVisc*varB 

                        varA = ep11 + ep22 + ep33
                        varB = ep11**2._wp + ep22**2._wp + ep33**2._wp
                        !varB = varA**2._wp
                        varC = (2._wp/3._wp)*(ep11**2._wp + ep22**2._wp + ep33**2._wp - ep11*ep22 - ep11*ep33 - ep22*ep33) + &
                                                                                2._wp*(ep12**2._wp + ep13**2._wp + ep23**2._wp)
                        intensity_ac = intensity_ac + bulkVisc*varB + 2._wp*shearVisc*varC
                        !intensity_ac = 2._wp * intensity_ac
                        
                        ! varA = ep11 + ep22 + ep33
                        ! varB = 2._wp*(ep12**2._wp + ep13**2._wp + ep23**2._wp) + (ep11**2._wp + ep22**2._wp + ep33**2._wp)
                        ! varB = varB - (1._wp/3._wp)*(varA**2._wp)
                        ! intensity_ac = intensity_ac + bulkVisc*(varA**2._wp) + 2._wp*shearVisc*varB

                        q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l) &
                                                                                                        + hdid      ! Update total sampling time
                        q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l) &
                                                                                            + intensity_ac*hdid     ! Sampling acoustic intensity
                        q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l) = intensity_ac_prms * &
                                                                       q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l)    ! Sampling acoustic intensity (prms)
                        
                        ! Checking for NaNs
                        if (q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l) /= q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l)) then
                            print*, 'Acoustic intensity is NaN', j, k, l, hdid, intensity_ac
                            print*, 'viscosities (bulk & shear)', bulkVisc, shearVisc
                            print*, 'var: A, B', varA, varB, ep11, ep22, ep33, ep13
                            print*, 'ep22:', vel_h(2), y_cc(k), rho_h
                            abortFlag = 1
                        end if

                        if (q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l) /= q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l)) then
                            print*, 'Acoustic intensity PRMS is NaN', j, k, l, hdid, intensity_ac_prms, &
                                                        absCoef,q_hifu%vf(hifu_params%P_idx)%sf(j, k, l), &
                                                        hifu_params%atmPres, rho_h, cson_h
                            abortFlag = 1
                        end if

                        abortFlag_max = max(abortFlag_max, abortFlag)

                        !Update average velocities for streaming
                        q_hifu%vf(hifu_params%u_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%u_idx)%sf(j, k, l) + vel_h(1)*hdid ! Sampling x-vel
                        q_hifu%vf(hifu_params%v_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%v_idx)%sf(j, k, l) + vel_h(2)*hdid ! Sampling y-vel

                        !Get focal intensity and velocities
                        axialCondition = (dy(k) > abs(y_cc(k)) .and. abs(y_cc(k)) >= 0._wp)
                        if (p>0) axialCondition = axialCondition .and. (dz(l) > abs(z_cc(l)) .and. abs(z_cc(l)) >= 0._wp)
                        radialCondition = (x_cb(j - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < x_cb(j))
                        !if (p>0) radialCondition = (z_cb(l - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < z_cb(l))
                        condition = (axialCondition .and. radialCondition)
                        if (condition) then
                            ! print*, j, k, l, sys_size, hdid, dt, proc_rank
                            print*, cson_h, pres_h, rho_h
                            print*, absCoef, shearVisc, bulkVisc
                            print*, 'q_ac', (absCoef*(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h))*q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l), q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l)
                            print*, 'visc terms focus', varB, varC, (2._wp/3._wp)*(ep11**2._wp + ep22**2._wp + ep33**2._wp - ep11*ep22 - ep11*ep33 - ep22*ep33) + 2._wp*(ep12**2._wp + ep13**2._wp + ep23**2._wp)
                            ! print*, 'strain', ep11, ep22, ep33, ep13, ep12, ep23
                            ! print*, 'prim', q_prim_vf(1)%sf(j, k, l), q_prim_vf(2)%sf(j, k, l), q_prim_vf(3)%sf(j, k, l), q_prim_vf(4)%sf(j, k, l), q_prim_vf(5)%sf(j, k, l), q_prim_vf(6)%sf(j, k, l), q_prim_vf(7)%sf(j, k, l), q_prim_vf(8)%sf(j, k, l)
                            ! print*, 'cons', q_cons_vf(1)%sf(j, k, l), q_cons_vf(2)%sf(j, k, l), q_cons_vf(3)%sf(j, k, l), q_cons_vf(4)%sf(j, k, l), q_cons_vf(5)%sf(j, k, l), q_cons_vf(6)%sf(j, k, l), q_cons_vf(7)%sf(j, k, l), q_cons_vf(8)%sf(j, k, l)
                            focalIntensity_ac = max(focalIntensity_ac, q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l))
                            focalIntensity_ac_prms = max(focalIntensity_ac_prms, q_hifu%vf(hifu_params%qus_prms_idx)%sf(j, k, l))
                        end if

                        if (proc_rank==96 .and. j==5 .and. k==5 .and. l==5) then
                            ! print*, j, k, l, sys_size, hdid, dt, proc_rank
                            ! print*, cson_h, pres_h, rho_h
                            print*, 'q_ac', (absCoef*(q_hifu%vf(hifu_params%P_idx)%sf(j, k, l) - hifu_params%atmPres)**2._wp/(rho_h*cson_h))*q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l), q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l)
                            print*, 'visc terms', bulkVisc*varB, 2._wp*shearVisc*varC
                            ! print*, 'strain', ep11, ep22, ep33, ep13, ep12, ep23
                            !call s_mpi_abort('Debug sampler 3D')
                        end if

                        !Intensity summation through the domain
                        sumIntensity_ac = sumIntensity_ac + q_hifu%vf(hifu_params%qus_idx)%sf(j, k, l)

                    end do
                end do
            end do

            if (abortFlag_max > 0) stop "NaNs in Acoustic intensity"

            if (num_procs > 1) then
                tmp = sumIntensity_ac
                call s_mpi_allreduce_sum(tmp, sumIntensity_ac)

                tmp = focalIntensity_ac
                call s_mpi_allreduce_max(tmp, focalIntensity_ac)

                tmp = focalIntensity_ac_prms
                call s_mpi_allreduce_max(tmp, focalIntensity_ac_prms)
            end if

            !$acc update host(q_hifu%vf(hifu_params%tsamp_idx)%sf)

            if (proc_rank == 0) write (99, '(6x,5E24.8)') &
                                        mytime, &
                                        q_hifu%vf(hifu_params%tsamp_idx)%sf(0, 0, 0), &
                                        focalIntensity_ac, &
                                        focalIntensity_ac_prms, &
                                        sumIntensity_ac

        else
            call s_mpi_abort('Getting HIFU samples (stage 2) works only with axisymmetric assumption so far!')
        end if

    end subroutine s_update_HIFU_vars_sampling


    subroutine s_space_derivative(q_var, i, j, k, dumdn, mtd_idx)
        !$acc routine seq

        type(scalar_field), intent(in) :: q_var
        integer, intent(in) :: i, j, k, mtd_idx
        real(wp), dimension(3), intent(out) :: dumdn

        if (mtd_idx == 1) then 
            !> First order centered difference approximation
            dumdn(1) = (q_var%sf(i + 1, j, k) - q_var%sf(i - 1, j, k))/ &
                                                                (x_cc(i + 1) - x_cc(i - 1))
            dumdn(2) = (q_var%sf(i, j + 1, k) - q_var%sf(i, j - 1, k))/ &
                                                                (y_cc(j + 1) - y_cc(j - 1))
            if (p > 0) dumdn(3) = (q_var%sf(i, j, k + 1) - q_var%sf(i, j, k - 1))/ &
                                                                (z_cc(k + 1) - z_cc(k - 1))
        else if (mtd_idx == 2) then 
            !> Second order centered difference approximation
            dumdn(1) = q_var%sf(i, j, k)*(dx(i + 1) - dx(i - 1)) &
                                    + q_var%sf(i + 1, j, k)*(dx(i) + dx(i - 1)) &
                                    - q_var%sf(i - 1, j, k)*(dx(i) + dx(i + 1))
            dumdn(1) = dumdn(1) / ((dx(i) + dx(i - 1))*(dx(i) + dx(i + 1)))

            dumdn(2) = q_var%sf(i, j, k)*(dy(j + 1) - dy(j - 1)) &
                                    + q_var%sf(i, j + 1, k)*(dy(j) + dy(j - 1)) &
                                    - q_var%sf(i, j - 1, k)*(dy(j) + dy(j + 1))
            dumdn(2) = dumdn(2) / ((dy(j) + dy(j - 1))*(dy(j) + dy(j + 1)))
            if (p > 0) then 
                dumdn(3) = q_var%sf(i, j, k)*(dz(k + 1) - dz(k - 1)) &
                                        + q_var%sf(i, j, k + 1)*(dz(k) + dz(k - 1)) &
                                        - q_var%sf(i, j, k - 1)*(dz(k) + dz(k + 1))
                dumdn(3) = dumdn(3) / ((dz(k) + dz(k - 1))*(dz(k) + dz(k + 1)))
            end if
        end if

    end subroutine s_space_derivative

    !> The purpose of this procedure is to write the maximum and minimum pressure through time
        !!      along the axisymmetric and radial axes
        !! @param save_count File identifier
    subroutine s_write_Pmax(save_count)

        integer, intent(in) :: save_count

        integer :: j, k, l
        logical :: axialCondition, radialCondition, condition

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    ! Specify enough conditions for axial and radial probe lines
                    axialCondition = (dy(k) > abs(y_cc(k)) .and. abs(y_cc(k)) >= 0._wp)
                    if (p>0) axialCondition = axialCondition .and. (dz(l) > abs(z_cc(l)) .and. abs(z_cc(l)) >= 0._wp)
                    radialCondition = (x_cb(j - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < x_cb(j))
                    ! if (p > 0) radialCondition = (z_cb(l - 1) < acoustic_bc_params%focLen .and. acoustic_bc_params%focLen < z_cb(l))
                    if (p > 0) radialCondition = radialCondition .and. l==0
                    condition = (axialCondition .or. radialCondition)
                    if (condition) then
                        if (p>0) then
                            write (100, '(6x,6E24.8)') &
                                mytime, &
                                x_cc(j), &
                                y_cc(k), &
                                z_cc(l), &
                                q_hifu%vf(hifu_params%P_idx)%sf(j, k, l), &
                                q_hifu%vf(hifu_params%P_idx + 1)%sf(j, k, l)

                        else
                            write (100, '(6x,5E24.8)') &
                                mytime, &
                                x_cc(j), &
                                y_cc(k), &
                                q_hifu%vf(hifu_params%P_idx)%sf(j, k, l), &
                                q_hifu%vf(hifu_params%P_idx + 1)%sf(j, k, l)
                        end if
                        
                    end if
                end do
            end do
        end do

    end subroutine s_write_Pmax

    !Initilazile 3d domain to solve heat equation.
    subroutine s_initialize_from_2d_to_3d()

        integer :: i, j, k, l
        real(wp) :: dz_val, z_max

        ! Deallocate some old vars
        do i = 1, sys_size_HIFU
            MPI_IO_HIFU_DATA%var(i)%sf => null()
        end do
        
        !call s_finalize_mpi_proxy_module()

        if (hifu_params%cartesian) then
            !> 3D CARTESIAN COORDS
            
            if (proc_rank==0) print*, 'Initializing 3D Cartesian'
            
            !Run heat transfer solver in serial only, but can read from parallel stg2!
            m_hf = hifu_params%m
            n_hf = hifu_params%n
            p_hf = hifu_params%p
            num_dims = 3
            grid_geometry = 1
            cyl_coord = .false.
            sys_size_HIFU = hifu_params%qth_idx
            !$acc update device(m_hf, n_hf, p_hf, num_dims, grid_geometry, cyl_coord, sys_size_HIFU)

            !Deallocate variables to free memory
            if (bubbles_lagrange) call s_free_memory_stg3()

            !> all bc are assumed to be ghost cell extrapolation
            if (proc_rank==0) then
                bc_x%beg = -6; bc_x%end = -6
                bc_y%beg = -6; bc_y%end = -6
                bc_z%beg = -6; bc_z%end = -6
                !$acc update device(bc_x, bc_y, bc_z)
            end if

            !> Allocate variables
            !x
            @:ALLOCATE(x_cb_hf(-1 - buff_size:m_hf + buff_size))
            @:ALLOCATE(x_cc_hf(-buff_size:m_hf + buff_size))
            @:ALLOCATE(dx_hf(-buff_size:m_hf + buff_size))
            !y
            @:ALLOCATE(y_cb_hf(-1 - buff_size:n_hf + buff_size))
            @:ALLOCATE(y_cc_hf(-buff_size:n_hf + buff_size))
            @:ALLOCATE(dy_hf(-buff_size:n_hf + buff_size))
            !z
            @:ALLOCATE(z_cb_hf(-1 - buff_size:p_hf + buff_size))
            @:ALLOCATE(z_cc_hf(-buff_size:p_hf + buff_size))
            @:ALLOCATE(dz_hf(-buff_size:p_hf + buff_size))

            !> Generate cartesian mesh
            call s_generate_cartesian_mesh()

            !> New mpi vars
            sys_size_hyd = sys_size
            sys_size = 1
            do i = 1, sys_size_hifu
                allocate (MPI_IO_HIFU_DATA%var(i)%sf(0:m_hf, 0:n_hf, 0:p_hf))
                MPI_IO_HIFU_DATA%var(i)%sf => null()
            end do

            deallocate (start_idx)
            allocate (start_idx(1:num_dims))
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            !call s_initialize_mpi_proxy_module()

            !> Allocate 3d q_hifu
            @:ALLOCATE(q_hifu_3d%vf(1:sys_size_hifu))

            do i = 1, sys_size_hifu
                @:ALLOCATE(q_hifu_3d%vf(i)%sf(-buff_size:m_hf + buff_size, &
                    -buff_size:n_hf + buff_size, &
                    -buff_size:p_hf + buff_size))
            end do
            @:ACC_SETUP_VFs(q_hifu_3d)

            !> Populate grid (Interpolate values)
            call s_populate_cartesian_3D()

            !> Smear qvis and qth from the bubbles
            if (bubbles_lagrange) then 
                call s_smoothfunction(nBubs, intfc_rad, intfc_vel, &
                                                mtn_s, mtn_posPrev, q_hifu_3d, bub_qvis, bub_qth)
            end if

            !> Unify into a general domain
            call s_unify_cartesian_sources()
            call s_compare_interpolation()


        else
            !> 3D CYLINDRICAL COORDS

            ! Azimutal direction range
            !z_max = min(2._wp*pi, abs(hifu_params%z_max))
            z_max = 2._wp*pi

            p = hifu_params%p_cyl
            p_glb = p
            num_dims = 3

            if (bc_x%beg == -20) bc_x%beg = -6    ! from -20: acoustic bc
            bc_z%beg = -1; bc_z%end = -1          ! Assume entire cylindrical ring is taking care by one processor
            if (bc_y%beg == -2) bc_y%beg = -14    ! from -2: reflective boundary -> -21: cilyndrical sector (or)
                                                !                                 -14: full cylinder
            grid_geometry = 3

            bc_pole = bc_y%beg

            idwbuff(3)%beg = -buff_size
            idwbuff(3)%end = p - idwbuff(3)%beg
            !$acc update device(p, num_dims, bc_x, bc_y, bc_z, grid_geometry, idwbuff, bc_pole)

            !> Allocate variables
            !Theta axis
            @:ALLOCATE(z_cb(-1 - buff_size:p + buff_size))
            @:ALLOCATE(z_cc(-buff_size:p + buff_size))
            @:ALLOCATE(dz(-buff_size:p + buff_size))

            !New mpi vars
            sys_size_hyd = sys_size
            sys_size = 1
            do i = 1, sys_size
                allocate (MPI_IO_HIFU_DATA%var(i)%sf(0:m, 0:n, 0:p))
                MPI_IO_HIFU_DATA%var(i)%sf => null()
            end do
            !call s_initialize_mpi_proxy_module()

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
                            q_hifu_3d%vf(i)%sf(j, k, l) = q_hifu%vf(i)%sf(j, k, 0)
                            !if (i==1 .and. j==0 .and. k==0 .and. l==0) print*, q_hifu_3d%vf(i)%sf(j, k, l), q_hifu%vf(i)%sf(j, k, 0), hifu_params%Tref
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
        end if

    end subroutine s_initialize_from_2d_to_3d

    subroutine s_initialize_3d

    end subroutine

    !> Identifies the correct angle to rotate the bubble cloud to obtain the maximum, minimum or average temperature profiles.
    !>      Every hifu_params%z_max rads, it calculates the summation of the viscous dissipation of each bubble weighted with the distance
    !>      bewtwen the transducer's focal point and the bubble's location
    subroutine s_generate_cartesian_mesh()

        character(LEN=path_len + name_len) :: file_loc !<
        !! Generic string used to store the address of a file

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status

        integer :: i, j !< Generic loop integers

        real(wp) :: dx_tmp, dy_tmp, dz_tmp
        real(wp) :: dmin, dmax

        logical :: file_exist

        real(wp), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        allocate (x_cb_glb(-1:m_hf))
        allocate (y_cb_glb(-1:n_hf))
        allocate (z_cb_glb(-1:p_hf))

         ! Grid generation in the x-direction
        dx_tmp = abs(hifu_params%xe - hifu_params%xb)/real(m_hf + 1, wp)
        do i = 0, m_hf
            x_cb_hf(i - 1) = hifu_params%xb + dx_tmp*real(i, wp)
            x_cb_glb(i - 1) = x_cb_hf(i - 1)
        end do
        x_cb_hf(m_hf) = hifu_params%xe
        x_cb_glb(m_hf) = x_cb_hf(m_hf)
        dx_hf(0:m_hf) = x_cb_hf(0:m_hf) - x_cb_hf(-1:m_hf - 1)  ! Computing the cell width distribution
        x_cc_hf(0:m_hf) = x_cb_hf(-1:m_hf - 1) + dx_hf(0:m_hf)/2._wp ! Computing the cell center locations

         ! Grid generation in the y-direction
        dy_tmp = abs(hifu_params%ye + hifu_params%ye)/real(n_hf + 1, wp)
        do i = 0, n_hf
            y_cb_hf(i - 1) = -hifu_params%ye + dy_tmp*real(i, wp)
            y_cb_glb(i - 1) = y_cb_hf(i - 1)
        end do
        y_cb_hf(n_hf) = hifu_params%ye
        y_cb_glb(n_hf) = y_cb_hf(n_hf)
        dy_hf(0:n_hf) = y_cb_hf(0:n_hf) - y_cb_hf(-1:n_hf - 1)  ! Computing the cell width distribution
        y_cc_hf(0:n_hf) = y_cb_hf(-1:n_hf - 1) + dy_hf(0:n_hf)/2._wp ! Computing the cell center locations

         ! Grid generation in the z-direction
        dz_tmp = abs(hifu_params%ye + hifu_params%ye)/real(p_hf + 1, wp)
        do i = 0, p_hf
            z_cb_hf(i - 1) = -hifu_params%ye + dz_tmp*real(i, wp)
            z_cb_glb(i - 1) = z_cb_hf(i - 1)
        end do
        z_cb_hf(p_hf) = hifu_params%ye
        z_cb_glb(p_hf) = z_cb_hf(p_hf)
        dz_hf(0:p_hf) = z_cb_hf(0:p_hf) - z_cb_hf(-1:p_hf - 1)  ! Computing the cell width distribution
        z_cc_hf(0:p_hf) = z_cb_hf(-1:p_hf - 1) + dz_hf(0:p_hf)/2._wp ! Computing the cell center locations

        if (proc_rank == 0) then
            ! Write cell boundary locations to grid data files
            !x
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb_hf.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            data_size = m_hf + 2
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                            mpi_info_int, ifile, ierr)
            call MPI_FILE_WRITE(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
            !y
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb_hf.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            data_size = n_hf + 2
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                                mpi_info_int, ifile, ierr)
            call MPI_FILE_WRITE(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
            !z
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb_hf.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            data_size = p_hf + 2
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                                mpi_info_int, ifile, ierr)
            call MPI_FILE_WRITE(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        ! Populate grid buffers
        do i = 1, buff_size
            !x-beg
            dx_hf(-i) = dx_hf(0)
            x_cb_hf(-1 - i) = x_cb_hf(-i) - dx_hf(-i)
            x_cc_hf(-i) = x_cc_hf(1 - i) - (dx_hf(1 - i) + dx_hf(-i))/2._wp
            !x-end
            dx_hf(m_hf + i) = dx_hf(m_hf)
            x_cb_hf(m_hf + i) = x_cb_hf(m_hf + (i - 1)) + dx_hf(m_hf + i)
            x_cc_hf(m_hf + i) = x_cc_hf(m_hf + (i - 1)) + (dx_hf(m_hf + (i - 1)) + dx_hf(m_hf + i))/2._wp

            !y-beg
            dy_hf(-i) = dy_hf(0)
            y_cb_hf(-1 - i) = y_cb_hf(-i) - dy_hf(-i)
            y_cc_hf(-i) = y_cc_hf(1 - i) - (dy_hf(1 - i) + dy_hf(-i))/2._wp
            !y-end
            dy_hf(n_hf + i) = dy_hf(n_hf)
            y_cb_hf(n_hf + i) = y_cb_hf(n_hf + (i - 1)) + dy_hf(n_hf + i)
            y_cc_hf(n_hf + i) = y_cc_hf(n_hf + (i - 1)) + (dy_hf(n_hf + (i - 1)) + dy_hf(n_hf + i))/2._wp

            !z-beg
            dz_hf(-i) = dz_hf(0)
            z_cb_hf(-1 - i) = z_cb_hf(-i) - dz_hf(-i)
            z_cc_hf(-i) = z_cc_hf(1 - i) - (dz_hf(1 - i) + dz_hf(-i))/2._wp
            !z-end
            dz_hf(p_hf + i) = dz_hf(p_hf)
            z_cb_hf(p_hf + i) = z_cb_hf(p_hf + (i - 1)) + dz_hf(p_hf + i)
            z_cc_hf(p_hf + i) = z_cc_hf(p_hf + (i - 1)) + (dz_hf(p_hf + (i - 1)) + dz_hf(p_hf + i))/2._wp

        end do

        ! Report max and min dx, dy, dz
        !x
        dmin = abs(dflt_real)
        dmax = -abs(dflt_real)
        do i = 0, m_hf
            dmin = min(dmin, dx_hf(i))
            dmax = max(dmax, dx_hf(i))
        end do
        if (proc_rank==0) print*, 'New mesh: x-dir:', dmin, dmax, x_cb_hf(-1), x_cb_hf(m_hf), m_hf
        !y
        dmin = abs(dflt_real)
        dmax = -abs(dflt_real)
        do i = 0, n_hf
            dmin = min(dmin, dy_hf(i))
            dmax = max(dmax, dy_hf(i))
        end do
        if (proc_rank==0) print*, 'New mesh: y-dir:', dmin, dmax, y_cb_hf(-1), y_cb_hf(n_hf), n_hf
        !z
        dmin = abs(dflt_real)
        dmax = -abs(dflt_real)
        do i = 0, p_hf
            dmin = min(dmin, dz_hf(i))
            dmax = max(dmax, dz_hf(i))
        end do
        if (proc_rank==0) print*, 'New mesh: z-dir:', dmin, dmax, z_cb_hf(-1), z_cb_hf(p_hf), p_hf

        !$acc update device(x_cb_hf, y_cb_hf, z_cb_hf, x_cc_hf, y_cc_hf, z_cc_hf, dx_hf, dy_hf, dz_hf)

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

    end subroutine s_generate_cartesian_mesh

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

    subroutine s_populate_cartesian_3D()

        integer :: j, k, l
        integer :: cellx, celly, cellz

        ! 3d q_hifu: Temperature
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = -buff_size, p_hf + buff_size
            do k = -buff_size, n_hf + buff_size
                do j = -buff_size, m_hf + buff_size
                    q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) = 1._wp
                    q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l) = q_hifu%vf(hifu_params%tsamp_idx)%sf(0, 0, 0)
                    q_hifu_3d%vf(hifu_params%qus_idx)%sf(j, k, l) = 0._wp
                    q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l) = 0._wp
                    q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l) = 0._wp
                end do
            end do
        end do

        ! 3d q_hifu: Acoustic intensity
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p_hf 
            do k = 0, n_hf
                do j = 0, m_hf
                    cellx = j; celly = k; cellz = l
                    q_hifu_3d%vf(hifu_params%qus_idx)%sf(j, k, l) = f_interpolate_qus(cellx, celly, cellz)
                end do
            end do
        end do

        if (proc_rank==0) print*, 'Grid populated: Initial temp & qus'


    end subroutine s_populate_cartesian_3D

    function f_interpolate_qus(j, k, l)
#ifdef _CRAYFTN
    !DIR$ INLINEALWAYS f_interpolate_qus
#else
    !$acc routine seq
#endif
        integer, intent(in) :: j, k, l
        real(wp) :: f_interpolate_qus, r_cc, minDist, minDist_old
        real(wp) :: nSamples, valCloseCell, valAvg
        real(wp) :: r_cb_1, r_cb_2, tmp1, tmp2, tmp3, tmp4, tmp
        integer :: i, q
        integer :: cell2D_x, cell2D_r
        integer :: cell2D_xb, cell2D_xe
        integer :: cell2D_rb, cell2D_re

        r_cc = sqrt(y_cc_hf(k)**2._wp + z_cc_hf(l)**2._wp)

        tmp1 = sqrt(y_cb_hf(k - 1)**2._wp + z_cb_hf(l - 1)**2._wp)
        tmp2 = sqrt(y_cb_hf(k)**2._wp + z_cb_hf(l)**2._wp)
        tmp3 = sqrt(y_cb_hf(k - 1)**2._wp + z_cb_hf(l)**2._wp)
        tmp4 = sqrt(y_cb_hf(k)**2._wp + z_cb_hf(l - 1)**2._wp)
        r_cb_1 = min(min(min(tmp1, tmp2), tmp3), tmp4)
        r_cb_2 = max(max(max(tmp1, tmp2), tmp3), tmp4)
        if ((y_cb_hf(k)*y_cb_hf(k - 1) <= 0._wp) .and. (z_cb_hf(l)*z_cb_hf(l - 1) <= 0._wp)) r_cb_1 = 0._wp
        !if (abs(r_cb_1 - r_cb_2) < 0.5_wp*dz_hf(l)) r_cb_1 = 0._wp

        !x-dir-beg
        cell2D_xb = 0
        do while(.true.)
            if( x_cc(cell2D_xb) >= x_cb_hf(j-1)) exit
            cell2D_xb = cell2D_xb + 1
            if (cell2D_xb > m + buff_size) return 
        end do
        !x-dir-end
        cell2D_xe = -buff_size
        do while(.true.)
            if( x_cc(cell2D_xe) >= x_cb_hf(j)) exit
            cell2D_xe = cell2D_xe + 1
            if (cell2D_xe > m + buff_size) return
        end do
        cell2D_xe = cell2D_xe - 1

        !r-dir-beg
        cell2D_rb = 0
        do while(.true.)
            if( y_cc(cell2D_rb) >= r_cb_1) exit
            cell2D_rb = cell2D_rb + 1
            if (cell2D_rb > n + buff_size) return
        end do
        !r-dir-end
        cell2D_re = 0
        do while(.true.)
            if( y_cc(cell2D_re) >= r_cb_2) exit
            cell2D_re = cell2D_re + 1
            if (cell2D_re > n + buff_size) return
        end do
        cell2D_re = cell2D_re - 1

        nSamples = 0._wp
        valAvg = 0._wp
        !$acc loop seq
        do i = cell2D_xb, cell2D_xe
            !$acc loop seq
            do q = cell2D_rb, cell2D_re
                valAvg = valAvg + q_hifu%vf(hifu_params%qus_idx)%sf(i, q, 0)
                nSamples = nSamples + 1._wp
            end do
        end do

        !if (nSamples >= 100._wp) print*, nSamples, cell2D_xb, cell2D_xe, cell2D_rb, cell2D_re

        if (nSamples <= 0.0_wp) then

            !x-dir
            cell2D_x = 0
            do while(.true.)
                if (x_cc_hf(j) >= x_cb(cell2D_x-1) .and. x_cc_hf(j) < x_cb(cell2D_x)) exit
                cell2D_x = cell2D_x + 1
                if (cell2D_x > m + buff_size) return
            end do

            !r-dir
            cell2D_r = 0
            do while(.true.)
                if (r_cc >= y_cb(cell2D_r-1) .and. r_cc < y_cb(cell2D_r)) exit
                cell2D_r = cell2D_r + 1
                if (cell2D_r > n + buff_size) return
            end do

            f_interpolate_qus = q_hifu%vf(hifu_params%qus_idx)%sf(cell2D_x, cell2D_r, 0)

        else

            f_interpolate_qus = valAvg/nSamples

        end if

    end function f_interpolate_qus

    
    subroutine s_unify_cartesian_sources()

        integer :: i, j, k, l
        real(wp) :: tmp_local, tmp_global
        
        if (num_procs == 1) return

        do i = 1, sys_size_hifu
            !$acc update host(q_hifu_3d%vf(i)%sf)
        end do

        ! Done by the CPUs only

        do l = 0, p_hf
            do k = 0, n_hf
                do j = 0, m_hf
                    !qus
                    tmp_local = q_hifu_3d%vf(hifu_params%qus_idx)%sf(j, k, l)
                    call s_mpi_allreduce_max(tmp_local, tmp_global)
                    q_hifu_3d%vf(hifu_params%qus_idx)%sf(j, k, l) = tmp_global

                    !qvis
                    tmp_local = q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l)
                    call s_mpi_allreduce_sum(tmp_local, tmp_global)
                    q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l) = tmp_global

                    !qth
                    tmp_local = q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l)
                    call s_mpi_allreduce_sum(tmp_local, tmp_global)
                    q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l) = tmp_global

                end do
            end do
        end do

        do i = 1, sys_size_hifu
            !$acc update device(q_hifu_3d%vf(i)%sf)
        end do

         if (proc_rank==0) print*, 'Sources got unified'


    end subroutine s_unify_cartesian_sources

    subroutine s_unify_temperature_field()

        integer :: i, j, k, l
        real(wp) :: tmp_local, tmp_global
        
        if (num_procs == 1) return

        !$acc update host(q_hifu_3d%vf(hifu_params%T_idx)%sf)

        ! Done by the CPUs only

        do l = 0, p_hf
            do k = 0, n_hf
                do j = 0, m_hf

                    tmp_local = 0._wp
                    if (proc_rank==0) tmp_local = q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l)
                    call s_mpi_allreduce_max(tmp_local, tmp_global)
                    q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) = tmp_global

                end do
            end do
        end do

         !$acc update device(q_hifu_3d%vf(hifu_params%T_idx)%sf)

         if (proc_rank==0) print*, 'Temperature field got unified'

    end subroutine s_unify_temperature_field

    subroutine s_compare_interpolation()

        real(wp) :: max_old, max_new
        real(wp) :: min_old, min_new
        real(wp) :: max_qvis
        real(wp) :: max_qvis_smooth
        real(wp) :: tmp_local, tmp_global, sampledTime
        integer :: j, k, l, i
        character(LEN=path_len + 2*name_len) :: file_loc

        ! 2D field acoustic intensity
        max_old = -abs(dflt_real)
        min_old = abs(dflt_real)

        !$acc parallel loop collapse(2) gang vector default(present) reduction(MAX: max_old) reduction(MIN: min_old) copy(max_old, min_old)
        do k = 0, n
            do j = 0, m
                max_old = max(max_old, q_hifu%vf(hifu_params%qus_idx)%sf(j, k, 0)/q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, 0))
                min_old = min(min_old, q_hifu%vf(hifu_params%qus_idx)%sf(j, k, 0)/q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, 0))
            end do
        end do

        if (num_procs>1) then
            tmp_local = max_old
            call s_mpi_allreduce_max(tmp_local, tmp_global)
            max_old = tmp_global

            tmp_local = min_old
            call s_mpi_allreduce_min(tmp_local, tmp_global)
            min_old = tmp_global
        end if
        
        if (proc_rank == 0) then

            ! 3D interpolated field acoustic intensity and smoothened viscous damping
            max_new = -abs(dflt_real)
            min_new = abs(dflt_real)
            sampledTime = 0._wp
            max_qvis_smooth = -abs(dflt_real)

            !$acc parallel loop collapse(3) gang vector default(present) reduction(MAX: max_new, sampledTime, max_qvis_smooth) &
            !$acc reduction(MIN: min_new) copy(max_new, min_new, sampledTime, max_qvis_smooth)
            do l = 0, p_hf
                do k = 0, n_hf
                    do j = 0, m_hf
                        max_new = max(max_new, q_hifu_3d%vf(hifu_params%qus_idx)%sf(j, k, l)/q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l))
                        min_new = min(min_new, q_hifu_3d%vf(hifu_params%qus_idx)%sf(j, k, l)/q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l))
                        sampledTime = max(sampledTime, q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l))
                        max_qvis_smooth = max(max_qvis_smooth, q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l)/q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l))
                    end do
                end do
            end do

            ! Report stats
            print*, 'Min avg qus 2D:', min_old
            print*, 'Min avg qus 3D:', min_new
            print*, 'Max avg qus 2D:', max_old
            print*, 'Max avg qus 3D:', max_new

        end if

        ! Printing viscous and thermal intensities in W for all bubbles in a separate filE
        max_qvis = -abs(dflt_real)

        !$acc parallel loop gang vector default(present) reduction(MAX: max_qvis) copy(max_qvis)
        do i = 1, nBubs
            max_qvis = max(max_qvis, bub_qvis(i)/q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(0,0,0))
        end do

        if (num_procs>1) then
            tmp_local = max_qvis
            call s_mpi_allreduce_max(tmp_local, tmp_global)
            max_qvis = tmp_global

            tmp_local = max_qvis_smooth
            call s_mpi_allreduce_max(tmp_local, tmp_global)
            max_qvis_smooth = tmp_global
        end if

        ! Report general stats
        if (proc_rank==0) then
            print*, 'Max avg q_vis (lagrange):', max_qvis
            print*, 'Max avg q_vis (euler):', max_qvis_smooth
            print*, 'Sampled time:', sampledTime
            print*, 'Host viscosity:', mul0
        end if

        if (proc_rank == 0) print*, 'Printing avg viscous and thermal intensities in W for all bubbles in ./D/ file'

        write (file_loc, '(A,I0,A)') 'bubble_intensities_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        open (13, FILE=trim(file_loc), FORM='formatted', position='rewind')
        write (13, *) 'bub_id, xPos, rPos, thetaPos, R0, avg_qvis, avg_qth'

        do k = 1, nBubs
            write (13, '(6X,I24.8,6e24.8)') &
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                bub_R0(k), &
                bub_qvis(k)/sampledTime, &
                bub_qth(k)/sampledTime
        end do

        close (13)

    end subroutine s_compare_interpolation

    !Finalize 3d cylindrical domain to solve heat equation.
    subroutine s_restore_initial_setup()

        integer :: i

        if ((.not. hifu_params%stg1) .and. (.not. hifu_params%stg2)) return

        if (hifu_params%stg3_3d) then !From 2D to 3D

            if (hifu_params%cartesian) then

                !Deallocate vars
                @:DEALLOCATE(x_cb_hf, x_cc_hf, dx_hf)
                @:DEALLOCATE(y_cb_hf, y_cc_hf, dy_hf)
                @:DEALLOCATE(z_cb_hf, z_cc_hf, dz_hf)

                do i = 1, sys_size_hifu
                    @:DEALLOCATE(q_hifu_3d%vf(i)%sf)
                end do
                @:DEALLOCATE(q_hifu_3d%vf)

                ! Restore 2D params
                num_dims = 2
                grid_geometry = 2
                !$acc update device(num_dims, grid_geometry)
                if (proc_rank==0) then
                    bc_x%beg = -20; bc_x%end = 10
                    bc_y%beg = -2; bc_y%end = 10
                    bc_z%beg = dflt_int; bc_z%end = dflt_int
                    !$acc update device(bc_x, bc_y, bc_z)
                end if

            else

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
            end if

        else !full 2D or 3D

            if (bc_x%beg == -6) bc_x%beg = -20

        end if

    end subroutine s_restore_initial_setup

    subroutine s_initialize_pure_3D()

        call s_print_hifu_source_stats(hifu_params%qus_idx)

        if (bubbles_lagrange) then
            if (proc_rank == 0) print*, 'Adding bubbles in pure 3D domain'
            call s_smoothfunction(nBubs, intfc_rad, intfc_vel, &
                                    mtn_s, mtn_posPrev, q_hifu, bub_qvis, bub_qth)
            call s_print_hifu_source_stats(hifu_params%qvis_idx)
            call s_print_hifu_source_stats(hifu_params%qth_idx)
        end if

    end subroutine s_initialize_pure_3D

    subroutine s_print_hifu_source_stats(idx)

        integer, intent(in) :: idx

        real(wp) :: max_val, min_val, tmp_local, tmp_global, val_test
        integer :: j, k, l

        max_val = -abs(dflt_real)
        min_val = abs(dflt_real)

        !$acc parallel loop collapse(3) gang vector default(present) reduction(MAX: max_val) &
        !$acc reduction(MIN: min_val) copy(max_val, min_val) copyin(idx)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    val_test = q_hifu%vf(idx)%sf(j, k, l)/&
                                           q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l)
                    max_val = max(max_val, val_test)
                    min_val = min(min_val, val_test)
                end do
            end do
        end do

        if (num_procs>1) then
            tmp_local = max_val
            call s_mpi_allreduce_max(tmp_local, tmp_global)
            max_val = tmp_global

            tmp_local = min_val
            call s_mpi_allreduce_min(tmp_local, tmp_global)
            min_val = tmp_global
        end if

        if (proc_rank == 0) then
            if (idx == hifu_params%qus_idx) print*, 'q_us (min, max):', min_val, max_val
            if (idx == hifu_params%qvis_idx) print*, 'q_vis smeared (min, max):', min_val, max_val
            if (idx == hifu_params%qth_idx) print*, 'q_th smeared (min, max):', min_val, max_val
        end if

    end subroutine s_print_hifu_source_stats

    !Calculate the rhs value from heat transfer eqn discretized with finite volumes.
    subroutine s_rhs_heatEqn(q_cons_vf, pb, mv, t_step, bc_type)

        type(scalar_field), dimension(sys_size_hyd), intent(in) :: q_cons_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: t_step
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        real(wp) :: CFL_heat, val_tmp, CFL_heat_old, CFL_heat_2, CFL_heat_3, CFL_heat_4, CFL_heat_max
        real(wp) :: dTdx_L, dTdx_R, dTdr_L, dTdr_R
        real(wp) :: dTdx, dTdr, dTdz
        real(wp) :: Tx_L, Tx_R, Tr_L, Tr_R
        real(wp) :: dTdz_L, dTdz_R, Tz_L, Tz_R
        real(wp) :: Ux_L, Ux_R, Ur_L, Ur_R
        real(wp) :: rho_cp, tdiff, alpha
        integer :: i, j, k, l, qus_hifu_idx_ht
        real(wp) :: abortFlag, abortFlag_max

        if (hifu_params%intPrms) then
            qus_hifu_idx_ht = hifu_params%qus_prms_idx
        else
            qus_hifu_idx_ht = hifu_params%qus_idx
        end if

        CFL_heat_max = -100_wp
        abortFlag_max = 0._wp

        if (hifu_params%cartesian) then !< From axisymmetric to 3D cartesian

            if (proc_rank == 0) then

                ! call s_populate_variables_buffers(q_hifu_3d%vf, pb, mv, bc_type)
                ! Assume boundaries are far away from the heating and that do not undergo any heating up.
                ! Buffers are equal to Tref as set in the initial condition

                !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht, t_step) &
                !$acc reduction(MAX: abortFlag_max, CFL_heat_max) copy(abortFlag_max, CFL_heat_max)
                do l = 0, p_hf
                    do k = 0, n_hf
                        do j = 0, m_hf

                            abortFlag = 0._wp
                            CFL_heat = -100._wp

                            !<  Zeroing RHS_heat
                            q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) = 0._wp

                            !> Temperature derivatives at the cell center. METHOD: Second order centered difference approximation
                            dTdx = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j + 1, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j - 1, k, l))/(x_cc_hf(j + 1) - x_cc_hf(j - 1))
                            dTdx_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j - 2, k, l))/(x_cc_hf(j) - x_cc_hf(j - 2))
                            dTdx_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j + 2, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(x_cc_hf(j + 2) - x_cc_hf(j))

                            dTdr = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k + 1, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k - 1, l))/(y_cc_hf(k + 1) - y_cc_hf(k - 1))
                            dTdr_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k - 2, l))/(y_cc_hf(k) - y_cc_hf(k - 2))
                            dTdr_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k + 2, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(y_cc_hf(k + 2) - y_cc_hf(k))

                            dTdz = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + 1) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l - 1))/(z_cc_hf(l + 1) - z_cc_hf(l - 1))
                            dTdz_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l - 2))/(z_cc_hf(l) - z_cc_hf(l - 2))
                            dTdz_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + 2) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(z_cc_hf(l + 2) - z_cc_hf(l))

                            !> Find temperature derivatives at the faces of the cell
                            dTdx_L = (dTdx*(x_cc_hf(j) - x_cb_hf(j - 1)) + dTdx_L*(x_cb_hf(j-1) - x_cc_hf(j-1)))/(x_cc_hf(j) - x_cc_hf(j-1))
                            dTdr_L = (dTdr*(y_cc_hf(k) - y_cb_hf(k - 1)) + dTdr_L*(y_cb_hf(k-1) - y_cc_hf(k-1)))/(y_cc_hf(k) - y_cc_hf(k-1))
                            dTdz_L = (dTdz*(z_cc_hf(l) - z_cb_hf(l - 1)) + dTdz_L*(z_cb_hf(l-1) - z_cc_hf(l-1)))/(z_cc_hf(l) - z_cc_hf(l-1))

                            dTdx_R = (dTdx*(x_cb_hf(j) - x_cc_hf(j)) + dTdx_R*(x_cc_hf(j+1) - x_cb_hf(j)))/(x_cc_hf(j+1) - x_cc_hf(j))
                            dTdr_R = (dTdr*(y_cb_hf(k) - y_cc_hf(k)) + dTdr_R*(y_cc_hf(k+1) - y_cb_hf(k)))/(y_cc_hf(k+1) - y_cc_hf(k))
                            dTdz_R = (dTdz*(z_cb_hf(l) - z_cc_hf(l)) + dTdz_R*(z_cc_hf(l+1) - z_cb_hf(l)))/(z_cc_hf(l+1) - z_cc_hf(l))

                            !> Get thermal properties (Assume host is num_fluids-1)
                            rho_cp = rho_cp_fluids(num_fluids - 1)
                            tdiff = tdiff_fluids(num_fluids - 1)

                            if (f_is_default(rho_cp) .or. f_is_default(tdiff)) then
                                print *, 'alpha, rho_cp, tdiff', alpha, rho_cp, tdiff
                                print*, "HIFU: Check thermal properties!"
                                abortFlag = 1._wp
                            end if

                            !> Obtain rhs (see notes)
                            q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) + tdiff*( &
                                                                            (1._wp/dx_hf(j))*(dTdx_R - dTdx_L) + &
                                                                            (1._wp/dy_hf(k))*(dTdr_R - dTdr_L) + &
                                                                            (1._wp/dz_hf(l))*(dTdz_R - dTdz_L))


                            if ((q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l) > 0._wp) .and. (t_step < hifu_params%stepStopSource)) then
                                !> Adding the heat source terms
                                q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) + &
                                                                                (1._wp/(rho_cp*q_hifu_3d%vf(hifu_params%tsamp_idx)%sf(j, k, l)))*( &
                                                                                q_hifu_3d%vf(qus_hifu_idx_ht)%sf(j, k, l) + &         !Acoustic intensity
                                                                                q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l) + &    !Viscous intensity
                                                                                q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l))        !Thermal intensity


                                if (hifu_params%streaming) then
                                    !> Convected heat flux
                                    print*, "HIFU: No streaming valid for 3D heat solver!"
                                    abortFlag = 1._wp
                                end if
                            end if

                            !Checking NaNs
                            if (q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) /= q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l)) then
                                print*, 'NaNs in q hifu rhs', q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l), j, k, l
                                print*, 'dx, dy, dz', dx_hf(j), dy_hf(k), dz_hf(l)
                                print*, 'x, y, z', x_cc_hf(j), y_cc_hf(k), z_cc_hf(l)
                                print*, 'Heat sources: ', q_hifu_3d%vf(qus_hifu_idx_ht)%sf(j, k, l), q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l), &
                                                            q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l)
                                print*, "Reduce dt!!"
                                abortFlag = 1._wp
                            end if

                            ! Calculate min CFL
                            if (t_step == 0) then
                                CFL_heat = max(CFL_heat, tdiff*dt/(dx_hf(j)**2_wp))
                                CFL_heat = max(CFL_heat, tdiff*dt/(dy_hf(k)**2_wp))
                                CFL_heat = max(CFL_heat, tdiff*dt/(dz_hf(l)**2_wp))
                                CFL_heat_max = max(CFL_heat_max, CFL_heat)
                            end if

                            abortFlag_max = max(abortFlag_max, abortFlag)

                        end do
                    end do
                end do
                
            end if

        else

            if (p == 0 .and. cyl_coord) then !< Axisymmetric rhs

                ! call s_populate_HIFU_variables_buffers(q_hifu)
                call s_populate_variables_buffers(q_hifu%vf, pb, mv, bc_type)

                !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht, t_step) &
                !$acc reduction(MAX: abortFlag_max, CFL_heat_max) copy(abortFlag_max, CFL_heat_max)
                do l = 0, p
                    do j = 0, m
                        do k = 0, n

                            abortFlag = 0._wp
                            CFL_heat = -100._wp

                            !<  Zeroing RHS_heat
                            q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) = 0._wp

                            !> Find temperature derivatives at the faces of the cell
                            dTdx_L = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j - 1, k, l))/(x_cc(j) - x_cc(j - 1))
                            dTdx_R = (q_hifu%vf(hifu_params%T_idx)%sf(j + 1, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k, l))/(x_cc(j + 1) - x_cc(j))
                            dTdr_L = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k - 1, l))/(y_cc(k) - y_cc(k - 1))
                            dTdr_R = (q_hifu%vf(hifu_params%T_idx)%sf(j, k + 1, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k, l))/(y_cc(k + 1) - y_cc(k))

                            !> Find temperature and streaming velocities at the faces of the cell
                            Tx_L = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%T_idx)%sf(j - 1, k, l))/2._wp
                            Ux_L = (q_hifu%vf(hifu_params%u_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%u_idx)%sf(j - 1, k, l))/2._wp
                            Tx_R = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%T_idx)%sf(j + 1, k, l))/2._wp
                            Ux_R = (q_hifu%vf(hifu_params%u_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%u_idx)%sf(j + 1, k, l))/2._wp
                            Tr_L = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%T_idx)%sf(j, k - 1, l))/2._wp
                            Ur_L = (q_hifu%vf(hifu_params%v_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%v_idx)%sf(j, k - 1, l))/2._wp
                            Tr_R = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%T_idx)%sf(j, k + 1, l))/2._wp
                            Ur_R = (q_hifu%vf(hifu_params%v_idx)%sf(j, k, l) + q_hifu%vf(hifu_params%v_idx)%sf(j, k + 1, l))/2._wp

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
                                print*, "HIFU: Check thermal properties!"
                                abortFlag = 1._wp
                            end if

                            !> Obtain rhs (see notes)
                            q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) + &
                                                                        tdiff*(1._wp/dx(j))*(dTdx_R - dTdx_L) + &
                                                                        tdiff*(1._wp/(2._wp*y_cc(k)*dy(k)))*((2._wp*y_cc(k) + dy(k))*dTdr_R - (2._wp*y_cc(k) - dy(k))*dTdr_L)

                            if ((q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l) > 0._wp) .and. (t_step < hifu_params%stepStopSource)) then
                                !> Adding the heat source terms
                                q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) + &
                                                                            (1._wp/(rho_cp))*(1._wp/q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l))*q_hifu%vf(qus_hifu_idx_ht)%sf(j, k, l) + &  !Acoustic intensity
                                                                            (1._wp/(rho_cp))*(1._wp/q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l))*q_hifu%vf(hifu_params%qvis_idx)%sf(j, k, l) + &    !Viscous intensity
                                                                            (1._wp/(rho_cp))*(1._wp/q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l))*q_hifu%vf(hifu_params%qth_idx)%sf(j, k, l)         !Thermal intensity

                                if (hifu_params%streaming) then
                                    !> Convected heat flux
                                    q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) = q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) - &
                                                                                (1._wp/q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l))*( & !Double check this
                                                                                (1._wp/dx(j))*(Ux_R*Tx_R - Ux_L*Tx_L) + &
                                                                                (1._wp/(2._wp*y_cc(k)*dy(k)))*((2._wp*y_cc(k) + dy(k))*Ur_R*Tr_R - (2._wp*y_cc(k) - dy(k))*Ur_L*Tr_L))
                                end if
                            end if

                            !Checking NaNs
                            if (q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) /= q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l)) then
                                print*, 'NaNs in q hifu rhs', q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l), j, k, l
                                print*, 'Current courant number', dx(j), dy(k), dz(l)
                                print*, 'Therm. properties', tdiff, rho_cp, dt
                                print*, "NaNs in q hifu rhs"
                                abortFlag = 1._wp
                            end if

                            ! Calculate max CFL
                            if (t_step == 0) then
                                CFL_heat = max(CFL_heat, tdiff*dt/(dx(j)**2_wp))
                                CFL_heat = max(CFL_heat, tdiff*dt/(dy(k)**2_wp))
                                CFL_heat_max = max(CFL_heat_max, CFL_heat)
                            end if

                            abortFlag_max = max(abortFlag_max, abortFlag)

                        end do
                    end do
                end do

                if (proc_rank==0 .and. t_step == 0) print*, 'Max CFL:', CFL_heat

            else 

                if (cyl_coord) then !< from axisymmetric to 3D Cylindrical

                    call s_populate_variables_buffers(q_hifu_3d%vf, pb, mv, bc_type)

                    !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht, t_step) &
                    !$acc reduction(MAX: abortFlag_max, CFL_heat_max) copy(abortFlag_max, CFL_heat_max)
                    do l = 0, p
                        do j = hifu_params%mb, hifu_params%me
                            do k = 0, hifu_params%ne

                                abortFlag = 0._wp
                                CFL_heat = -100._wp

                                !<  Zeroing RHS_heat
                                q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l) = 0._wp

                                !> Temperature derivatives at the cell center. METHOD: Second order centered difference approximation
                                dTdx = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j + 1, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j - 1, k, l))/(x_cc(j + 1) - x_cc(j - 1))
                                dTdx_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j - 2, k, l))/(x_cc(j) - x_cc(j - 2))
                                dTdx_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j + 2, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(x_cc(j + 2) - x_cc(j))

                                dTdr = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k + 1, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k - 1, l))/(y_cc(k + 1) - y_cc(k - 1))
                                dTdr_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k - 2, l))/(y_cc(k) - y_cc(k - 2))
                                dTdr_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k + 2, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(y_cc(k + 2) - y_cc(k))

                                dTdz = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + 1) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                                dTdz_L = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l - 2))/(z_cc(l) - z_cc(l - 2))
                                dTdz_R = (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + 2) - q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l))/(z_cc(l + 2) - z_cc(l))

                                !> Find temperature derivatives at the faces of the cell
                                dTdx_L = (dTdx*(x_cc(j) - x_cb(j - 1)) + dTdx_L*(x_cb(j-1) - x_cc(j-1)))/(x_cc(j) - x_cc(j-1))
                                dTdr_L = (dTdr*(y_cc(k) - y_cb(k - 1)) + dTdr_L*(y_cb(k-1) - y_cc(k-1)))/(y_cc(k) - y_cc(k-1))
                                dTdz_L = (dTdz*(z_cc(l) - z_cb(l - 1)) + dTdz_L*(z_cb(l-1) - z_cc(l-1)))/(z_cc(l) - z_cc(l-1))

                                dTdx_R = (dTdx*(x_cb(j) - x_cc(j)) + dTdx_R*(x_cc(j+1) - x_cb(j)))/(x_cc(j+1) - x_cc(j))
                                dTdr_R = (dTdr*(y_cb(k) - y_cc(k)) + dTdr_R*(y_cc(k+1) - y_cb(k)))/(y_cc(k+1) - y_cc(k))
                                dTdz_R = (dTdz*(z_cb(l) - z_cc(l)) + dTdz_R*(z_cc(l+1) - z_cb(l)))/(z_cc(l+1) - z_cc(l))

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
                                    print*, "HIFU: Check thermal properties!"
                                    abortFlag = 1._wp
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
                                    print*, 'dx, dy, dz', dx(j), dy(k), dz(l)
                                    print*, 'x, y, z', x_cc(j), y_cc(k), z_cc(l)
                                    print*, 'Heat sources: ', q_hifu_3d%vf(qus_hifu_idx_ht)%sf(j, k, l), q_hifu_3d%vf(hifu_params%qvis_idx)%sf(j, k, l), &
                                                                q_hifu_3d%vf(hifu_params%qth_idx)%sf(j, k, l)
                                    print*, 'Diff components:', (1._wp/dx(j))*(dTdx_R - dTdx_L), (1._wp/(y_cc(k)*dy(k)))*((y_cc(k) + 0.5_wp*dy(k))*dTdr_R - &
                                                                        (y_cc(k) - 0.5_wp*dy(k))*dTdr_L), (1._wp/(dz(l)*y_cc(k)**2._wp))*(dTdz_R - dTdz_L)
                                    print*, 'T flux:', dTdx_L, dTdx_R, dTdr_L, dTdr_R, dTdz_L, dTdz_R
                                    print*, 'T field:', q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l), q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l - 1), &
                                                                                                    q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + 1), &
                                                                                                    q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k - 1, l), &
                                                                                                    q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k + 1, l), &
                                                                                                    q_hifu_3d%vf(hifu_params%T_idx)%sf(j - 1, k, l), &
                                                                                                    q_hifu_3d%vf(hifu_params%T_idx)%sf(j + 1, k, l)

                                    print*, "Apply s_pole_correction merging more cells at the pole or reduce dt!!"
                                    abortFlag = 1._wp
                                end if

                                ! Calculate min CFL
                                if (t_step == 0) then
                                    CFL_heat_old = CFL_heat
                                    CFL_heat = max(CFL_heat, tdiff*dt/(dx(j)**2_wp))
                                    CFL_heat = max(CFL_heat, tdiff*dt/(dy(k)**2_wp))
                                    CFL_heat = max(CFL_heat, tdiff*dt/((y_cc(k)*dz(l))**2_wp))
                                    if (CFL_heat_old /= CFL_heat) then ! Find second layers azimuthal Courant number
                                        CFL_heat_2 = tdiff*dt/((y_cc(1)*dz(l))**2_wp)
                                        CFL_heat_3 = tdiff*dt/((y_cc(2)*dz(l))**2_wp)
                                        CFL_heat_4 = tdiff*dt/((y_cc(3)*dz(l))**2_wp)
                                    end if
                                    CFL_heat_max = max(CFL_heat_max, CFL_heat)
                                end if

                                abortFlag_max = max(abortFlag_max, abortFlag)

                            end do
                        end do
                    end do


                    if (t_step == 0) then
                        if (num_procs > 1) then
                            val_tmp = CFL_heat
                            call s_mpi_allreduce_max(val_tmp, CFL_heat)
                            val_tmp = CFL_heat_2
                            call s_mpi_allreduce_max(val_tmp, CFL_heat_2)
                            val_tmp = CFL_heat_3
                            call s_mpi_allreduce_max(val_tmp, CFL_heat_3)
                            val_tmp = CFL_heat_4
                            call s_mpi_allreduce_max(val_tmp, CFL_heat_4)
                        end if
                        if (proc_rank == 0) print*, 'max. CFL. First layer:', CFL_heat, CFL_heat/4._wp, &
                                                                                                CFL_heat/16._wp, CFL_heat/64._wp, CFL_heat/128._wp 
                        ! singe cell at the pole
                        ! merging 2 cells at the pole
                        ! merging 4 cells at the pole
                        ! merging 8 cells at the pole
            
                        if (proc_rank == 0) print*, 'max. CFL. Second layer:', CFL_heat_2, CFL_heat_2/4._wp, &
                                                                            CFL_heat_2/16._wp, CFL_heat_2/64._wp, CFL_heat_2/128._wp
                        if (proc_rank == 0) print*, 'max. CFL. Third layer:', CFL_heat_3, CFL_heat_3/4._wp, &
                                                                            CFL_heat_3/16._wp, CFL_heat_3/64._wp, CFL_heat_3/128._wp
                        if (proc_rank == 0) print*, 'max. CFL. Fourth layer:', CFL_heat_4, CFL_heat_4/4._wp, &
                                                                            CFL_heat_4/16._wp, CFL_heat_4/64._wp, CFL_heat_4/128._wp
                    end if

                else ! 3D cartesian (all stages)

                    call s_populate_variables_buffers(q_hifu%vf, pb, mv, bc_type)

                    !$acc parallel loop collapse(3) gang vector default(present) copyin(qus_hifu_idx_ht, t_step) &
                    !$acc reduction(MAX: abortFlag_max, CFL_heat_max) copy(abortFlag_max, CFL_heat_max)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m

                                abortFlag = 0._wp
                                CFL_heat = -100._wp

                                !<  Zeroing RHS_heat
                                q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) = 0._wp

                                !> Temperature derivatives at the cell center. METHOD: Second order centered difference approximation
                                dTdx = (q_hifu%vf(hifu_params%T_idx)%sf(j + 1, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j - 1, k, l))/(x_cc(j + 1) - x_cc(j - 1))
                                dTdx_L = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j - 2, k, l))/(x_cc(j) - x_cc(j - 2))
                                dTdx_R = (q_hifu%vf(hifu_params%T_idx)%sf(j + 2, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k, l))/(x_cc(j + 2) - x_cc(j))

                                dTdr = (q_hifu%vf(hifu_params%T_idx)%sf(j, k + 1, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k - 1, l))/(y_cc(k + 1) - y_cc(k - 1))
                                dTdr_L = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k - 2, l))/(y_cc(k) - y_cc(k - 2))
                                dTdr_R = (q_hifu%vf(hifu_params%T_idx)%sf(j, k + 2, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k, l))/(y_cc(k + 2) - y_cc(k))

                                dTdz = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l + 1) - q_hifu%vf(hifu_params%T_idx)%sf(j, k, l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                                dTdz_L = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l) - q_hifu%vf(hifu_params%T_idx)%sf(j, k, l - 2))/(z_cc(l) - z_cc(l - 2))
                                dTdz_R = (q_hifu%vf(hifu_params%T_idx)%sf(j, k, l + 2) - q_hifu%vf(hifu_params%T_idx)%sf(j, k, l))/(z_cc(l + 2) - z_cc(l))

                                !> Find temperature derivatives at the faces of the cell
                                dTdx_L = (dTdx*(x_cc(j) - x_cb(j - 1)) + dTdx_L*(x_cb(j-1) - x_cc(j-1)))/(x_cc(j) - x_cc(j-1))
                                dTdr_L = (dTdr*(y_cc(k) - y_cb(k - 1)) + dTdr_L*(y_cb(k-1) - y_cc(k-1)))/(y_cc(k) - y_cc(k-1))
                                dTdz_L = (dTdz*(z_cc(l) - z_cb(l - 1)) + dTdz_L*(z_cb(l-1) - z_cc(l-1)))/(z_cc(l) - z_cc(l-1))

                                dTdx_R = (dTdx*(x_cb(j) - x_cc(j)) + dTdx_R*(x_cc(j+1) - x_cb(j)))/(x_cc(j+1) - x_cc(j))
                                dTdr_R = (dTdr*(y_cb(k) - y_cc(k)) + dTdr_R*(y_cc(k+1) - y_cb(k)))/(y_cc(k+1) - y_cc(k))
                                dTdz_R = (dTdz*(z_cb(l) - z_cc(l)) + dTdz_R*(z_cc(l+1) - z_cb(l)))/(z_cc(l+1) - z_cc(l))

                                !> Get thermal properties
                                rho_cp = 0._wp
                                tdiff = 0._wp
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha = q_cons_vf(advxb + i - 1)%sf(j, k, l)
                                    rho_cp = rho_cp + alpha * rho_cp_fluids(i)
                                    tdiff = tdiff + alpha * tdiff_fluids(i)
                                end do

                                if (f_is_default(rho_cp) .or. f_is_default(tdiff)) then
                                    print *, 'alpha, rho_cp, tdiff', alpha, rho_cp, tdiff
                                    print *, "HIFU: Check thermal properties!"
                                    abortFlag = 1._wp
                                end if

                                !> Obtain rhs (see notes)
                                q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) = & 
                                            q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) + tdiff*( &
                                            (1._wp/dx(j))*(dTdx_R - dTdx_L) + &
                                            (1._wp/dy(k))*(dTdr_R - dTdr_L) + &
                                            (1._wp/dz(l))*(dTdz_R - dTdz_L))


                                if ((q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l) > 0._wp) .and. &
                                                            (t_step < hifu_params%stepStopSource)) then
                                    !> Adding the heat source terms
                                    q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) = &
                                            q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) + &
                                            (1._wp/(rho_cp*q_hifu%vf(hifu_params%tsamp_idx)%sf(j, k, l)))*( &
                                            q_hifu%vf(qus_hifu_idx_ht)%sf(j, k, l) + &         !Acoustic intensity
                                            q_hifu%vf(hifu_params%qvis_idx)%sf(j, k, l) + &    !Viscous intensity
                                            q_hifu%vf(hifu_params%qth_idx)%sf(j, k, l))        !Thermal intensity


                                    if (hifu_params%streaming) then
                                        !> Convected heat flux
                                        print*, "HIFU: No streaming valid for 3D heat solver!"
                                        abortFlag = 1._wp
                                    end if
                                end if

                                !Checking NaNs
                                if (q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l) /= &
                                                q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l)) then
                                    print*, 'NaNs in q hifu rhs', q_hifu%vf(hifu_params%T_idx + 1)%sf(j, k, l), j, k, l
                                    print*, 'dx, dy, dz', dx(j), dy(k), dz(l)
                                    print*, 'x, y, z', x_cc(j), y_cc(k), z_cc(l)
                                    print*, 'Heat sources: ', q_hifu%vf(qus_hifu_idx_ht)%sf(j, k, l), &
                                                            q_hifu%vf(hifu_params%qvis_idx)%sf(j, k, l), &
                                                            q_hifu%vf(hifu_params%qth_idx)%sf(j, k, l)
                                    print*, "Reduce dt!!"
                                    abortFlag = 1._wp
                                end if

                                ! Calculate min CFL
                                if (t_step == 0) then
                                    CFL_heat = max(CFL_heat, tdiff*dt/(dx(j)**2_wp))
                                    CFL_heat = max(CFL_heat, tdiff*dt/(dy(k)**2_wp))
                                    CFL_heat = max(CFL_heat, tdiff*dt/(dz(l)**2_wp))
                                    CFL_heat_max = max(CFL_heat_max, CFL_heat)
                                end if

                                abortFlag_max = max(abortFlag_max, abortFlag)

                            end do
                        end do
                    end do

                end if
            end if
        end if

        if (num_procs > 1) then
            val_tmp = abortFlag_max
            call s_mpi_allreduce_max(val_tmp, abortFlag_max)
            val_tmp = CFL_heat_max
            call s_mpi_allreduce_max(val_tmp, CFL_heat_max)
        end if

        if (t_step == 0 .and. proc_rank == 0) print*, 'Max CFL:', CFL_heat_max

        if (abortFlag_max > 0._wp) call s_mpi_abort("Errors found in s_rhs_heatEqn")

        if (proc_rank == 0 .and. t_step == hifu_params%stepStopSource) print *, 'WARNING :: Turn off HIFU source'

    end subroutine s_rhs_heatEqn ! =============================================

!     subroutine s_pole_correction(rhs_heat, j, k, l, t_step)
! #ifdef _CRAYFTN
!     !DIR$ INLINEALWAYS s_get_char_vol
! #else
!     !$acc routine seq
! #endif
!         real(wp), intent(inout) :: rhs_heat
!         integer, intent(in) :: j, k, l, t_step

!         integer :: sub_id, nCells
!         integer :: q, qq
!         real(wp) :: volCell, totVol, Nr

!         if (bc_pole /= -14) return

!         ! Number of cells to merge
!         Nr = ceiling(2._wp*pi*y_cc(k)/(y_cb(k) - y_cb(k - 1)))
!         Nr = ceiling(log(Nr)/log(2._wp))
!         nCells = int(ceiling((p+1)/(2**Nr)))

!         if (nCells == 1) return

!         ! Find cells to merge
!         sub_id = 0
!         do while (.true.)
!             if (mod(l + sub_id, nCells) == 0) then
!                 exit
!             end if
!             sub_id = sub_id + 1
!         end do

!         rhs_heat = 0._wp
!         totVol = 0._wp

!         if (sub_id==0) then

!             if (j==0 .and. l==0 .and. t_step==0) print*, k, nCells, proc_rank

!             !$acc loop seq
!             do q = 0, nCells-1

!                 ! Check temperature is the same in the set of cells
!                 if (q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l) /= &
!                     q_hifu_3d%vf(hifu_params%T_idx)%sf(j, k, l + q)) then

!                     print*, j, k, l, nCells
!                     stop "Different temperatures in the unified cells!!"
!                 end if

!                 ! Calculate volume
!                 volCell = dx(j)*dy(k)*y_cc(k)*dz(l + q)

!                 ! Corrected RHS
!                 rhs_heat = rhs_heat + &
!                             q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l + q) * volCell
!                 totVol = totVol + volCell
                
!             end do

!         else
            
!             qq = sub_id - nCells

!             !$acc loop seq
!             do q = 0, nCells-1

!                 ! Calculate volume
!                 volCell = dx(j)*dy(k)*y_cc(k)*dz(l + qq)

!                 ! Corrected RHS
!                 rhs_heat = rhs_heat + &
!                             q_hifu_3d%vf(hifu_params%T_idx + 1)%sf(j, k, l + qq) * volCell
!                 totVol = totVol + volCell

!                 qq = qq + 1 
                
!             end do

!         end if         
        
!         rhs_heat = rhs_heat / totVol

!     end subroutine s_pole_correction

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
            @:DEALLOCATE(q_hifu%vf(i)%sf)
        end do
        @:DEALLOCATE(q_hifu%vf)

        @:DEALLOCATE(shear_viscous_fluids)
        @:DEALLOCATE(bulk_viscous_fluids)
        @:DEALLOCATE(abs_coef_fluids)
        @:DEALLOCATE(rho_cp_fluids)
        @:DEALLOCATE(tdiff_fluids)

    end subroutine s_finalize_HIFU_module

end module m_hifu
