!>
!! @file m_bubbles_EL.fpp
!! @brief Contains module m_bubbles_EL

#:include 'macros.fpp'

!> @brief This module is used to to compute the volume-averaged bubble model
module m_bubbles_EL

    use m_global_parameters             !< Definitions of the global parameters

    use m_mpi_proxy                     !< Message passing interface (MPI) module proxy

    use m_bubbles_EL_kernels            !< Definitions of the kernel functions

    use m_bubbles                       !< General bubble dynamics procedures

    use m_variables_conversion          !< State variables type conversion procedures

    use m_compile_specific

    use m_boundary_common

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_sim_helpers

    use m_helper

    implicit none

    !(nBub)
    integer, allocatable, dimension(:, :) :: lag_id          !< Global and local IDs
    real(wp), allocatable, dimension(:) :: bub_R0            !< Initial bubble radius
    real(wp), allocatable, dimension(:) :: Rmax_stats        !< Maximum radius
    real(wp), allocatable, dimension(:) :: Rmin_stats        !< Minimum radius
    $:GPU_DECLARE(create='[lag_id, bub_R0, Rmax_stats, Rmin_stats]')

    real(wp), allocatable, dimension(:) :: gas_mg            !< Bubble's gas mass
    real(wp), allocatable, dimension(:) :: gas_betaT         !< heatflux model (Preston et al., 2007)
    real(wp), allocatable, dimension(:) :: gas_betaC         !< massflux model (Preston et al., 2007)
    real(wp), allocatable, dimension(:) :: bub_dphidt        !< subgrid velocity potential (Maeda & Colonius, 2018)
    $:GPU_DECLARE(create='[gas_mg, gas_betaT, gas_betaC, bub_dphidt]')

    !(nBub, 1 -> actual val or 2 -> temp val)
    real(wp), allocatable, dimension(:, :) :: gas_p          !< Pressure in the bubble (Polytropic 1-> actual, 2-> initial)
    real(wp), allocatable, dimension(:, :) :: gas_mv         !< Vapor mass in the bubble
    real(wp), allocatable, dimension(:, :) :: intfc_rad      !< Bubble radius
    real(wp), allocatable, dimension(:, :) :: intfc_vel      !< Velocity of the bubble interface
    real(wp), allocatable, dimension(:, :) :: intfc_ac       !< Acceleration of the bubble interface
    $:GPU_DECLARE(create='[gas_p, gas_mv, intfc_rad, intfc_vel, intfc_ac]')

    !(nBub, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(wp), allocatable, dimension(:, :, :) :: mtn_pos     !< Bubble's position
    real(wp), allocatable, dimension(:, :, :) :: mtn_posPrev !< Bubble's previous position
    real(wp), allocatable, dimension(:, :, :) :: mtn_vel     !< Bubble's velocity
    real(wp), allocatable, dimension(:, :, :) :: mtn_s       !< Bubble's computational cell position in real format
    $:GPU_DECLARE(create='[mtn_pos, mtn_posPrev, mtn_vel, mtn_s]')

    !(nBub, 1-> x or 2->y or 3 ->z, time-stage)
    real(wp), allocatable, dimension(:, :) :: intfc_draddt   !< Time derivative of bubble's radius
    real(wp), allocatable, dimension(:, :) :: intfc_dveldt   !< Time derivative of bubble's interface velocity
    real(wp), allocatable, dimension(:, :) :: gas_dpdt       !< Time derivative of gas pressure
    real(wp), allocatable, dimension(:, :) :: gas_dmvdt      !< Time derivative of the vapor mass in the bubble
    ! real(wp), allocatable, dimension(:, :, :) :: mtn_dposdt  !< Time derivative of the bubble's position
    ! real(wp), allocatable, dimension(:, :, :) :: mtn_dveldt  !< Time derivative of the bubble's velocity
    $:GPU_DECLARE(create='[intfc_draddt, intfc_dveldt, gas_dpdt, gas_dmvdt]')

    real(wp), allocatable, dimension(:) :: bub_interact        !< Scattered pressure from each bubble
    real(wp), allocatable, dimension(:, :) :: bub_int_ids     !< Ids of the neighboring bubbles for pout interaction
    !real(wp), allocatable, dimension(:) :: bub_lambda_c      !< Mean inter-bubble distance (p' white noise)
    !real(wp), allocatable, dimension(:, :) :: bub_rnd_phase  !< Random phases (1:num_noise) per bubble (p' white noise)
    $:GPU_DECLARE(create='[bub_interact, bub_int_ids]')

    integer, private :: lag_num_ts      !<  Number of time stages in the time-stepping scheme
    integer :: nBubs                    !< Number of bubbles in the local domain
    real(wp) :: Rmax_glb, Rmin_glb, Rmean_glb   !< Global stats of bubbe size in the local domain
    type(vector_field) :: q_beta        !< Projection of the lagrangian particles in the Eulerian framework
    integer :: q_beta_idx               !< Size of the q_beta vector field
    $:GPU_DECLARE(create='[nBubs,lag_num_ts,Rmax_glb,Rmin_glb,Rmean_glb,q_beta,q_beta_idx]')

    real(wp), allocatable, dimension(:,:) :: moments_bubs   !< Moments of volume, qac, qvis, qth_pos, qth_neg
    real(wp), allocatable, dimension(:) :: acPw_bubs         !< Acoustic power of the bubbles (HIFU)
    real(wp), allocatable, dimension(:, :) :: mrmtnt_shell  !< Lipid shell indicator (Marmotant model)
    real(wp), allocatable, dimension(:) :: mrmtnt_Rbuck     !< Buckling radius (Marmotant model)
    real(wp), allocatable, dimension(:) :: mrmtnt_Rrupt     !< Rupture radius (Marmotant model)
    real(wp), allocatable, dimension(:) :: bub_qvis         !< Time-averaged viscous intensity (HIFU)
    real(wp), allocatable, dimension(:) :: bub_qth          !< Time-averaged thermal intensity (HIFU)
    real(wp), allocatable, dimension(:) :: bub_hifu_rad     !< Time-averaged radius
    $:GPU_DECLARE(create='[mrmtnt_shell, mrmtnt_Rbuck, mrmtnt_Rrupt, bub_qvis, bub_qth, bub_hifu_rad, moments_bubs, acPw_bubs]')

    

contains

    !> Initializes the lagrangian subgrid bubble solver
        !! @param q_cons_vf Initial conservative variables
    impure subroutine s_initialize_bubbles_EL_module(q_cons_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        integer :: nBubs_glb, i, int_var

        ! Setting number of time-stages for selected time-stepping scheme
        lag_num_ts = time_stepper

        ! Allocate space for the Eulerian fields needed to map the effect of the bubbles
        if (lag_params%solver_approach == 1) then
            ! One-way coupling
            q_beta_idx = 3
        elseif (lag_params%solver_approach == 2) then
            ! Two-way coupling
            q_beta_idx = 4
            if (p == 0) then
                !Subgrid noise model for 2D approximation
                q_beta_idx = 6
            end if
        else
            call s_mpi_abort('Please check the lag_params%solver_approach input')
        end if

        $:GPU_UPDATE(device='[lag_num_ts, q_beta_idx]')

        @:ALLOCATE(q_beta%vf(1:q_beta_idx))

        do i = 1, q_beta_idx
            @:ALLOCATE(q_beta%vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
        end do

        @:ACC_SETUP_VFs(q_beta)

        ! Allocating space for lagrangian variables
        nBubs_glb = lag_params%nBubs_glb

        @:ALLOCATE(lag_id(1:nBubs_glb, 1:2))
        @:ALLOCATE(bub_R0(1:nBubs_glb))
        @:ALLOCATE(Rmax_stats(1:nBubs_glb))
        @:ALLOCATE(Rmin_stats(1:nBubs_glb))
        @:ALLOCATE(gas_mg(1:nBubs_glb))
        @:ALLOCATE(gas_betaT(1:nBubs_glb))
        @:ALLOCATE(gas_betaC(1:nBubs_glb))
        @:ALLOCATE(bub_dphidt(1:nBubs_glb))
        @:ALLOCATE(gas_p(1:nBubs_glb, 1:2))
        @:ALLOCATE(gas_mv(1:nBubs_glb, 1:2))
        @:ALLOCATE(intfc_rad(1:nBubs_glb, 1:2))
        @:ALLOCATE(intfc_vel(1:nBubs_glb, 1:2))
        @:ALLOCATE(intfc_ac(1:nBubs_glb, 1:2))
        @:ALLOCATE(mtn_pos(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_posPrev(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_vel(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_s(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(intfc_draddt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(intfc_dveldt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(gas_dpdt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(gas_dmvdt(1:nBubs_glb, 1:lag_num_ts))
        ! @:ALLOCATE(mtn_dposdt(1:nBubs_glb, 1:3, 1:lag_num_ts))
        ! @:ALLOCATE(mtn_dveldt(1:nBubs_glb, 1:3, 1:lag_num_ts))
        ! Marmotant model
        @:ALLOCATE(mrmtnt_shell(1:nBubs_glb, 1:2))
        @:ALLOCATE(mrmtnt_Rbuck(1:nBubs_glb))
        @:ALLOCATE(mrmtnt_Rrupt(1:nBubs_glb))

        ! hifu
        @:ALLOCATE(bub_qvis(1:nBubs_glb))
        @:ALLOCATE(bub_qth(1:nBubs_glb))
        @:ALLOCATE(bub_hifu_rad(1:nBubs_glb))

        if (hifu_params%moments) then
            @:ALLOCATE(moments_bubs(1:4, 1:4))
        end if
        if (hifu_params%power_balance) then
            @:ALLOCATE(acPw_bubs(1:3))
        end if

        ! Interbubble interaction
        ! 1: emitted Pout, 2: sum of Pouts from volume of influence (self-inclusive)
        @:ALLOCATE(bub_interact(1:nBubs_glb))
        ! 1: number of interacting bubbles (self-inclusive), 2:nBubs_glb+1: IDs in the volume of influence (self-inclusive)
        if (lag_params%pressure_corrector .and. any(lag_params%interaction_model == (/2, 3/))) then
            int_var = min(max_bub_int, nBubs_glb+1)
            @:ALLOCATE(bub_int_ids(1:nBubs_glb, 1:int_var))
        end if
        !@:ALLOCATE(bub_lambda_c(1:nBubs_glb))
        !@:ALLOCATE(bub_rnd_phase(1:nBubs_glb, 1:num_noise))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

        call s_initialize_bubbles_EL_kernels()

        ! Starting bubbles
        call s_start_lagrange_inputs()
        call s_read_input_bubbles(q_cons_vf, bc_type)

    end subroutine s_initialize_bubbles_EL_module

    !> The purpose of this procedure is to start lagrange bubble parameters applying nondimensionalization if needed
    impure subroutine s_start_lagrange_inputs()

        integer :: id_bubbles, id_host
        real(wp) :: rho0, c0, T0, x0, p0

        id_bubbles = num_fluids
        id_host = num_fluids - 1

        !Reference values
        rho0 = lag_params%rho0
        c0 = lag_params%c0
        T0 = lag_params%T0
        x0 = lag_params%x0
        p0 = rho0*c0*c0

        !Update inputs
        Tw = lag_params%Thost/T0
        pv = fluid_pp(id_host)%pv/p0
        gamma_v = fluid_pp(id_host)%gamma_v
        gamma_n = fluid_pp(id_bubbles)%gamma_v
        k_vl = fluid_pp(id_host)%k_v*(T0/(x0*rho0*c0*c0*c0))
        k_nl = fluid_pp(id_bubbles)%k_v*(T0/(x0*rho0*c0*c0*c0))
        cp_v = fluid_pp(id_host)%cp_v*(T0/(c0*c0))
        cp_n = fluid_pp(id_bubbles)%cp_v*(T0/(c0*c0))
        R_v = (R_uni/fluid_pp(id_host)%M_v)*(T0/(c0*c0))
        R_n = (R_uni/fluid_pp(id_bubbles)%M_v)*(T0/(c0*c0))
        lag_params%diffcoefvap = lag_params%diffcoefvap/(x0*c0)
        ss = fluid_pp(id_host)%ss/(rho0*x0*c0*c0)
        mul0 = fluid_pp(id_host)%mul0/(rho0*x0*c0)

        !Marmotant model
        lag_params%ss0_ctdBub = lag_params%ss0_ctdBub/(rho0*x0*c0*c0)
        lag_params%srfDilVsc_ctdBub = lag_params%srfDilVsc_ctdBub/(rho0*x0*x0*c0)
        lag_params%srfElast_ctdBub = lag_params%srfElast_ctdBub/(rho0*x0*c0*c0)

        ! Parameters used in bubble_model
        Web = 1._wp/ss
        Re_inv = mul0

        if (polytropic) then
          Ca = (p0-pv)/(rho0*c0*c0)
          gamma_m = gamma_n
          if (thermal == 2) gamma_m = 1._wp ! Isothermal
        end if

        ! Need improvements to accept polytropic gas compression, isothermal and adiabatic thermal models, and
        ! the Gilmore and RP bubble models.
        ! polytropic = .false.    ! Forcing no polytropic model
        ! thermal = 3             ! Forcing constant transfer coefficient model based on Preston et al., 2007
        ! If Keller-Miksis model is not selected, then no radial motion

        !GPU vars get updated in initialize_gpu_vars

    end subroutine s_start_lagrange_inputs

    !> The purpose of this procedure is to obtain the initial bubbles' information
        !! @param q_cons_vf Conservative variables
    impure subroutine s_read_input_bubbles(q_cons_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        real(wp), dimension(8) :: inputBubble
        real(wp) :: qtime
        integer :: id, bub_id, save_count
        integer :: i, ios
        logical :: file_exist, read_flag, indomain
        real(wp) :: safeStop, tmp_val
        character(LEN=path_len + 2*name_len) :: path_D_dir, file_loc !<

        ! Initialize number of particles
        bub_id = 0
        id = 0
        safeStop = 0._wp

        ! Read the input lag_bubble file or restart point
        if (cfl_dt) then
            save_count = n_start
            qtime = n_start*t_save
        else
            save_count = t_step_start
            qtime = t_step_start*dt
        end if

        ! Read input file in the middle of a pure Euler simulation
        write (file_loc, '(a,i0,a)') 'lag_bubbles_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=file_exist)
        read_flag = .true.
        lag_params%initial_corrector = .true.
        if (file_exist) then
            read_flag = .false.
            lag_params%initial_corrector = .false.
        end if

        if (read_flag) then
            if (proc_rank == 0) print *, 'Reading lagrange bubbles input file.'
            call s_mpi_barrier()
            inquire (file='input/lag_bubbles.dat', exist=file_exist)
            if (file_exist) then
                open (94, file='input/lag_bubbles.dat', form='formatted', iostat=ios)
                do while (ios == 0)
                    read (94, *, iostat=ios) (inputBubble(i), i=1, 8)
                    if (ios /= 0) cycle
                    indomain = particle_in_domain_physical(inputBubble(1:3))
                    id = id + 1
                    if (indomain) then
                        bub_id = bub_id + 1
                        if (bub_id > lag_params%nBubs_glb) then
                            safeStop = 1._wp*bub_id
                        else
                            call s_add_bubbles(inputBubble, q_cons_vf, bub_id)
                            lag_id(bub_id, 1) = id      !global ID
                            lag_id(bub_id, 2) = bub_id  !local ID
                            nBubs = bub_id              ! local number of bubbles
                        end if
                    end if
                end do
                close (94)
            else
                call s_mpi_abort("Initialize the lagrange bubbles in input/lag_bubbles.dat")
            end if
        else
            if (proc_rank == 0) print *, 'Restarting lagrange bubbles at save_count: ', save_count
            call s_mpi_barrier()
            call s_restart_bubbles(bub_id, save_count)
        end if

        print '("Lagrange bubbles running, in proc ", I8, " number: ", I8, " / ", I8)', proc_rank, bub_id, id

        call s_mpi_barrier()
        if (num_procs > 1) then
            call s_mpi_allreduce_max(safeStop, tmp_val)
            safeStop = tmp_val
        end if
        
        if (int(safeStop) > lag_params%nBubs_glb) then
            if (proc_rank == 0) print '("Maximum number of bubbbles per processor is: ", I8)', int(safeStop)
            call s_mpi_abort('Current number of bubbles is larger than nBubs_glb.')
        else
            safeStop = nBubs
            if (num_procs > 1) then
                call s_mpi_allreduce_max(safeStop, tmp_val)
                safeStop = tmp_val
            end if
            if (proc_rank == 0) print '("Maximum number of bubbbles per processor is: ", I8)', int(safeStop)
        end if

        $:GPU_UPDATE(device='[bubbles_lagrange, lag_params]')

        $:GPU_UPDATE(device='[lag_id,bub_R0,Rmax_stats,Rmin_stats,gas_mg, &
            & gas_betaT,gas_betaC,bub_dphidt,gas_p,gas_mv, &
            & intfc_rad,intfc_vel,mtn_pos,mtn_posPrev,mtn_vel, &
            & mtn_s,intfc_draddt,intfc_dveldt,gas_dpdt,gas_dmvdt, nBubs]')
            ! & mtn_dposdt,mtn_dveldt]')
        
        $:GPU_UPDATE(device='[intfc_ac, mrmtnt_shell, mrmtnt_Rbuck, &
            & mrmtnt_Rrupt, bub_qvis, bub_qth, bub_hifu_rad]')


        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)
        Rmean_glb = 0._wp
        $:GPU_UPDATE(device='[Rmax_glb, Rmin_glb, Rmean_glb]')

        $:GPU_UPDATE(device='[dx,dy,dz,x_cb,x_cc,y_cb,y_cc,z_cb,z_cc]')

        !Populate temporal variables
        call s_transfer_data_to_tmp
        call s_start_bubble_interaction
        call s_smear_voidfraction(bc_type)

        if (read_flag) then
            ! Create ./D directory
            write (path_D_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
            call my_inquire(path_D_dir, file_exist)
            if (.not. file_exist) call s_create_directory(trim(path_D_dir))
            call s_write_restart_lag_bubbles(save_count) ! Needed for post_processing
        end if

        call s_calculate_lag_bubble_stats()
        if (lag_params%write_bubbles) call s_write_lag_particles(qtime, replace=.false.)
        call s_write_void_evol(qtime, replace=.false.)

    end subroutine s_read_input_bubbles

    !> The purpose of this procedure is to obtain the information of the bubbles when starting fresh
        !! @param inputBubble Bubble information
        !! @param q_cons_vf Conservative variables
        !! @param bub_id Local id of the bubble
    impure subroutine s_add_bubbles(inputBubble, q_cons_vf, bub_id)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(8), intent(in) :: inputBubble
        integer, intent(in) :: bub_id
        integer :: i

        real(wp) :: pliq, volparticle, concvap, totalmass, kparticle, cpparticle
        real(wp) :: omegaN_local, PeG, PeT, rhol, pcrit, qv, gamma, pi_inf, dynP
        integer, dimension(3) :: cell
        real(wp), dimension(2) :: Re
        real(wp) :: massflag, heatflag, Re_trans, Im_trans, Web_mod

        bub_R0(bub_id) = inputBubble(7)
        Rmax_stats(bub_id) = min(dflt_real, -dflt_real)
        Rmin_stats(bub_id) = max(dflt_real, -dflt_real)
        bub_dphidt(bub_id) = 0._wp
        intfc_rad(bub_id, 1) = inputBubble(7)
        intfc_vel(bub_id, 1) = inputBubble(8)
        intfc_ac(bub_id, 1) = 0._wp
        mtn_pos(bub_id, 1:3, 1) = inputBubble(1:3)
        mtn_posPrev(bub_id, 1:3, 1) = mtn_pos(bub_id, 1:3, 1)
        mtn_vel(bub_id, 1:3, 1) = inputBubble(4:6)
        bub_qvis(bub_id) = 0._wp
        bub_qth(bub_id) = 0._wp
        bub_hifu_rad(bub_id) = 0._wp
        bub_interact(bub_id) = 0._wp

        if (cyl_coord .and. p == 0) then
            mtn_pos(bub_id, 2, 1) = sqrt(mtn_pos(bub_id, 2, 1)**2._wp + &
                                         mtn_pos(bub_id, 3, 1)**2._wp)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            mtn_pos(bub_id, 3, 1) = atan2(inputBubble(3), inputBubble(2))
            !mtn_posPrev(bub_id, 1:3, 1) = mtn_pos(bub_id, 1:3, 1) ! Need 3D coords for hifu heat solver
        end if

        cell = -buff_size
        call s_locate_cell(mtn_pos(bub_id, 1:3, 1), cell, mtn_s(bub_id, 1:3, 1))

        ! Check if the bubble is located in the ghost cell of a symmetric boundary
        if ((bc_x%beg == BC_REFLECTIVE .and. cell(1) < 0) .or. &
            (bc_x%end == BC_REFLECTIVE .and. cell(1) > m) .or. &
            (bc_y%beg == BC_REFLECTIVE .and. cell(2) < 0) .or. &
            (bc_y%end == BC_REFLECTIVE .and. cell(2) > n)) then
            call s_mpi_abort("Lagrange bubble is in the ghost cells of a symmetric boundary.")
        end if

        if (p > 0) then
            if ((bc_z%beg == BC_REFLECTIVE .and. cell(3) < 0) .or. &
                (bc_z%end == BC_REFLECTIVE .and. cell(3) > p)) then
                call s_mpi_abort("Lagrange bubble is in the ghost cells of a symmetric boundary.")
            end if
        end if

        ! If particle is in the ghost cells, find the closest non-ghost cell
        cell(1) = min(max(cell(1), 0), m)
        cell(2) = min(max(cell(2), 0), n)
        if (p > 0) cell(3) = min(max(cell(3), 0), p)
        call s_convert_to_mixture_variables(q_cons_vf, cell(1), cell(2), cell(3), &
                                            rhol, gamma, pi_inf, qv, Re)
        dynP = 0._wp
        do i = 1, num_dims
            dynP = dynP + 0.5_wp*q_cons_vf(contxe + i)%sf(cell(1), cell(2), cell(3))**2/rhol
        end do
        if (.not. f_is_default(acoustic_bc_params%Pbase)) then
            pliq = acoustic_bc_params%Pbase
        else
            pliq = (q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - dynP - pi_inf)/gamma
        end if
        if (pliq < 0) print *, "Negative pressure", proc_rank, &
            q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)), pi_inf, gamma, pliq, cell, dynP

        ! Activate or deactivate the mass model
        massflag = 0._wp
        !if (lag_params%coatedBub_model .or. lag_params%massTransfer_model) then
        if (lag_params%massTransfer_model) then
            !Assume vapor and gas is present in bubble
            massflag = 1._wp
        end if

        ! Marmotant model parameters
        mrmtnt_shell(bub_id, 1) = 0._wp
        if (lag_params%coatedBub_model) mrmtnt_shell(bub_id, 1) = 1._wp
        mrmtnt_Rbuck(bub_id) = mrmtnt_shell(bub_id, 1)*bub_R0(bub_id)/sqrt(1._wp + &
                                                                           lag_params%ss0_ctdBub/lag_params%srfElast_ctdBub)
        mrmtnt_Rrupt(bub_id) = mrmtnt_Rbuck(bub_id)*sqrt(1._wp + ss/lag_params%srfElast_ctdBub)

        ! Initial particle pressure
        gas_p(bub_id, 1) = pliq + 2._wp*(1._wp/Web)/bub_R0(bub_id)
        if (lag_params%coatedBub_model) then
            gas_p(bub_id, 1) = pliq + 2._wp*(lag_params%ss0_ctdBub)/bub_R0(bub_id)
            !print *, 'Rbuck and Rrupt', mrmtnt_Rbuck(bub_id), mrmtnt_Rrupt(bub_id), bub_id
        end if
        if (pv*(massflag) > gas_p(bub_id, 1)) then
            print*, proc_rank, gas_p(bub_id, 1), pv*(massflag), pliq
            call s_mpi_abort("Lagrange bubble initially located in a region with pressure below the vapor pressure.")
        end if

        ! Initial particle mass
        volparticle = 4._wp/3._wp*pi*bub_R0(bub_id)**3._wp ! volume
        gas_mv(bub_id, 1) = pv*volparticle*(1._wp/(R_v*Tw))*(massflag) ! vapermass
        gas_mg(bub_id) = (gas_p(bub_id, 1) - pv*(massflag))*volparticle*(1._wp/(R_n*Tw)) ! gasmass
        if (gas_mg(bub_id) <= 0._wp) then
            call s_mpi_abort("The initial mass of gas inside the bubble is negative. Check the initial conditions.")
        end if
        totalmass = gas_mg(bub_id) + gas_mv(bub_id, 1) ! totalmass

        ! Bubble natural frequency
        concvap = gas_mv(bub_id, 1)/(gas_mv(bub_id, 1) + gas_mg(bub_id))
        omegaN_local = (3._wp*(gas_p(bub_id, 1) - pv*(massflag)) + 4._wp*(1._wp/Web)/bub_R0(bub_id))/rhol
        if (lag_params%coatedBub_model) then
            omegaN_local = (3._wp*(gas_p(bub_id, 1) - pv*(massflag)) + 4._wp*(lag_params%ss0_ctdBub)/bub_R0(bub_id))/rhol
        end if
        omegaN_local = sqrt(omegaN_local/bub_R0(bub_id)**2._wp)

        cpparticle = concvap*cp_v + (1._wp - concvap)*cp_n
        kparticle = concvap*k_vl + (1._wp - concvap)*k_nl

        ! Mass and heat transfer coefficients (based on Preston 2007)
        PeT = totalmass/volparticle*cpparticle*bub_R0(bub_id)**2._wp*omegaN_local/kparticle
        call s_transcoeff(1._wp, PeT, Re_trans, Im_trans)
        gas_betaT(bub_id) = Re_trans*kparticle

        PeG = bub_R0(bub_id)**2._wp*omegaN_local/lag_params%diffcoefvap
        call s_transcoeff(1._wp, PeG, Re_trans, Im_trans)
        gas_betaC(bub_id) = Re_trans*lag_params%diffcoefvap

        if (polytropic) then
            gas_p(bub_id, 2) = gas_p(bub_id, 1)
        else
            if (gas_betaT(bub_id) /= gas_betaT(bub_id) .or. gas_betaC(bub_id) /= gas_betaC(bub_id)) then
                print *, bub_id, gas_betaT(bub_id), gas_betaC(bub_id)
                call s_mpi_abort("NaN mass and heat transfer coefficients")
            end if
        end if


    end subroutine s_add_bubbles

    subroutine s_initial_pressure_correction(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        real(wp) :: pinf, aux1, aux2, massflag, volparticle
        real(wp) :: concvap, totalmass, kparticle, cpparticle
        real(wp) :: omegaN, PeG, PeT, cson, rhol, Re_trans, Im_trans
        real(wp) :: gamma, pi_inf, qv, myRcell
        real(wp), dimension(contxe) :: myalpha_rho, myalpha
        real(wp), dimension(2) :: Re
        integer, dimension(3) :: cell
        integer :: i, k
        complex(wp) :: imag, trans, c1, c2, c3
        real(wp) :: qtime
        integer :: save_count

        lag_params%initial_corrector = .false.

        if (lag_params%cluster_type /= 1) then

            if (proc_rank == 0) print *, 'Performing s_initial_pressure_correction '

            $:GPU_PARALLEL_LOOP(private='[k, myalpha_rho, myalpha, Re, cell]')
            do k = 1, nBubs
                ! Obtaining driving pressure
                call s_get_pinf(k, q_prim_vf, 1, pinf, cell, aux1, aux2, myRcell)

                ! Obtain liquid density and computing speed of sound from pinf
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    myalpha_rho(i) = q_prim_vf(i)%sf(cell(1), cell(2), cell(3))
                    myalpha(i) = q_prim_vf(E_idx + i)%sf(cell(1), cell(2), cell(3))
                end do
                call s_convert_species_to_mixture_variables_acc(rhol, gamma, pi_inf, qv, myalpha, &
                                                                myalpha_rho, Re)
                call s_compute_cson_from_pinf(q_prim_vf, pinf, cell, rhol, gamma, pi_inf, cson)

                ! Activate or deactivate the mass model
                massflag = 0._wp
                if (lag_params%massTransfer_model) then
                    !Assume vapor and gas is present in bubble
                    massflag = 1._wp
                end if

                ! Marmotant model parameters: no need correction

                ! Initial particle pressure
                gas_p(k, 1) = pinf + 2._wp*(1._wp/Web)/bub_R0(k)
                if (lag_params%coatedBub_model) then
                    gas_p(k, 1) = pinf + 2._wp*(lag_params%ss0_ctdBub)/bub_R0(k)
                end if
                if (polytropic) gas_p(k, 2) = gas_p(k, 1)

                ! Initial particle mass
                volparticle = 4._wp/3._wp*pi*bub_R0(k)**3._wp ! volume
                gas_mv(k, 1) = pv*volparticle*(1._wp/(R_v*Tw))*(massflag) ! vapermass
                gas_mg(k) = (gas_p(k, 1) - pv*(massflag))*volparticle*(1._wp/(R_n*Tw)) ! gasmass
                totalmass = gas_mg(k) + gas_mv(k, 1) ! totalmass

                ! Bubble natural frequency
                concvap = gas_mv(k, 1)/(gas_mv(k, 1) + gas_mg(k))
                omegaN = (3._wp*(gas_p(k, 1) - pv*(massflag)) + 4._wp*(1._wp/Web)/bub_R0(k))/rhol
                if (lag_params%coatedBub_model) then
                    omegaN = (3._wp*(gas_p(k, 1) - pv*(massflag)) + 4._wp*(lag_params%ss0_ctdBub)/bub_R0(k))/rhol
                end if

                omegaN = sqrt(omegaN/bub_R0(k)**2._wp)
                cpparticle = concvap*cp_v + (1._wp - concvap)*cp_n
                kparticle = concvap*k_vl + (1._wp - concvap)*k_nl

                ! Mass and heat transfer coefficients (based on Preston 2007)
                PeT = totalmass/volparticle*cpparticle*bub_R0(k)**2._wp*omegaN/kparticle
                imag = (0._wp, 1._wp)
                c1 = imag*PeT
                c2 = sqrt(c1)
                c3 = (exp(c2) - exp(-c2))/(exp(c2) + exp(-c2)) ! tanh(c2)
                trans = ((c2/c3 - 1._wp)**(-1) - 3._wp/c1)**(-1) ! transfer function
                Re_trans = trans
                Im_trans = aimag(trans)
                gas_betaT(k) = Re_trans*kparticle

                PeG = bub_R0(k)**2._wp*omegaN/lag_params%diffcoefvap
                c1 = imag*PeG
                c2 = sqrt(c1)
                c3 = (exp(c2) - exp(-c2))/(exp(c2) + exp(-c2)) ! tanh(c2)
                trans = ((c2/c3 - 1._wp)**(-1) - 3._wp/c1)**(-1) ! transfer function
                Re_trans = trans
                Im_trans = aimag(trans)
                gas_betaC(k) = Re_trans*lag_params%diffcoefvap

            end do

            $:GPU_UPDATE(device='[lag_id,bub_R0,Rmax_stats,Rmin_stats,gas_mg, &
                & gas_betaT,gas_betaC,bub_dphidt,gas_p,gas_mv, &
                & intfc_rad,intfc_vel,mtn_pos,mtn_posPrev,mtn_vel, &
                & mtn_s,intfc_draddt,intfc_dveldt,gas_dpdt,gas_dmvdt, nBubs]')
                ! & mtn_dposdt,mtn_dveldt]')
            
            $:GPU_UPDATE(device='[intfc_ac, mrmtnt_shell, mrmtnt_Rbuck, &
                & mrmtnt_Rrupt, bub_qvis, bub_qth, bub_hifu_rad]')

            call s_transfer_data_to_tmp
            call s_smear_voidfraction(bc_type)

            !Replace files
            if (cfl_dt) then
                save_count = n_start
                qtime = n_start*t_save
            else
                save_count = t_step_start
                qtime = t_step_start*dt
            end if

            call s_calculate_lag_bubble_stats()
            if (lag_params%write_bubbles) call s_write_lag_particles(qtime, replace=.true.)
            call s_write_restart_lag_bubbles(save_count) ! Needed for post_processing
            call s_write_void_evol(qtime, replace=.true.)

            !call s_mpi_barrier()

        end if

    end subroutine s_initial_pressure_correction

    !> The purpose of this procedure is to obtain the information of the bubbles from a restart point.
        !! @param bub_id Local ID of the particle
        !! @param save_count File identifier
    impure subroutine s_restart_bubbles(bub_id, save_count)

        integer, intent(inout) :: bub_id, save_count

        character(LEN=path_len + 2*name_len) :: file_loc

#ifdef MFC_MPI
        real(wp), dimension(27) :: inputvals
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer :: view

        integer, dimension(3) :: cell
        logical :: indomain, particle_file, file_exist

        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer :: ifile, ierr, tot_data, id
        integer :: i
        integer :: varsExtra = 7

        real(wp):: savedTime, saved_dt

        write (file_loc, '(a,i0,a)') 'lag_bubbles_mpi_io_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=file_exist)

        if (file_exist) then
            if (proc_rank == 0) then
                open (9, file=trim(file_loc), form='unformatted', status='unknown')
                read (9) tot_data, mytime, saved_dt
                close (9)
            end if
        else
            print '(a)', trim(file_loc)//' is missing. exiting.'
            call s_mpi_abort
        end if

        call MPI_BCAST(tot_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        gsizes(1) = tot_data
        gsizes(2) = 21 + varsExtra
        lsizes(1) = tot_data
        lsizes(2) = 21 + varsExtra
        start_idx_part(1) = 0
        start_idx_part(2) = 0

        call MPI_type_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        ! Open the file to write all flow variables
        write (file_loc, '(a,i0,a)') 'lag_bubbles_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=particle_file)

        if (particle_file) then
            call MPI_FILE_open(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)
            disp = 0._wp
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, &
                                   'native', mpi_info_null, ierr)
            allocate (MPI_IO_DATA_lag_bubbles(tot_data, 1:(21 + varsExtra)))
            call MPI_FILE_read_ALL(ifile, MPI_IO_DATA_lag_bubbles, (21 + varsExtra)*tot_data, &
                                   mpi_p, status, ierr)
            do i = 1, tot_data
                id = int(MPI_IO_DATA_lag_bubbles(i, 1))
                inputvals(1:(20 + varsExtra)) = MPI_IO_DATA_lag_bubbles(i, 2:(21 + varsExtra))
                indomain = particle_in_domain_physical(inputvals(1:3))
                if (indomain .and. (id > 0)) then
                    bub_id = bub_id + 1
                    nBubs = bub_id                  ! local number of bubbles
                    lag_id(bub_id, 1) = id          ! global ID
                    lag_id(bub_id, 2) = bub_id      ! local ID
                    mtn_pos(bub_id, 1:3, 1) = inputvals(1:3)
                    mtn_posPrev(bub_id, 1:3, 1) = inputvals(4:6)
                    mtn_vel(bub_id, 1:3, 1) = inputvals(7:9)
                    intfc_rad(bub_id, 1) = inputvals(10)
                    intfc_vel(bub_id, 1) = inputvals(11)

                    bub_R0(bub_id) = inputvals(12)
                    Rmax_stats(bub_id) = inputvals(13)
                    Rmin_stats(bub_id) = inputvals(14)
                    bub_dphidt(bub_id) = inputvals(15)
                    gas_p(bub_id, 1) = inputvals(16)
                    gas_mv(bub_id, 1) = inputvals(17)
                    gas_mg(bub_id) = inputvals(18)
                    gas_betaT(bub_id) = inputvals(19)
                    gas_betaC(bub_id) = inputvals(20)
                    ! Marmotant
                    mrmtnt_shell(bub_id, 1) = inputvals(21)
                    mrmtnt_Rbuck(bub_id) = inputvals(22)
                    mrmtnt_Rrupt(bub_id) = inputvals(23)
                    ! hifu
                    bub_qvis(bub_id) = inputvals(24)
                    bub_qth(bub_id) = inputvals(25)
                    intfc_ac(bub_id, 1) = inputvals(26)
                    bub_hifu_rad(bub_id) = inputvals(27)

                    bub_interact(bub_id) = 0._wp
                    if (polytropic) then
                        gas_p(bub_id, 2) = gas_p(bub_id, 1)
                        gas_p(bub_id, 1) = pv + (gas_p(bub_id, 2) - pv)*(bub_R0(bub_id)/intfc_rad(bub_id, 1))**(3._wp*gamma_m)
                    end if
                    cell = -buff_size
                    call s_locate_cell(mtn_pos(bub_id, 1:3, 1), cell, mtn_s(bub_id, 1:3, 1))
                end if
            end do
            deallocate (MPI_IO_DATA_lag_bubbles)
        end if
        call MPI_FILE_CLOSE(ifile, ierr)
#endif

    end subroutine s_restart_bubbles

    !>  3D modeling of bubble interaction and 2D approximation.
        !!      3D: Obtain the identifiers of the bubbles within the volume of influence. Assume that the bubbles
        !!      are smaller than the grid size always. (no need to re-run during simulation).
        !!      2D: Compute the local number density of each bubble for the p'_cell model as white noise.
        !!      It is the number of bubbles per volume of mixture in the physical domain (not buffers).
    subroutine s_start_bubble_interaction()

        integer :: i, j, k
        integer :: nb_local
        real(wp) :: xb_smear, xe_smear
        real(wp) :: yb_smear, ye_smear
        real(wp) :: zb_smear, ze_smear
        real(wp), dimension(3) :: scoord
        integer, dimension(3) :: cell
        integer :: smear_idx
        real(wp) :: num_rn1, num_rn2, num_rn
        real(wp) :: st_dev_rn, mean_rn
        real(wp) :: safeStop, tmp_val

        if (.not. lag_params%pressure_corrector) return

        if (proc_rank == 0) print '("Influence volume (bubble interaction) in # of surrounding cells is ", I0)', lag_params%influence
        ! mean_rn = 0.5_wp*pi
        ! st_dev_rn = 1._wp

        safeStop = 0._wp
        $:GPU_PARALLEL_LOOP(private='[j, cell, scoord]', &
        & reduction='[[safeStop]]',reductionOp='[MAX]',copy='[safeStop]')
        do j = 1, nBubs

            ! Is the bubble in the physical domain?
            !if (particle_in_domain_physical(mtn_pos(j, 1:3, 1))) then

            ! Find the cell location
            scoord = mtn_s(j, 1:3, 1)
            cell(:) = int(scoord(:))
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_dims
                if (scoord(i) < 0._wp) cell(i) = cell(i) - 1
            end do

            ! Define smearing boundaries
            !   Assuming that the cell is always larger than the bubble, then
            !   the smearing volume is constant (influence+1+influence)x(3+1+3).

            smear_idx = cell(1) - lag_params%influence - 1
            if (smear_idx < -buff_size - 1) then
                do while (smear_idx < -buff_size - 1)
                    smear_idx = smear_idx + 1
                end do
            end if
            xb_smear = x_cb(smear_idx)

            smear_idx = cell(1) + lag_params%influence
            if (smear_idx > m + buff_size) then
                do while (smear_idx > m + buff_size)
                    smear_idx = smear_idx - 1
                end do
            end if
            xe_smear = x_cb(smear_idx)

            smear_idx = cell(2) - lag_params%influence - 1
            if (smear_idx < -buff_size - 1 .and. .not. cyl_coord) then
                do while (smear_idx < -buff_size - 1)
                    smear_idx = smear_idx + 1
                end do
            end if
            if (smear_idx < -1 .and. cyl_coord) smear_idx = -1
            yb_smear = y_cb(smear_idx)

            smear_idx = cell(2) + lag_params%influence
            if (smear_idx > n + buff_size) then
                do while (smear_idx > n + buff_size)
                    smear_idx = smear_idx - 1
                end do
            end if
            ye_smear = y_cb(smear_idx)

            if (p > 0) then
                smear_idx = cell(3) - lag_params%influence - 1
                if (smear_idx < -buff_size - 1) then
                    do while (smear_idx < -buff_size - 1)
                        smear_idx = smear_idx + 1
                    end do
                end if
                zb_smear = z_cb(smear_idx)

                smear_idx = cell(3) + lag_params%influence
                if (smear_idx > p + buff_size) then
                    do while (smear_idx > p + buff_size)
                        smear_idx = smear_idx - 1
                    end do
                end if
                ze_smear = z_cb(smear_idx)
            end if

            if (any(lag_params%interaction_model == (/2, 3/)) .and. p == 0) then

                yb_smear = mtn_posPrev(j, 2, 1) - abs(xe_smear - xb_smear)
                ye_smear = mtn_posPrev(j, 2, 1) + abs(xe_smear - xb_smear)

                zb_smear = mtn_posPrev(j, 3, 1) - abs(xe_smear - xb_smear)
                ze_smear = mtn_posPrev(j, 3, 1) + abs(xe_smear - xb_smear)

            end if

            ! Find bubbles inside the boundaries
            nb_local = 0
            if (any(lag_params%interaction_model == (/2, 3/))) then
                $:GPU_LOOP(parallelism='[seq]')
                do k = 1, nBubs
                    if ((mtn_posPrev(k, 1, 1) < xe_smear) .and. (mtn_posPrev(k, 1, 1) >= xb_smear) .and. &
                        (mtn_posPrev(k, 2, 1) < ye_smear) .and. (mtn_posPrev(k, 2, 1) >= yb_smear) .and. &
                        (mtn_posPrev(k, 3, 1) < ze_smear) .and. (mtn_posPrev(k, 3, 1) >= zb_smear)) then

                        if (k /= j) then
                            nb_local = nb_local + 1
                            if (nb_local <= max_bub_int) then
                                bub_int_ids(j, nb_local + 1) = k
                            end if
                            safeStop = max(safeStop, 1._wp*nb_local)
                        end if

                    end if
                end do
            else
                $:GPU_LOOP(parallelism='[seq]')
                do k = 1, nBubs
                    if ((mtn_pos(k, 1, 1) < xe_smear) .and. (mtn_pos(k, 1, 1) >= xb_smear) .and. &
                        (mtn_pos(k, 2, 1) < ye_smear) .and. (mtn_pos(k, 2, 1) >= yb_smear)) then

                        if (p > 0) then
                            ! if ((mtn_pos(k, 3, 1) < ze_smear) .and. (mtn_pos(k, 3, 1) >= zb_smear)) then
                            !     nb_local = nb_local + 1
                            !     if (nb_local <= max_bub_int) then
                            !         !bub_int_ids(j, nb_local + 1) = k
                            !     end if
                            !     safeStop = max(safeStop, 1._wp*nb_local)
                            ! end if
                        else
                            nb_local = nb_local + 1
                        end if
                    end if
                end do
            end if

            if (any(lag_params%interaction_model == (/2, 3/))) then
                ! Total number of interacting bubbles
                bub_int_ids(j, 1) = nb_local
            else
                ! Compute and update the mean inter-bubble distance lambda_c
                !bub_lambda_c(j) = 1._wp/(nb_local**(1._wp/3._wp))
            end if

            ! ! Populate random phases (1:num_noise)
            ! ! Should they remain constant at every calculation? or should they be replaced every t_step?
            ! i = 1
            ! do while (.true.)

            !     call random_number(num_rn1)
            !     num_rn1 = 1._wp - num_rn1
            !     call random_number(num_rn2)
            !     num_rn2 = 1._wp - num_rn2

            !     num_rn = st_dev_rn*sqrt(-2._wp*log(num_rn1))*cos(2._wp*pi*num_rn2) + mean_rn

            !     if (num_rn >= 0._wp .and. num_rn <= 2._wp*pi) then
            !         bub_rnd_phase(j, i) = num_rn
            !         if (i == num_noise) exit
            !         i = i + 1
            !     end if

            ! end do

            ! if (lag_id(j, 1) == 1) then
            !     print*, 'bub_lambda_c', bub_lambda_c(j)
            !     ! print*, 'bub_rnd_phase', bub_rnd_phase(j, 1:num_noise)
            ! end if

        end do

        if (num_procs > 1) then
            call s_mpi_allreduce_max(safeStop, tmp_val)
            safeStop = tmp_val
        end if

        if (proc_rank==0) print '("Maximum number of interacting bubbles is ", I0)', int(safeStop)

        if (safeStop > max_bub_int) then
            call s_mpi_abort('Failed getting interacting bubbles.')
        end if

        if (any(lag_params%interaction_model == (/2, 3/))) then
            $:GPU_UPDATE(host='[bub_int_ids]')
            if (lag_params%nBubs_glb < 100) then
                do j = 1, nBubs
                    if (bub_int_ids(j, 1) /= 0) then
                        print '(" (proc: ", I3, ") Bubble ", I5, " interacts with ", I5, " bubbles.")', &
                            proc_rank, &
                            j, &
                            int(bub_int_ids(j, 1))

                    end if
                end do
            end if
        end if

    end subroutine s_start_bubble_interaction

    !>  Contains the bubble dynamics subroutines.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !! @param rhs_vf Calculated change of conservative variables
        !! @param t_step Current time step
        !! @param stage Current stage in the time-stepper algorithm
    subroutine s_compute_bubble_EL_dynamics(q_prim_vf, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in) :: stage

        real(wp) :: myVapFlux
        real(wp) :: preterm1, term2, paux, pint, Romega, term1_fac, Rb
        real(wp) :: myConc_v, myR_m, mygamma_m, myPb, myMass_n, myMass_v, myPout, myInt
        real(wp) :: myR, myV, myBeta_c, myBeta_t, myR0, myPbdot, myShell, myRbuck, myMvdot
        real(wp) :: myPinf, aux1, aux2, myCson, myRho, myRrupt, myAc
        real(wp) :: myQth, myQvis, myRcell, myRmean, myKe
        real(wp) :: gamma, pi_inf, qv
        real(wp), dimension(contxe) :: myalpha_rho, myalpha
        real(wp), dimension(2) :: Re
        integer, dimension(3) :: cell

        real(wp) :: myTzPcell, myNoise_constant, myLambda_c, myloc, mydk, myPnoise, myLag_time
        real(wp), dimension(num_noise) :: myPhase

        integer :: adap_dt_stop_max, adap_dt_stop !< Fail-safe exit if max iteration count reached
        real(wp) :: dmalf, dmntait, dmBtait, dm_bub_adv_src, dm_divu !< Dummy variables for unified subgrid bubble subroutines

        integer :: i, k, l

        real(wp), dimension(1:4) :: mom_vol, mom_qvis, mom_qth_p, mom_qth_n
        real(wp) :: fVol, fxb_Rc

        integer :: total_ids, bub_idx

        logical :: flg_bub_in_cv
        real(wp) :: acPw_qvis, acPw_qth, acPW_nbubs, acPw_ke

        call nvtxStartRange("LAGRANGE-BUBBLE-DYNAMICS")

        !< BUBBLE DYNAMICS
        if (hifu_params%moments) then
            mom_vol(1:4) = 0._wp; mom_qvis(1:4) = 0._wp
            mom_qth_p(1:4) = 0._wp; mom_qth_n(1:4) = 0._wp
        end if
        if (hifu_params%power_balance) then
            acPw_qvis = 0._wp; acPw_qth = 0._wp
            acPW_nbubs = 0._wp; acPw_ke = 0._wp
        end if

        ! Subgrid p_inf model based on Maeda and Colonius (2018).
        if (lag_params%pressure_corrector) then
            call s_calculate_scattered_pressure(q_prim_vf)
            ! Calculate velocity potentials (valid for one bubble per cell)
            ! $:GPU_PARALLEL_LOOP(private='[k,cell]')
            ! do k = 1, nBubs
            !     call s_get_pinf(k, q_prim_vf, 2, paux, cell, preterm1, term2, Romega)
            !     myR0 = bub_R0(k)
            !     myR = intfc_rad(k, 2)
            !     myV = intfc_vel(k, 2)
            !     myPb = gas_p(k, 2)
            !     pint = f_cpbw_KM(myR0, myR, myV, myPb)
            !     pint = pint + 0.5_wp*myV**2._wp
            !     if (lag_params%cluster_type == 2) then
            !         bub_dphidt(k) = (paux - pint) + term2
            !         ! Accounting for the potential induced by the bubble averaged over the control volume
            !         ! Note that this is based on the incompressible flow assumption near the bubble.
            !         term1_fac = 3._wp/2._wp*(myR*(Romega**2._wp - myR**2._wp))/(Romega**3._wp - myR**3._wp)
            !         bub_dphidt(k) = bub_dphidt(k)/(1._wp - term1_fac)
            !     end if
            ! end do
        end if

        ! Radial motion
        adap_dt_stop_max = 0
        $:GPU_PARALLEL_LOOP(private='[k,i,myalpha_rho,myalpha,Re,cell,myPinf]', &
            & reduction='[[adap_dt_stop_max],[mom_vol(1:4),mom_qvis(1:4),mom_qth_p(1:4),mom_qth_n(1:4)],[acPw_qvis,acPw_qth,acPW_nbubs,acPw_ke]]', &
            & reductionOp='[MAX,+,+]', &
            & copy='[adap_dt_stop_max,mom_vol(1:4),mom_qvis(1:4),mom_qth_p(1:4),mom_qth_n(1:4),acPw_qvis,acPw_qth,acPW_nbubs,acPw_ke]', &
            & copyin='[stage]')
        do k = 1, nBubs
            ! Keller-Miksis model

            ! Current bubble state
            myPb = gas_p(k, 2)
            myMass_n = gas_mg(k)
            myMass_v = gas_mv(k, 2)
            myR = intfc_rad(k, 2)
            myV = intfc_vel(k, 2)
            myBeta_c = gas_betaC(k)
            myBeta_t = gas_betaT(k)
            myR0 = bub_R0(k)
            myShell = mrmtnt_shell(k, 2)
            myRbuck = mrmtnt_Rbuck(k)
            myRrupt = mrmtnt_Rrupt(k)
            if (myR > myRrupt) myShell = 0._wp
            myLag_time = mytime - dt
            myPout = 0._wp !Self-scaterred pressure
            if (lag_params%pressure_corrector .and. any(lag_params%interaction_model == (/1, 3/)) .and. &
                .not. adap_dt) myPout = bub_interact(k)
            myInt = 0._wp !Interaction term from surrounding bubbles
            if (lag_params%pressure_corrector .and. any(lag_params%interaction_model == (/2, 3/))) myInt = bub_interact(k)

            ! Vapor and heat fluxes
            if (.not. polytropic) then
                call s_vflux(myR, myV, myPb, myMass_v, k, myVapFlux, myMass_n, myBeta_c, myR_m, mygamma_m, myShell)
                myPbdot = f_bpres_dot(myVapFlux, myR, myV, myPb, myMass_v, k, myBeta_t, myR_m, mygamma_m, myShell)
                myMvdot = 4._wp*pi*myR**2._wp*myVapFlux
            else
                myVapFlux = 0._wp; myPbdot = 0._wp; myMvdot = 0._wp
            end if

            ! Retrieving driving pressure
            call s_get_pinf(k, q_prim_vf, 1, myPinf, cell, aux1, aux2, myRcell)

            ! Obtain liquid density and computing speed of sound from pinf
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                myalpha_rho(i) = q_prim_vf(i)%sf(cell(1), cell(2), cell(3))
                myalpha(i) = q_prim_vf(E_idx + i)%sf(cell(1), cell(2), cell(3))
            end do
            call s_convert_species_to_mixture_variables_acc(myRho, gamma, pi_inf, qv, myalpha, &
                                                            myalpha_rho, Re)
            if (lag_params%pressure_corrector .and. any(lag_params%interaction_model == (/1, 3/)) .and. &
                .not. adap_dt) then
                !Kazuki's model to adjust Pinf
                if (p > 0) then
                    myPinf = myPinf + myPout
                    ! else
                    !     ! White noise for 2D reduced model (myTzPcell is myPinf)
                    !     myLambda_c = bub_lambda_c(k)
                    !     myloc = mtn_s(k, 3, 2)
                    !     !myPhase = bub_rnd_phase(k, 1:num_noise)
                    !     call s_white_noise_constants(k, myLambda_c, q_prim_vf, cell, myPinf, myNoise_constant, mydk)
                    !     call s_compute_cson_from_pinf(q_prim_vf, myPinf, cell, myRho, gamma, pi_inf, myCson)
                    !     myPnoise = f_pres_stochastic(myPinf, myNoise_constant, myLambda_c, mydk, myloc, myLag_time, myCson)
                    !     myPinf = myPinf + myPnoise*lag_params%pnoise_scale
                end if
            end if
            call s_compute_cson_from_pinf(q_prim_vf, myPinf, cell, myRho, gamma, pi_inf, myCson)

            ! Adaptive time stepping
            if (adap_dt) then

                if (stage == 3) myLag_time = mytime - 0.5_wp*dt

                call s_advance_step(myRho, myPinf, myR, myV, myR0, myPb, myPbdot, dmalf, &
                                    dmntait, dmBtait, dm_bub_adv_src, dm_divu, &
                                    k, myMass_v, myMass_n, myBeta_c, &
                                    myBeta_t, myCson, myInt, myShell, myRbuck, myRrupt, myRcell, &
                                    myNoise_constant, myLambda_c, mydk, myloc, myLag_time, myAc, & !myPhase, &
                                    myQvis, myQth, myKe, myRmean, adap_dt_stop)

                ! Update bubble state
                intfc_rad(k, 1) = myR
                intfc_vel(k, 1) = myV
                intfc_ac(k, 1) = myAc
                gas_p(k, 1) = myPb
                if (polytropic) gas_p(k, 1) = pv + (myPb - pv)*(bub_R0(k)/myR)**(3._wp*gamma_m)
                gas_mv(k, 1) = myMass_v
                mrmtnt_shell(k, 1) = myShell
                if (hifu_params%sampling) then
                    bub_qvis(k) = bub_qvis(k) + myQvis  !> Viscous damping of the bubble (Watts*second)
                    bub_qth(k) = bub_qth(k) + myQth     !> Thermal damping of the bubble (Watts*second)
                    bub_hifu_rad(k) = bub_hifu_rad(k) + myRmean !> Mean radius (m*second)
                    ! if (k == 1) print *, 'Sampling qvis and qth (adap dt)', stage, bub_qvis(k), bub_qth(k)
                    if (hifu_params%moments) then
                        fxb_Rc = (mtn_pos(k, 1, 1)-hifu_params%cloud_center(1))/hifu_params%R_cloud
                        fVol = (4._wp/3._wp)*pi*myR**3._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 4
                            mom_vol(i) = mom_vol(i) + fVol*(fxb_Rc)**(i-1)
                            mom_qvis(i) = mom_qvis(i) + (myQvis/(0.5_wp*dt))*(fxb_Rc)**(i-1)
                            if (myQth < 0._wp) then
                                mom_qth_p(i) = mom_qth_p(i) + (myQth/(0.5_wp*dt))*(fxb_Rc)**(i-1)
                            else
                                mom_qth_n(i) = mom_qth_n(i) + (myQth*(0.5_wp*dt))*(fxb_Rc)**(i-1)
                            end if
                        end do
                    end if
                    if (hifu_params%power_balance) then
                        flg_bub_in_cv = f_bub_in_cv(mtn_pos(k, 1:3, 1))
                        if (flg_bub_in_cv) then
                            acPw_qvis = acPw_qvis + myQvis/(0.5_wp*dt)  !(Watts)
                            acPw_qth = acPw_qth + myQth/(0.5_wp*dt)     !(Watts)
                            acPW_nbubs = acPW_nbubs + 1._wp
                            acPw_ke = acPw_ke + myKe/(0.5_wp*dt)  !(Watts)
                        end if
                    end if
                end if

            else

                ! Radial acceleration from bubble models
                intfc_dveldt(k, stage) = f_rddot(myRho, myPinf, myR, myV, myR0, &
                                                 myPb, myPbdot, dmalf, dmntait, dmBtait, &
                                                 dm_bub_adv_src, dm_divu, &
                                                 myCson, myInt, myShell, myRbuck, myRcell)
                intfc_draddt(k, stage) = myV
                gas_dmvdt(k, stage) = myMvdot
                gas_dpdt(k, stage) = myPbdot
                mrmtnt_shell(k, 2) = myShell

                ! Bubble translation
                ! $:GPU_LOOP(parallelism='[seq]')
                ! do l = 1, 3
                !     mtn_dposdt(k, l, stage) = 0._wp
                !     mtn_dveldt(k, l, stage) = 0._wp
                ! end do

            end if

            adap_dt_stop_max = max(adap_dt_stop_max, adap_dt_stop)

        end do

        if (adap_dt .and. adap_dt_stop_max > 0) call s_mpi_abort("Adaptive time stepping failed to converge.")

        if (hifu_params%sampling .and. adap_dt) then
            if (hifu_params%moments) then
                if (stage == 3) then
                    do i=1,4
                        moments_bubs(1, i) = moments_bubs(1, i) + mom_qvis(i)
                        moments_bubs(2, i) = moments_bubs(2, i) + mom_qth_p(i)
                        moments_bubs(3, i) = moments_bubs(3, i) + mom_qth_n(i)
                        moments_bubs(4, i) = moments_bubs(4, i) + mom_vol(i)
                    end do

                    do i = 1, 4
                        call s_write_moments(moments_bubs(i, 1:4), idx=i)
                    end do

                else
                    do i=1,4
                        moments_bubs(1, i) = mom_qvis(i)
                        moments_bubs(2, i) = mom_qth_p(i)
                        moments_bubs(3, i) = mom_qth_n(i)
                        moments_bubs(4, i) = mom_vol(i)
                    end do
                end if
            end if
            if (hifu_params%power_balance) then
                if (stage == 3) then
                    acPw_bubs(1) = 0.5_wp * acPw_bubs(1) + 0.5_wp * acPw_qvis
                    acPw_bubs(2) = 0.5_wp * acPw_bubs(2) + 0.5_wp * acPw_qth
                    acPw_bubs(3) = 0.5_wp * acPw_bubs(3) + 0.5_wp * acPw_ke
                    call s_write_power_balance_bubs(acPw_bubs(1), acPw_bubs(2), acPw_bubs(3), acPW_nbubs, dt)
                else
                    acPw_bubs(1) = acPw_qvis
                    acPw_bubs(2) = acPw_qth
                    acPw_bubs(3) = acPw_ke
                end if
            end if
        end if

        call nvtxEndRange

    end subroutine s_compute_bubble_EL_dynamics

    function f_bub_in_cv(pos_part)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(3), intent(in) :: pos_part
        logical :: f_bub_in_cv

        f_bub_in_cv =  ((pos_part(1) < hifu_params%cv_xe) .and. (pos_part(1) >= hifu_params%cv_xb) .and. &
                        (pos_part(2) < hifu_params%cv_ye) .and. (pos_part(2) >= hifu_params%cv_yb) .and. &
                        (pos_part(3) < hifu_params%cv_ze) .and. (pos_part(3) >= hifu_params%cv_zb))

    end function f_bub_in_cv

    subroutine s_write_power_balance_bubs(acPw_qvis, acPw_qth, acPw_ke, acPW_nbubs, hdid)
        real(wp), intent(inout) :: acPw_qvis, acPw_qth, acPw_ke, acPW_nbubs
        real(wp) :: hdid
        real(wp) :: var_glb
        integer :: i

        if (num_procs > 1) then
            call s_mpi_allreduce_sum(acPw_qvis, var_glb)
            acPw_qvis = var_glb

            call s_mpi_allreduce_sum(acPw_qth, var_glb)
            acPw_qth = var_glb

            call s_mpi_allreduce_sum(acPw_ke, var_glb)
            acPw_ke = var_glb

            call s_mpi_allreduce_sum(acPW_nbubs, var_glb)
            acPW_nbubs = var_glb
        end if

        if (proc_rank == 0) then 
          
          write (89, '(*(E24.8,:,","))') &
                mytime, hdid, &
                acPW_nbubs, &
                acPw_qvis, &
                acPw_qth, &
                acPw_ke
                
        end if

    end subroutine s_write_power_balance_bubs

    !>  The purpose of this subroutine is to obtain the bubble source terms based on Maeda and Colonius (2018)
        !!      and add them to the RHS scalar field.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Conservative variables
        !! @param rhs_vf Time derivative of the conservative variables
    subroutine s_compute_bubbles_EL_source(q_cons_vf, q_prim_vf, rhs_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        integer :: i, j, k, l

        if (.not. adap_dt) call s_smear_voidfraction(bc_type)

        if (lag_params%solver_approach == 2) then

            if (p == 0 .and. .not. lag_params%newModel_2D) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            do l = 1, E_idx
                                if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                    rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                            q_cons_vf(l)%sf(i, j, k)*(q_beta%vf(2)%sf(i, j, k) + &
                                                                                      q_beta%vf(5)%sf(i, j, k))

                                end if
                            end do
                        end do
                    end do
                end do
            else
                $:GPU_PARALLEL_LOOP(collapse=4)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            do l = 1, E_idx
                                if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                    rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                            q_cons_vf(l)%sf(i, j, k)/q_beta%vf(1)%sf(i, j, k)* &
                                                            q_beta%vf(2)%sf(i, j, k)
                                end if
                            end do
                        end do
                    end do
                end do
            end if

            do l = 1, num_dims

                call s_gradient_dir(q_prim_vf(E_idx), q_beta%vf(3), l)

                $:GPU_PARALLEL_LOOP(collapse=3)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                rhs_vf(contxe + l)%sf(i, j, k) = rhs_vf(contxe + l)%sf(i, j, k) - &
                                                                 (1._wp - q_beta%vf(1)%sf(i, j, k))/ &
                                                                 q_beta%vf(1)%sf(i, j, k)* &
                                                                 q_beta%vf(3)%sf(i, j, k)
                            end if
                        end do
                    end do
                end do

                !source in energy
                $:GPU_PARALLEL_LOOP(collapse=3)
                do k = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(2)%beg, idwbuff(2)%end
                        do i = idwbuff(1)%beg, idwbuff(1)%end
                            q_beta%vf(3)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)*q_prim_vf(contxe + l)%sf(i, j, k)
                        end do
                    end do
                end do

                call s_gradient_dir(q_beta%vf(3), q_beta%vf(4), l)

                $:GPU_PARALLEL_LOOP(collapse=3)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) - &
                                                            q_beta%vf(4)%sf(i, j, k)*(1._wp - q_beta%vf(1)%sf(i, j, k))/ &
                                                            q_beta%vf(1)%sf(i, j, k)
                            end if
                        end do
                    end do
                end do
            end do

        end if

    end subroutine s_compute_bubbles_EL_source

    !>  This procedure computes the speed of sound from a given driving pressure
        !! @param bub_id Bubble id
        !! @param q_prim_vf Primitive variables
        !! @param pinf Driving pressure
        !! @param cell Bubble cell
        !! @param rhol Liquid density
        !! @param gamma Liquid specific heat ratio
        !! @param pi_inf Liquid stiffness
        !! @param cson Calculated speed of sound
    pure subroutine s_compute_cson_from_pinf(q_prim_vf, pinf, cell, rhol, gamma, pi_inf, cson)
        $:GPU_ROUTINE(function_name='s_compute_cson_from_pinf', &
            & parallelism='[seq]', cray_inline=True)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(in) :: pinf, rhol, gamma, pi_inf
        integer, dimension(3), intent(in) :: cell
        real(wp), intent(out) :: cson

        real(wp) :: E, H
        real(wp), dimension(num_dims) :: vel
        integer :: i

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_dims
            vel(i) = q_prim_vf(i + contxe)%sf(cell(1), cell(2), cell(3))
        end do
        E = gamma*pinf + pi_inf + 0.5_wp*rhol*dot_product(vel, vel)
        H = (E + pinf)/rhol
        cson = sqrt((H - 0.5_wp*dot_product(vel, vel))/gamma)

    end subroutine s_compute_cson_from_pinf

    !>  The purpose of this subroutine is to smear the effect of the bubbles in the Eulerian framework
    subroutine s_smear_voidfraction(bc_type)

        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        integer :: i, j, k, l

        call nvtxStartRange("BUBBLES-LAGRANGE-KERNELS")

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, q_beta_idx
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_beta%vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do

        if (lag_params%newModel_2D) then
            call s_smoothfunction(nBubs, intfc_rad, intfc_vel, &
                                  mtn_s, mtn_posPrev, q_beta)
        else
            call s_smoothfunction(nBubs, intfc_rad, intfc_vel, &
                                  mtn_s, mtn_pos, q_beta)
        end if

        ! Add effect of bubbles across processors
        if (num_procs > 0) call s_populate_EL_buffers(q_beta, bc_type, q_beta_idx, .false.)

        !Store 1-beta
        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    q_beta%vf(1)%sf(j, k, l) = 1._wp - q_beta%vf(1)%sf(j, k, l)
                    ! Limiting void fraction given max value
                    q_beta%vf(1)%sf(j, k, l) = max(q_beta%vf(1)%sf(j, k, l), &
                                                   1._wp - lag_params%valmaxvoid)
                end do
            end do
        end do

        call nvtxEndRange

    end subroutine s_smear_voidfraction

!     !> The purpose of this procedure is obtain the pressure from the Eulerian field
!         !! @param bub_id Particle identifier
!         !! @param q_prim_vf  Primitive variables
!     subroutine s_white_noise_constants(bub_id, l_c, q_prim_vf, cell, TzPcell, noise_constant, dk)
! #ifdef _CRAYFTN
!         !DIR$ INLINEALWAYS s_white_noise_constants
! #else
!         !$acc routine seq
! #endif
!         integer, intent(in) :: bub_id
!         real(wp), intent(in) :: l_c
!         type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
!         integer, dimension(3), intent(in) :: cell
!         real(wp), intent(out) :: TzPcell, noise_constant, dk

!         real(wp) :: vol
!         real(wp) :: denom
!         real(wp) :: charvol, charpres, charvol2, charpres2, charbeta
!         real(wp) :: charpres_sqrd, charpres2_sqrd, cell_count, c_t
!         real(wp) :: chardist, k
!         integer, dimension(3) :: cellaux
!         integer :: i, j
!         integer :: smearGrid
!         logical :: celloutside

!         if (lag_params%cluster_type == 2) then
!             ! Stochastic closure from Maeda and Colonius (2018)
!             ! Only valid for 2D reduced model

!             ! Conditions:
!             if (lag_params%smooth_type /= 1) stop "lag_params%cluster_type: 2 requires lag_params%smooth_type: 1."

!             ! Include the cell that contains the bubble (mapCells+1+mapCells)
!             ! Assume that the bubble radius is always smaller than the characteristic cell size.

!             smearGrid = mapCells - (-mapCells)

!             vol = 0._wp
!             charvol = 0._wp
!             charpres = 0._wp
!             charpres_sqrd = 0._wp
!             charvol2 = 0._wp
!             charpres2 = 0._wp
!             charpres2_sqrd = 0._wp
!             cell_count = 0._wp

!             $:GPU_LOOP(parallelism='[seq]')
!             do i = 0, smearGrid
!                 $:GPU_LOOP(parallelism='[seq]')
!                 do j = 0, smearGrid
!                     cellaux(1) = cell(1) + i - mapCells
!                     cellaux(2) = cell(2) + j - mapCells
!                     cellaux(3) = 0

!                     !< check if the current cell is outside the computational domain or not (including ghost cells)
!                     celloutside = .false.
!                     if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
!                         celloutside = .true.
!                     end if
!                     if ((cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
!                         celloutside = .true.
!                     end if
!                     if (.not. celloutside .and. cyl_coord) then
!                         if (cellaux(2) < 0) celloutside = .true.
!                     end if

!                     if (.not. celloutside) then
!                         !< Obtaining the cell volulme
!                         if (cyl_coord) then
!                             vol = dx(cellaux(1))*dy(cellaux(2))*y_cc(cellaux(2))*2._wp*pi
!                         else
!                             vol = dx(cellaux(1))*dy(cellaux(2))*lag_params%charwidth
!                         end if
!                         !< Update values
!                         cell_count = cell_count + 1._wp
!                         charvol = charvol + vol
!                         charpres = charpres + vol*q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3))
!                         charpres_sqrd = charpres_sqrd + vol*(q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)))**2._wp
!                         charvol2 = charvol2 + vol*q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
!                         charpres2 = charpres2 + vol*q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
!                                                    *q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
!                         charpres2_sqrd = charpres2_sqrd + vol*(q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3)) &
!                                                    *q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)))**2._wp
!                     end if
!                 end do
!             end do

!             TzPcell = charpres2/charvol2
!             noise_constant = 0._wp

!             if (lag_params%pressure_corrector) then
!                 ! find noise_constant: C_A

!                 if (cyl_coord) then
!                     c_t = ceiling(2._wp*pi*y_cc(cell(2))/(y_cb(cell(2)) - y_cb(cell(2) - 1)))
!                 else
!                     c_t = ceiling(lag_params%charwidth/(y_cb(cell(2)) - y_cb(cell(2) - 1)))
!                 end if
!                 c_t = c_t / cell_count

!                 chardist = sqrt(dx(cell(1))*dy(cell(2)))
!                 dk = pi/(num_noise*chardist)

!                 k = 0._wp
!                 denom = 0._wp
!                 $:GPU_LOOP(parallelism='[seq]')
!                 do i = 1, num_noise
!                     denom = denom + dk * exp(-0.5_wp*((2_wp*pi/k - l_c)/(0.5_wp*l_c))**2._wp)
!                     k = k + dk
!                 end do

!                 noise_constant = (0.5_wp*l_c) * sqrt(2._wp*pi) * c_t * ((charpres_sqrd/charvol)-(charpres/charvol)**2_wp) / denom
!                 !noise_constant = (0.5_wp*l_c) * sqrt(2._wp*pi) * c_t * ((charpres2_sqrd/charvol2)-(charpres2/charvol2)**2_wp) / denom

!                 !if (noise_constant == 0._wp) stop "noise_constant si zero. Exiting."
!                 if (lag_id(bub_id, 1) ==1) print*, noise_constant, c_t, ((charpres_sqrd/charvol)-(charpres/charvol)**2_wp), ((charpres2_sqrd/charvol2)-(charpres2/charvol2)**2_wp)

!                 !if (noise_constant <= 0._wp) noise_constant = 0._wp
!             end if

!         else

!             stop "Check white noise. Exiting."

!         end if

!     end subroutine s_white_noise_constants

    subroutine s_calculate_scattered_pressure(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        integer, dimension(3) :: cell
        real(wp) :: myR0, myR, myV, myPb, myShell, myRbuck, myRrupt, myDist, myA

        real(wp) :: Pcell, Rcell, Pw, myRho, myPout, sumPout
        real(wp) :: preterm1, term2, aux, denom, c1, c2, myInt
        integer :: bub_idx, total_ids
        integer :: i, k

        ! if (lag_params%interaction_model == 0) then !MFC for tests
        !     !$acc parallel loop gang vector default(present) private(k, cell)
        !     do k = 1, nBubs
        !         ! Current bubble state
        !         myR0 = bub_R0(k)
        !         myR = intfc_rad(k, 2)
        !         myV = intfc_vel(k, 2)
        !         myPb = gas_p(k, 2)
        !         myShell = mrmtnt_shell(k, 2)
        !         myRbuck = mrmtnt_Rbuck(k)
        !         myRrupt = mrmtnt_Rrupt(k)
        !         if (myR > myRrupt) myShell = 0._wp

        !         ! Calculate velocity potentials (valid for one bubble per cell)
        !         call s_get_pinf(k, q_prim_vf, 2, Pcell, cell, preterm1, term2, Rcell)

        !         ! Obtain liquid density and computing speed of sound from myPinf
        !         myRho = 0._wp
        !         $:GPU_LOOP(parallelism='[seq]')
        !         do i = 1, num_fluids
        !             myRho = myRho + q_prim_vf(i)%sf(cell(1), cell(2), cell(3))
        !         end do

        !         aux = Rcell**3._wp - myR**3._wp
        !         c2 = (3._wp/2._wp)*(myR**3._wp)*(1._wp - myR/Rcell)/aux
        !         c1 = 3._wp/2._wp*(myR*(Rcell**2._wp - myR**2._wp))/aux

        !         Pw = f_cpbw_KM(myR0, myR, myV, myPb, myShell, myRbuck)
        !         Pw = Pw + 0.5_wp*myV**2._wp
        !         bub_dphidt(k) = (Pcell - Pw) + c2*myV**2._wp
        !         ! Accounting for the potential induced by the bubble averaged over the control volume
        !         ! Note that this is based on the incompressible flow assumption near the bubble.
        !         bub_dphidt(k) = bub_dphidt(k)/(1._wp - c1)

        !         ! Scattered pressure
        !         myPout = c1*bub_dphidt(k) + c2*myV**2._wp

        !         !Update emitted Pout
        !         bub_interact(k) = myPout
        !         print*, 'myPout matching:', myPout, bub_dphidt(k), mytime

        !     end do

        ! end if

        if (any(lag_params%interaction_model == (/1, 3/)) .and. p > 0 .and. .not. adap_dt) then !Kazuki's model (DV version)
            $:GPU_PARALLEL_LOOP(private='[k, cell]')
            do k = 1, nBubs
                ! Current bubble state
                myR0 = bub_R0(k)
                myR = intfc_rad(k, 2)
                myV = intfc_vel(k, 2)
                myPb = gas_p(k, 2)
                myShell = mrmtnt_shell(k, 2)
                myRbuck = mrmtnt_Rbuck(k)
                myRrupt = mrmtnt_Rrupt(k)
                if (myR > myRrupt) myShell = 0._wp

                ! ! Calculate velocity potentials (valid for one bubble per cell)
                ! call s_get_pinf(k, q_prim_vf, 2, Pcell, cell, preterm1, term2, Rcell)

                ! ! Obtain liquid density and computing speed of sound from myPinf
                ! myRho = 0._wp
                ! $:GPU_LOOP(parallelism='[seq]')
                ! do i = 1, num_fluids
                !     myRho = myRho + q_prim_vf(i)%sf(cell(1), cell(2), cell(3))
                ! end do

                ! aux = Rcell**3._wp - myR**3._wp
                ! c2 = (3._wp/2._wp)*(myR**3._wp)*(1._wp - myR/Rcell)/aux
                ! c1 = (3._wp/2._wp)*(myR*(Rcell**2._wp - myR**2._wp))/aux

                ! Pw = f_cpbw_KM(myR0, myR, myV, myPb, myShell, myRbuck)
                ! Pw = Pw/myRho - 0.5_wp*myV**2._wp
                ! bub_dphidt(k) = (Pcell/myRho - Pw) - c2*myV**2._wp
                ! ! Accounting for the potential induced by the bubble averaged over the control volume
                ! ! Note that this is based on the incompressible flow assumption near the bubble.
                ! bub_dphidt(k) = bub_dphidt(k)/(1._wp - c1)

                ! ! Scattered pressure
                ! myPout = myRho*(c1*bub_dphidt(k) - c2*myV**2._wp)

                ! !Update emitted Pout
                ! bub_interact(k) = myPout
                ! !print*, 'myPout matching:', myPout

                ! Calculate velocity potentials (valid for one bubble per cell)
                call s_get_pinf(k, q_prim_vf, 2, Pcell, cell, preterm1, term2, Rcell)

                ! Obtain liquid density and computing speed of sound from myPinf
                myRho = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    myRho = myRho + q_prim_vf(i)%sf(cell(1), cell(2), cell(3))
                end do

                aux = Rcell**3._wp - myR**3._wp
                c2 = (3._wp/2._wp)*(myR**3._wp)*(1._wp - myR/Rcell)/aux
                c1 = 3._wp/2._wp*(myR*(Rcell**2._wp - myR**2._wp))/aux

                Pw = f_cpbw_KM(myR0, myR, myV, myPb, myShell, myRbuck)
                Pw = Pw + 0.5_wp*myV**2._wp
                bub_dphidt(k) = (Pcell - Pw) + c2*myV**2._wp
                ! Accounting for the potential induced by the bubble averaged over the control volume
                ! Note that this is based on the incompressible flow assumption near the bubble.
                bub_dphidt(k) = bub_dphidt(k)/(1._wp - c1)

                ! Scattered pressure
                myPout = c1*bub_dphidt(k) + c2*myV**2._wp

                !Update emitted Pout
                bub_interact(k) = myPout
            end do

        end if

        if (any(lag_params%interaction_model == (/2, 3/))) then !Aditya's model, Pout is going to be I term from eqn 3.19

            $:GPU_PARALLEL_LOOP(private='[k, cell]')
            do k = 1, nBubs

                ! Number of the bubbles in the smearing volume (Self-inclusive)
                total_ids = int(bub_int_ids(k, 1))
                sumPout = 0._wp

                if (total_ids + 1 >= 2) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 2, total_ids + 1
                        bub_idx = bub_int_ids(k, i)
                        ! Current interacting bubble state
                        myR = intfc_rad(bub_idx, 2)
                        myV = intfc_vel(bub_idx, 2)
                        myA = intfc_ac(bub_idx, 2)
                        myDist = (mtn_posPrev(bub_idx, 1, 2) - mtn_posPrev(k, 1, 2))**2._wp + &
                             (mtn_posPrev(bub_idx, 2, 2) - mtn_posPrev(k, 2, 2))**2._wp + &
                             (mtn_posPrev(bub_idx, 3, 2) - mtn_posPrev(k, 3, 2))**2._wp
                        myDist = sqrt(myDist)

                        !if (bub_idx /= k .and. .not. f_approx_equal(myDist, 0._wp)) then  ! non-inclusive for Aditya's model
                         myInt = (2._wp*myR*myV**2._wp + myA*myR**2._wp)/myDist
                        if (myInt /= myInt) then
                            print*, myR, myV, myA, myDist, 'Bub', k, 'with bub', bub_idx
                        else
                            sumPout = sumPout - myInt
                        end if

                        !if (k==50) print*, 'bub-50:', k, bub_int_ids(k, 1), bub_idx, sumPout, myR, myV, myA, myDist
                        !if (k==50) print*, 'bub-50: bub_int_ids', bub_int_ids(k, i)
                    end do
                end if

                ! I term: sum over bubbles
                bub_interact(k) = sumPout
                !print*, 'sum Pout:', bub_interact(k)

            end do

        end if

        !call s_mpi_barrier()
    end subroutine s_calculate_scattered_pressure

    !> The purpose of this procedure is obtain the bubble driving pressure p_inf (OLD VERSION)
        !! @param bub_id Particle identifier
        !! @param q_prim_vf  Primitive variables
        !! @param ptype 1: p at infinity, 2: averaged P at the bubble location
        !! @param f_pinfl Driving pressure
        !! @param cell Bubble cell
        !! @param Romega Control volume radius
    pure subroutine s_get_pinf(bub_id, q_prim_vf, ptype, f_pinfl, cell, preterm1, term2, Romega)
        $:GPU_ROUTINE(function_name='s_get_pinf',parallelism='[seq]', &
            & cray_inline=True)

        integer, intent(in) :: bub_id, ptype
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(out) :: f_pinfl
        integer, dimension(3), intent(out) :: cell
        real(wp), intent(out), optional :: preterm1, term2, Romega

        real(wp), dimension(3) :: scoord, psi
        real(wp) :: dc, vol, aux, chardist, dist_cc
        real(wp) :: volgas, term1, Rbeq, denom
        real(wp) :: charvol, charpres, charvol2, charpres2, charbeta
        integer, dimension(3) :: cellaux
        integer :: i, j, k
        integer :: mapCells_pinf, smearGrid, smearGridz
        logical :: celloutside, condition

        scoord = mtn_s(bub_id, 1:3, 2)
        f_pinfl = 0._wp

        !< Find current bubble cell
        cell(:) = int(scoord(:))
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_dims
            if (scoord(i) < 0._wp) cell(i) = cell(i) - 1
        end do

        if ((lag_params%cluster_type == 1)) then
            !< Getting p_cell in terms of only the current cell by interpolation

            !< Getting the cell volulme as Omega
            if (p > 0) then
                vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            else
                if (cyl_coord) then
                    vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
                else
                    vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                end if
            end if

            !< Obtain bilinear interpolation coefficients, based on the current location of the bubble.
            psi(1) = (scoord(1) - real(cell(1)))*dx(cell(1)) + x_cb(cell(1) - 1)
            if (cell(1) == (m + buff_size)) then
                cell(1) = cell(1) - 1
                psi(1) = 1._wp
            else if (cell(1) == (-buff_size)) then
                psi(1) = 0._wp
            else
                if (psi(1) < x_cc(cell(1))) cell(1) = cell(1) - 1
                psi(1) = abs((psi(1) - x_cc(cell(1)))/(x_cc(cell(1) + 1) - x_cc(cell(1))))
            end if

            psi(2) = (scoord(2) - real(cell(2)))*dy(cell(2)) + y_cb(cell(2) - 1)
            if (cell(2) == (n + buff_size)) then
                cell(2) = cell(2) - 1
                psi(2) = 1._wp
            else if (cell(2) == (-buff_size)) then
                psi(2) = 0._wp
            else
                if (psi(2) < y_cc(cell(2))) cell(2) = cell(2) - 1
                psi(2) = abs((psi(2) - y_cc(cell(2)))/(y_cc(cell(2) + 1) - y_cc(cell(2))))
            end if

            if (p > 0) then
                psi(3) = (scoord(3) - real(cell(3)))*dz(cell(3)) + z_cb(cell(3) - 1)
                if (cell(3) == (p + buff_size)) then
                    cell(3) = cell(3) - 1
                    psi(3) = 1._wp
                else if (cell(3) == (-buff_size)) then
                    psi(3) = 0._wp
                else
                    if (psi(3) < z_cc(cell(3))) cell(3) = cell(3) - 1
                    psi(3) = abs((psi(3) - z_cc(cell(3)))/(z_cc(cell(3) + 1) - z_cc(cell(3))))
                end if
            else
                psi(3) = 0._wp
            end if

            !< Perform bilinear interpolation
            if (p == 0) then  !2D
                f_pinfl = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))*(1._wp - psi(1))*(1._wp - psi(2))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1._wp - psi(2))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3))*(1._wp - psi(1))*psi(2)
            else              !3D
                f_pinfl = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))*(1._wp - psi(1))*(1._wp - psi(2))*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1._wp - psi(2))*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3))*(1._wp - psi(1))*psi(2)*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) + 1)*(1._wp - psi(1))*(1._wp - psi(2))*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3) + 1)*psi(1)*(1._wp - psi(2))*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3) + 1)*psi(1)*psi(2)*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3) + 1)*(1._wp - psi(1))*psi(2)*psi(3)
            end if

            !R_Omega
            dc = (3._wp*vol/(4._wp*pi))**(1._wp/3._wp)

        else if (lag_params%cluster_type == 2) then
            ! Bubble dynamic closure from Maeda and Colonius (2018)

            ! Include the cell that contains the bubble (mapCells+1+mapCells)
            smearGrid = mapCells - (-mapCells) + 1
            smearGridz = smearGrid
            if (p == 0) smearGridz = 1

            charvol = 0._wp
            charpres = 0._wp
            charvol2 = 0._wp
            charpres2 = 0._wp
            vol = 0._wp
            charbeta = 0._wp

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, smearGrid
                $:GPU_LOOP(parallelism='[seq]')
                do j = 1, smearGrid
                    $:GPU_LOOP(parallelism='[seq]')
                    do k = 1, smearGridz
                        cellaux(1) = cell(1) + i - (mapCells + 1)
                        cellaux(2) = cell(2) + j - (mapCells + 1)
                        cellaux(3) = cell(3) + k - (mapCells + 1)
                        if (p == 0) cellaux(3) = 0

                        !< check if the current cell is outside the computational domain or not (including ghost cells)
                        celloutside = .false.
                        if (num_dims == 2) then
                            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                celloutside = .true.
                            end if
                            if ((cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                                celloutside = .true.
                            end if
                        else
                            if ((cellaux(3) < -buff_size) .or. (cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                celloutside = .true.
                            end if

                            if ((cellaux(3) > p + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                                celloutside = .true.
                            end if
                        end if

                        if (lag_params%interaction_model == 2 .and. .not. celloutside) then
                            ! Liquid pressure from the cells around a virtual sphere of radius K*chardist that surrounds the bubble
                            if (p > 0) then
                                chardist = sqrt(dx(cell(1))*dy(cell(2))*dz(cell(3)))
                                dist_cc = sqrt((x_cc(cell(1)) - x_cc(cellaux(1)))**2._wp + &
                                               (y_cc(cell(2)) - y_cc(cellaux(2)))**2._wp + &
                                               (z_cc(cell(3)) - z_cc(cellaux(3)))**2._wp)
                                ! condition = abs(dist_cc-lag_params%scaleVirtualSphere*chardist) < 1.8_wp*chardist
                                condition = (cellaux(1) == cell(1) - mapCells .or. cellaux(1) == cell(1) + mapCells .or. &
                                             cellaux(2) == cell(2) - mapCells .or. cellaux(2) == cell(2) + mapCells .or. &
                                             cellaux(3) == cell(3) - mapCells .or. cellaux(3) == cell(3) + mapCells)
                                if (.not. condition) celloutside = .true.
                            else
                                chardist = sqrt(dx(cell(1))*dy(cell(2)))
                                dist_cc = sqrt((x_cc(cell(1)) - x_cc(cellaux(1)))**2._wp + &
                                               (y_cc(cell(2)) - y_cc(cellaux(2)))**2._wp)
                                ! condition = abs(dist_cc-lag_params%scaleVirtualSphere*chardist) < 1.5_wp*chardist
                                condition = (cellaux(1) == cell(1) - mapCells .or. cellaux(1) == cell(1) + mapCells .or. &
                                             cellaux(2) == cell(2) - mapCells .or. cellaux(2) == cell(2) + mapCells)
                                if (.not. condition) celloutside = .true.
                            end if

                        end if

                        if (.not. celloutside) then
                            !< Obtaining the cell volulme
                            if (p > 0) then
                                vol = dx(cellaux(1))*dy(cellaux(2))*dz(cellaux(3))
                            else
                                if (cyl_coord) then
                                    vol = dx(cellaux(1))*dy(cellaux(2))*y_cc(cellaux(2))*2._wp*pi
                                else
                                    vol = dx(cellaux(1))*dy(cellaux(2))*lag_params%charwidth
                                end if
                            end if

                            !< Update values
                            charvol = charvol + vol
                            charbeta = charbeta + q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                            charpres = charpres + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3))*vol
                            charvol2 = charvol2 + vol*q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                            charpres2 = charpres2 + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                        *vol*q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                            ! print*, vol, i, j, k

                        end if

                    end do
                end do
            end do

            f_pinfl = charpres2/charvol2
            if (lag_params%interaction_model == 2) f_pinfl = charpres/charvol
            vol = charvol
            dc = (3._wp*abs(vol)/(4._wp*pi))**(1._wp/3._wp)
        end if

        !Control volume radius
        Romega = dc
        ! print*, 'Rcell', Romega
        ! call s_mpi_abort('Debbuging DV')

        ! if (lag_params%pressure_corrector .and. p > 0) then

        !     !Valid if only one bubble exists per cell
        !     volgas = intfc_rad(bub_id, 2)**3._wp
        !     denom = intfc_rad(bub_id, 2)**2._wp
        !     term1 = bub_dphidt(bub_id)*intfc_rad(bub_id, 2)**2._wp
        !     term2 = intfc_vel(bub_id, 2)*intfc_rad(bub_id, 2)**2._wp

        !     Rbeq = volgas**(1._wp/3._wp) !surrogate bubble radius
        !     aux = dc**3._wp - Rbeq**3._wp
        !     term2 = term2/denom
        !     term2 = 3._wp/2._wp*term2**2._wp*Rbeq**3._wp*(1._wp - Rbeq/dc)/aux
        !     preterm1 = 3._wp/2._wp*Rbeq*(dc**2._wp - Rbeq**2._wp)/(aux*denom)

        !     !Control volume radius
        !     if (ptype == 2) Romega = dc

        !     ! Getting p_inf
        !     if (ptype == 1) then
        !         print*, 'P_inf subr:', preterm1*term1 + term2, f_pinfl + preterm1*term1 + term2
        !         f_pinfl = f_pinfl + preterm1*term1 + term2
        !     end if

        ! end if

    end subroutine s_get_pinf

    !> The purpose of this procedure is to calculate and store the time-averaged heat sources from the lagrange bubbles.
        !!      The heat sources model the viscous and thermal damping of the bubbles valid with the hifu solver.
    subroutine s_compute_bubble_heat_sources_HIFU(hdid)

        real(wp), intent(in) :: hdid

        real(wp) :: fpb_h, fmass_n_h, fmass_v_h, fR_h, fV_h, fbeta_t_h, fshell_h
        real(wp) :: conc_v_h, R_m_h, gamma_m_h, T_bar_h, grad_T_h, heatflux_h, fR0_h
        integer :: k, i
        integer :: abortFlag, abortFlag_max

        real(wp), dimension(1:4) :: mom_vol, mom_qvis, mom_qth_p, mom_qth_n
        real(wp) :: fxb_Rc, fqvis, fqth, fVol

        logical :: flg_bub_in_cv
        real(wp) :: acPw_qvis, acPw_qth, acPW_nbubs, acPw_ke
        
        if (hifu_params%moments) then
            mom_vol(1:4) = 0._wp; mom_qvis(1:4) = 0._wp
            mom_qth_p(1:4) = 0._wp; mom_qth_n(1:4) = 0._wp
        end if

        if (hifu_params%power_balance) then
            acPw_qvis = 0._wp; acPw_qth = 0._wp
            acPW_nbubs = 0._wp; acPw_ke = 0._wp
        end if


#ifdef MFC_DEBUG
        if (proc_rank == 0) print *, 'Computing bubble heat sources', mytime, hdid
#endif
        abortFlag_max = 0
        $:GPU_PARALLEL_LOOP(private='[k]',reduction='[[abortFlag_max],[acPw_qvis,acPw_qth,acPW_nbubs],[mom_vol(1:4),mom_qvis(1:4),mom_qth_p(1:4),mom_qth_n(1:4)]]', &
        & reductionOp='[MAX,+,+]',copy='[abortFlag_max, mom_vol(1:4),mom_qvis(1:4),mom_qth_p(1:4),mom_qth_n(1:4), acPw_qvis, acPw_qth, acPW_nbubs]')
        do k = 1, nBubs

            abortFlag = 0
            !> Current bubble state (no temporal values)
            fpb_h = gas_p(k, 1)
            fmass_n_h = gas_mg(k)
            fmass_v_h = gas_mv(k, 1)
            fR0_h = bub_R0(k)
            fR_h = intfc_rad(k, 1)
            fV_h = intfc_vel(k, 1)
            fbeta_t_h = gas_betaT(k)
            fshell_h = mrmtnt_shell(k, 1)
            if (hifu_params%moments) fxb_Rc = (mtn_pos(k, 1, 1)-hifu_params%cloud_center(1))/hifu_params%R_cloud

            ! Mixture properties in the bubble
            conc_v_h = 0._wp
            if (lag_params%massTransfer_model .and. (fshell_h == 0._wp)) then
                conc_v_h = 1._wp/(1._wp + (R_v/R_n)*(fpb_h/pv - 1._wp))
            end if
            R_m_h = fmass_n_h*R_n + fmass_v_h*R_v
            gamma_m_h = conc_v_h*gamma_v + (1._wp - conc_v_h)*gamma_n

            !> Viscous damping of the bubble (Watts)
            fqvis = (4._wp*pi*fR_h**2._wp)*(4._wp*mul0*(fV_h**2._wp)/(fR_h))
            bub_qvis(k) = bub_qvis(k) + hdid*fqvis

            !> Thermal damping of the bubble (Watts)
            if (.not. polytropic) then
                T_bar_h = fpb_h*(4._wp/3._wp*pi*fR_h**3._wp)/R_m_h
                grad_T_h = -fbeta_t_h*(T_bar_h - Tw)
                if (lag_params%heatTransfer_model .and. (fshell_h == 0._wp)) then
                    heatflux_h = (gamma_m_h - 1._wp)/gamma_m_h*grad_T_h/fR_h
                end if
            else
                T_bar_h = Tw * (fR0_h/fR_h)**(3._wp*(gamma_m-1._wp)) ! Polytropic temp
                heatflux_h = conc_v_h*k_vl + (1._wp - conc_v_h)*k_nl 
                heatflux_h = 3._wp*(1._wp-gamma_m)*T_bar_h/fR_h
            end if
            fqth = heatflux_h*4._wp*pi*fR_h**2._wp
            bub_qth(k) = bub_qth(k) + hdid*fqth

            !Mean radius
            bub_hifu_rad(k) = bub_hifu_rad(k) + hdid * fR_h

            fVol = (4._wp/3._wp)*pi*fR_h**3._wp

            ! Checking for NaNs and negative qvis
            if (bub_qvis(k) /= bub_qvis(k) .or. &
                bub_qth(k) /= bub_qth(k) .or. &
                bub_qvis(k) < 0._wp) then
                print *, 'Bubble intensity is NaN', k, bub_qvis(k), bub_qth(k), hdid
                print *, 'Viscous damping', fR_h, mul0, fV_h
                print *, 'Thermal damping', heatflux_h, fR_h
                abortFlag = 1
            end if

            abortFlag_max = max(abortFlag_max, abortFlag)

            if (hifu_params%moments) then
                fxb_Rc = (mtn_pos(k, 1, 1)-hifu_params%cloud_center(1))/hifu_params%R_cloud

                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, 4
                    mom_vol(i) = mom_vol(i) + fVol*(fxb_Rc)**(i-1)
                    mom_qvis(i) = mom_qvis(i) + (fqvis/hdid)*(fxb_Rc)**(i-1)
                    if (fqth < 0._wp) then
                        mom_qth_p(i) = mom_qth_p(i) + fqth*(fxb_Rc)**(i-1)
                    else
                        mom_qth_n(i) = mom_qth_n(i) + fqth*(fxb_Rc)**(i-1)
                    end if
                end do
            end if
            if (hifu_params%power_balance) then
                flg_bub_in_cv = f_bub_in_cv(mtn_pos(k, 1:3, 1))
                if (flg_bub_in_cv) then
                    acPw_qvis = acPw_qvis + fqvis !(Watts)
                    acPw_qth = acPw_qth + fqth    !(Watts)
                    acPW_nbubs = acPW_nbubs + 1._wp
                    acPw_ke = 0._wp ! Kinetic energy of the bubble, can be added if needed
                end if
            end if

        end do

        if (abortFlag_max > 0) stop "NaNs in viscous (or thermal) damping of the bubbles"

        if (hifu_params%moments) then
            call s_write_moments(mom_qvis, idx=1)
            call s_write_moments(mom_qth_p, idx=2)
            call s_write_moments(mom_qth_n, idx=3)
            call s_write_moments(mom_vol, idx=4)
        end if

        if (hifu_params%power_balance) call s_write_power_balance_bubs(acPw_qvis, acPw_qth, acPw_ke, acPW_nbubs, dt)

    end subroutine s_compute_bubble_heat_sources_HIFU

    subroutine s_write_moments(mom_all, idx)

        integer, intent(in) :: idx
        real(wp), dimension(4), intent(in) :: mom_all
        real(wp) :: total, moment1, moment2, moment3

        real(wp) :: val_tmp

        total = mom_all(1); moment1 = mom_all(2); 
        moment2 = mom_all(3); moment3 = mom_all(4)

        if (num_procs>1) then
            val_tmp = total
            call s_mpi_allreduce_sum(val_tmp, total)
            val_tmp = moment1
            call s_mpi_allreduce_sum(val_tmp, moment1)
            val_tmp = moment2
            call s_mpi_allreduce_sum(val_tmp, moment2)
            val_tmp = moment3
            call s_mpi_allreduce_sum(val_tmp, moment3)
        end if

        ! Write the heat statistics to file
        if (proc_rank == 0) then
            write (97-idx, '(4X,5e24.8)') &
                    mytime, &
                    moment1/total, &
                    moment2/total, &
                    moment3/total, &
                    total
        end if      
        
    end subroutine s_write_moments

    subroutine s_mean_radius_hifu(t_sampled)

        real(wp), intent(in) :: t_sampled
        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]',copyin='[t_sampled]')
        do k = 1, nBubs
            bub_hifu_rad(k) = bub_hifu_rad(k) / t_sampled ! meters
        end do

    end subroutine

    ! Compute the first, second, and third moments of the heat source from the bubbles' damping.
    subroutine s_write_heat_stats_bubbles(sampledTime)

        real(wp), intent(in) :: sampledTime

        real(wp) :: total_heat_vis, heat_moment1_vis, heat_moment2_vis, heat_moment3_vis
        real(wp) :: total_heat_th, heat_moment1_th, heat_moment2_th, heat_moment3_th
        real(wp) :: total_vol, moment1_vol, moment2_vol, moment3_vol
        real(wp) :: val_tmp, fR_h, fqvis, fqth, fxb_Rc, fVol
        integer :: i, j, k, l

        logical :: file_exist
        character(LEN=path_len + 2*name_len) :: file_loc

        total_heat_vis = 0._wp;   heat_moment1_vis = 0._wp
        heat_moment2_vis = 0._wp; heat_moment3_vis = 0._wp

        total_heat_th = 0._wp;    heat_moment1_th = 0._wp
        heat_moment2_th = 0._wp;  heat_moment3_th = 0._wp

        total_vol = 0._wp;    moment1_vol = 0._wp
        moment2_vol = 0._wp;  moment3_vol = 0._wp

        $:GPU_PARALLEL_LOOP(private='[k]', &
        & reduction='[[total_heat_vis, heat_moment1_vis, heat_moment2_vis, heat_moment3_vis, total_heat_th, heat_moment1_th, heat_moment2_th, heat_moment3_th, total_vol, moment1_vol, moment2_vol, moment3_vol]]', &
        & reductionOp='[MAX]', &
        & copy='[total_heat_vis, heat_moment1_vis, heat_moment2_vis, heat_moment3_vis, total_heat_th, heat_moment1_th, heat_moment2_th, heat_moment3_th, total_vol, moment1_vol, moment2_vol, moment3_vol]')
        do k = 1, nBubs

            fR_h = intfc_rad(k, 1)
            fqvis = bub_qvis(k)
            fqth = bub_qth(k)
            fxb_Rc = (mtn_pos(k, 1, 1)-hifu_params%cloud_center(1))/hifu_params%R_cloud
            fVol = (4._wp/3._wp)*pi*fR_h**3._wp

            total_heat_vis = total_heat_vis + fqvis
            heat_moment1_vis = heat_moment1_vis + fqvis*(fxb_Rc)
            heat_moment2_vis = heat_moment2_vis + fqvis*(fxb_Rc)**2._wp
            heat_moment3_vis = heat_moment3_vis + fqvis*(fxb_Rc)**3._wp

            total_heat_th = total_heat_th + fqth
            heat_moment1_th = heat_moment1_th + fqth*(fxb_Rc)
            heat_moment2_th = heat_moment2_th + fqth*(fxb_Rc)**2._wp
            heat_moment3_th = heat_moment3_th + fqth*(fxb_Rc)**3._wp

            total_vol = total_vol + fVol
            moment1_vol = moment1_vol + fVol*(fxb_Rc)
            moment2_vol = moment2_vol + fVol*(fxb_Rc)**2._wp
            moment3_vol = moment3_vol + fVol*(fxb_Rc)**3._wp

        end do

        if (num_procs>1) then
            val_tmp = total_heat_vis
            call s_mpi_allreduce_sum(val_tmp, total_heat_vis)
            val_tmp = heat_moment1_vis
            call s_mpi_allreduce_sum(val_tmp, heat_moment1_vis)
            val_tmp = heat_moment2_vis
            call s_mpi_allreduce_sum(val_tmp, heat_moment2_vis)
            val_tmp = heat_moment3_vis
            call s_mpi_allreduce_sum(val_tmp, heat_moment3_vis)

            val_tmp = total_heat_th
            call s_mpi_allreduce_sum(val_tmp, total_heat_th)
            val_tmp = heat_moment1_th
            call s_mpi_allreduce_sum(val_tmp, heat_moment1_th)
            val_tmp = heat_moment2_th
            call s_mpi_allreduce_sum(val_tmp, heat_moment2_th)
            val_tmp = heat_moment3_th
            call s_mpi_allreduce_sum(val_tmp, heat_moment3_th)

            val_tmp = total_vol
            call s_mpi_allreduce_sum(val_tmp, total_vol)
            val_tmp = moment1_vol
            call s_mpi_allreduce_sum(val_tmp, moment1_vol)
            val_tmp = moment2_vol
            call s_mpi_allreduce_sum(val_tmp, moment2_vol)
            val_tmp = moment3_vol
            call s_mpi_allreduce_sum(val_tmp, moment3_vol)

        end if
        
        ! Write the heat statistics to file
        if (proc_rank == 0) then

            write (file_loc, '(A,I0,A)') 'moments_qvis.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

             if (.not. file_exist) then
                open (11, FILE=trim(file_loc), FORM='formatted', position='rewind')
                write (11, *) 'sampledTime, normMomment_1, normMomment_2, normMomment_3, totalHeat_qvis'
            else
                open (11, FILE=trim(file_loc), FORM='formatted', position='append')
            end if

            write (11, '(4X,I24.8,4e24.8)') &
                sampledTime, &
                heat_moment1_vis/total_heat_vis, &
                heat_moment2_vis/total_heat_vis, &
                heat_moment3_vis/total_heat_vis, &
                total_heat_vis

            close(11)

            write (file_loc, '(A,I0,A)') 'moments_qth.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

             if (.not. file_exist) then
                open (11, FILE=trim(file_loc), FORM='formatted', position='rewind')
                write (11, *) 'sampledTime, normMomment_1, normMomment_2, normMomment_3, totalHeat_qth'
            else
                open (11, FILE=trim(file_loc), FORM='formatted', position='append')
            end if

            write (11, '(4X,I24.8,4e24.8)') &
                sampledTime, &
                heat_moment1_th/total_heat_th, &
                heat_moment2_th/total_heat_th, &
                heat_moment3_th/total_heat_th, &
                total_heat_th

            close(11)

            write (file_loc, '(A,I0,A)') 'moments_vol.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

             if (.not. file_exist) then
                open (11, FILE=trim(file_loc), FORM='formatted', position='rewind')
                write (11, *) 'sampledTime, normMomment_1, normMomment_2, normMomment_3, totalVolume'
            else
                open (11, FILE=trim(file_loc), FORM='formatted', position='append')
            end if

            write (11, '(4X,I24.8,4e24.8)') &
                sampledTime, &
                moment1_vol/total_vol, &
                moment2_vol/total_vol, &
                moment3_vol/total_vol, &
                total_vol

            close(11)
        end if

    end subroutine s_write_heat_stats_bubbles

    !>  This subroutine updates the Lagrange variables using the tvd RK time steppers.
        !!      The time derivative of the bubble variables must be stored at every stage to avoid precision errors.
        !! @param stage Current tvd RK stage
    impure subroutine s_update_lagrange_tdv_rk(stage)

        integer, intent(in) :: stage

        integer :: k

        if (time_stepper == 1) then ! 1st order TVD RK
            $:GPU_PARALLEL_LOOP(private='[k]')
            do k = 1, nBubs
                !u{1} = u{n} +  dt * RHS{n}
                intfc_rad(k, 1) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                intfc_vel(k, 1) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                ! mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                ! mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                gas_p(k, 1) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                gas_mv(k, 1) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                mrmtnt_shell(k, 1) = mrmtnt_shell(k, 2)
                intfc_ac(k, 1) = intfc_dveldt(k, 1)
                if (polytropic) gas_p(k, 1) = pv + (gas_p(k, 2) - pv)*(bub_R0(k)/intfc_rad(k, 1))**(3._wp*gamma_m)
            end do

            call s_transfer_data_to_tmp
            call s_calculate_lag_bubble_stats()
            if (lag_params%write_bubbles) then
                $:GPU_UPDATE(host='[gas_p,gas_mv,intfc_rad,intfc_vel]')
                call s_write_lag_particles(mytime, replace=.false.)
            end if
            call s_write_void_evol(mytime, replace=.false.)

        elseif (time_stepper == 2) then ! 2nd order TVD RK
            if (stage == 1) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{1} = u{n} +  dt * RHS{n}
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    ! mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                    ! mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    if (.not. polytropic) gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                end do

            elseif (stage == 2) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{1} = u{n} + (1/2) * dt * (RHS{n} + RHS{1})
                    intfc_rad(k, 1) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/2._wp
                    intfc_vel(k, 1) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/2._wp
                    ! mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/2._wp
                    ! mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/2._wp
                    gas_p(k, 1) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/2._wp
                    gas_mv(k, 1) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/2._wp
                    if (lag_params%coatedBub_model .and. (mrmtnt_shell(k, 2) == 0._wp)) then
                        if (intfc_rad(k, 1) < mrmtnt_Rrupt(k)) mrmtnt_shell(k, 2) = 1._wp ! No actual rupture happened during dt
                    end if
                    mrmtnt_shell(k, 1) = mrmtnt_shell(k, 2)
                    intfc_ac(k, 1) = (intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/2._wp
                    if (polytropic) gas_p(k, 1) = pv + (gas_p(k, 2) - pv)*(bub_R0(k)/intfc_rad(k, 1))**(3._wp*gamma_m)
                end do

                call s_transfer_data_to_tmp
                call s_calculate_lag_bubble_stats()
                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[gas_p,gas_mv,intfc_rad,intfc_vel]')
                    call s_write_lag_particles(mytime, replace=.false.)
                end if
                call s_write_void_evol(mytime, replace=.false.)

            end if

        elseif (time_stepper == 3) then ! 3rd order TVD RK
            if (stage == 1) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{1} = u{n} +  dt * RHS{n}
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    ! mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                    ! mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    if (.not. polytropic) gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                end do

            elseif (stage == 2) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{2} = u{n} + (1/4) * dt * [RHS{n} + RHS{1}]
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/4._wp
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/4._wp
                    ! mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/4._wp
                    ! mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/4._wp
                    if (.not. polytropic) gas_p(k, 2) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/4._wp
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/4._wp
                end do
            elseif (stage == 3) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{n+1} = u{n} + (2/3) * dt * [(1/4)* RHS{n} + (1/4)* RHS{1} + RHS{2}]
                    intfc_rad(k, 1) = intfc_rad(k, 1) + (2._wp/3._wp)*dt*(intfc_draddt(k, 1)/4._wp + intfc_draddt(k, 2)/4._wp + intfc_draddt(k, 3))
                    intfc_vel(k, 1) = intfc_vel(k, 1) + (2._wp/3._wp)*dt*(intfc_dveldt(k, 1)/4._wp + intfc_dveldt(k, 2)/4._wp + intfc_dveldt(k, 3))
                    ! mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dposdt(k, 1:3, 1)/4._wp + mtn_dposdt(k, 1:3, 2)/4._wp + mtn_dposdt(k, 1:3, 3))
                    ! mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dveldt(k, 1:3, 1)/4._wp + mtn_dveldt(k, 1:3, 2)/4._wp + mtn_dveldt(k, 1:3, 3))
                    gas_p(k, 1) = gas_p(k, 1) + (2._wp/3._wp)*dt*(gas_dpdt(k, 1)/4._wp + gas_dpdt(k, 2)/4._wp + gas_dpdt(k, 3))
                    gas_mv(k, 1) = gas_mv(k, 1) + (2._wp/3._wp)*dt*(gas_dmvdt(k, 1)/4._wp + gas_dmvdt(k, 2)/4._wp + gas_dmvdt(k, 3))
                    if (lag_params%coatedBub_model .and. (mrmtnt_shell(k, 2) == 0._wp)) then
                        if (intfc_rad(k, 1) < mrmtnt_Rrupt(k)) mrmtnt_shell(k, 2) = 1._wp ! No actual rupture happened during dt
                    end if
                    mrmtnt_shell(k, 1) = mrmtnt_shell(k, 2)
                    intfc_ac(k, 1) = (2._wp/3._wp)*(intfc_dveldt(k, 1)/4._wp + intfc_dveldt(k, 2)/4._wp + intfc_dveldt(k, 3))
                    if (polytropic) gas_p(k, 1) = pv + (gas_p(k, 2) - pv)*(bub_R0(k)/intfc_rad(k, 1))**(3._wp*gamma_m)
                end do

                call s_transfer_data_to_tmp
                call s_calculate_lag_bubble_stats()
                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[gas_p,gas_mv,intfc_rad,intfc_vel]')
                    call s_write_lag_particles(mytime, replace=.false.)
                end if
                call s_write_void_evol(mytime, replace=.false.)

            end if

        end if

    end subroutine s_update_lagrange_tdv_rk

    !> This subroutine returns the computational coordinate of the cell for the given position.
          !! @param pos Input coordinates
          !! @param cell Computational coordinate of the cell
          !! @param scoord Calculated particle coordinates
    pure subroutine s_locate_cell(pos, cell, scoord)

        real(wp), dimension(3), intent(in) :: pos
        real(wp), dimension(3), intent(out) :: scoord
        integer, dimension(3), intent(inout) :: cell

        integer :: i

        do while (pos(1) < x_cb(cell(1) - 1))
            cell(1) = cell(1) - 1
        end do

        do while (pos(1) > x_cb(cell(1)))
            cell(1) = cell(1) + 1
        end do

        do while (pos(2) < y_cb(cell(2) - 1))
            cell(2) = cell(2) - 1
        end do

        do while (pos(2) > y_cb(cell(2)))
            cell(2) = cell(2) + 1
        end do

        if (p > 0) then
            do while (pos(3) < z_cb(cell(3) - 1))
                cell(3) = cell(3) - 1
            end do
            do while (pos(3) > z_cb(cell(3)))
                cell(3) = cell(3) + 1
            end do
        end if

        ! The numbering of the cell of which left boundary is the domain boundary is 0.
        ! if comp.coord of the pos is s, the real coordinate of s is
        ! (the coordinate of the left boundary of the Floor(s)-th cell)
        ! + (s-(int(s))*(cell-width).
        ! In other words,  the coordinate of the center of the cell is x_cc(cell).

        !coordinates in computational space
        scoord(1) = cell(1) + (pos(1) - x_cb(cell(1) - 1))/dx(cell(1))
        scoord(2) = cell(2) + (pos(2) - y_cb(cell(2) - 1))/dy(cell(2))
        scoord(3) = 0._wp
        if (p > 0) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3) - 1))/dz(cell(3))
        cell(:) = int(scoord(:))
        do i = 1, num_dims
            if (scoord(i) < 0._wp) cell(i) = cell(i) - 1
        end do

    end subroutine s_locate_cell

    !> This subroutine transfer data into the temporal variables.
    impure subroutine s_transfer_data_to_tmp()

        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, nBubs
            if (.not. polytropic) gas_p(k, 2) = gas_p(k, 1)
            gas_mv(k, 2) = gas_mv(k, 1)
            intfc_rad(k, 2) = intfc_rad(k, 1)
            intfc_vel(k, 2) = intfc_vel(k, 1)
            intfc_ac(k, 2) = intfc_ac(k, 1)
            mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1)
            mtn_posPrev(k, 1:3, 2) = mtn_posPrev(k, 1:3, 1)
            mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1)
            mtn_s(k, 1:3, 2) = mtn_s(k, 1:3, 1)
            mrmtnt_shell(k, 2) = mrmtnt_shell(k, 1)
        end do

    end subroutine s_transfer_data_to_tmp

    !> The purpose of this procedure is to determine if the global coordinates of the bubbles
        !!      are present in the current MPI processor (including ghost cells).
        !! @param pos_part Spatial coordinates of the bubble
    ! pure function particle_in_domain(pos_part, restartFlag)

    !     logical :: particle_in_domain
    !     real(wp), dimension(3), intent(in) :: pos_part
    !     logical, intent(in) :: restartFlag

    !     real(wp) :: pos_part_radial

    !     ! 2D
    !     if (p == 0 .and. cyl_coord .neqv. .true.) then
    !         ! Defining a virtual z-axis that has the same dimensions as y-axis
    !         ! defined in the input file
    !         particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
    !                               (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
    !                               (pos_part(3) < lag_params%charwidth/2._wp) .and. (pos_part(3) >= -lag_params%charwidth/2._wp))
    !     else
    !         ! cyl_coord
    !         if (restartFlag) then
    !             pos_part_radial = pos_part(2)
    !         else
    !             pos_part_radial = sqrt(pos_part(2)**2._wp + pos_part(3)**2._wp)
    !         end if

    !         particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
    !                               (pos_part_radial < y_cb(n + buff_size)) .and. (pos_part_radial >= max(y_cb(-buff_size - 1), 0._wp)))
    !     end if

    !     ! 3D
    !     if (p > 0) then
    !         particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
    !                               (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
    !                               (pos_part(3) < z_cb(p + buff_size)) .and. (pos_part(3) >= z_cb(-buff_size - 1)))
    !     end if

    !     ! For symmetric boundary condition
    !     if (bc_x%beg == BC_REFLECTIVE) then
    !         particle_in_domain = (particle_in_domain .and. (pos_part(1) >= x_cb(-1)))
    !     end if
    !     if (bc_x%end == BC_REFLECTIVE) then
    !         particle_in_domain = (particle_in_domain .and. (pos_part(1) < x_cb(m)))
    !     end if
    !     if (bc_y%beg == BC_REFLECTIVE .and. (.not. cyl_coord)) then
    !         particle_in_domain = (particle_in_domain .and. (pos_part(2) >= y_cb(-1)))
    !     end if
    !     if (bc_y%end == BC_REFLECTIVE .and. (.not. cyl_coord)) then
    !         particle_in_domain = (particle_in_domain .and. (pos_part(2) < y_cb(n)))
    !     end if

    !     if (p > 0) then
    !         if (bc_z%beg == BC_REFLECTIVE) then
    !             particle_in_domain = (particle_in_domain .and. (pos_part(3) >= z_cb(-1)))
    !         end if
    !         if (bc_z%end == BC_REFLECTIVE) then
    !             particle_in_domain = (particle_in_domain .and. (pos_part(3) < z_cb(p)))
    !         end if
    !     end if

    ! end function particle_in_domain

    !> The purpose of this procedure is to determine if the lagrangian bubble is located in the
        !!       physical domain. The ghost cells are not part of the physical domain.
        !! @param pos_part Spatial coordinates of the bubble
    pure function particle_in_domain_physical(pos_part)

        logical :: particle_in_domain_physical
        real(wp), dimension(3), intent(in) :: pos_part

        particle_in_domain_physical = ((pos_part(1) < x_cb(m)) .and. (pos_part(1) >= x_cb(-1)) .and. &
                                       (pos_part(2) < y_cb(n)) .and. (pos_part(2) >= y_cb(-1)))

        if (p > 0) then
            particle_in_domain_physical = (particle_in_domain_physical .and. (pos_part(3) < z_cb(p)) .and. (pos_part(3) >= z_cb(-1)))
        end if

    end function particle_in_domain_physical

    !> The purpose of this procedure is to calculate the gradient of a scalar field along the x, y and z directions
        !!      following a second-order central difference considering uneven widths
        !! @param q Input scalar field
        !! @param dq Output gradient of q
        !! @param dir Gradient spatial direction
    pure subroutine s_gradient_dir(q, dq, dir)

        type(scalar_field), intent(inout) :: q
        type(scalar_field), intent(inout) :: dq
        integer, intent(in) :: dir

        integer :: i, j, k

        if (dir == 1) then
            ! Gradient in x dir.
            $:GPU_PARALLEL_LOOP(collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq%sf(i, j, k) = 0._wp
                        if (q%sf(i, j, k) /= q%sf(i + 1, j, k) .or. q%sf(i, j, k) /= q%sf(i - 1, j, k) .or. &
                            q%sf(i - 1, j, k) /= q%sf(i + 1, j, k)) then

                            dq%sf(i, j, k) = q%sf(i, j, k)*(dx(i + 1) - dx(i - 1)) &
                                             + q%sf(i + 1, j, k)*(dx(i) + dx(i - 1)) &
                                             - q%sf(i - 1, j, k)*(dx(i) + dx(i + 1))
                            dq%sf(i, j, k) = dq%sf(i, j, k)/ &
                                             ((dx(i) + dx(i - 1))*(dx(i) + dx(i + 1)))
                        end if
                    end do
                end do
            end do
        else
            if (dir == 2) then
                ! Gradient in y dir.
                $:GPU_PARALLEL_LOOP(collapse=3)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            dq%sf(i, j, k) = 0._wp
                            if (q%sf(i, j, k) /= q%sf(i, j + 1, k) .or. q%sf(i, j, k) /= q%sf(i, j - 1, k) .or. &
                                q%sf(i, j - 1, k) /= q%sf(i, j + 1, k)) then

                                dq%sf(i, j, k) = q%sf(i, j, k)*(dy(j + 1) - dy(j - 1)) &
                                                 + q%sf(i, j + 1, k)*(dy(j) + dy(j - 1)) &
                                                 - q%sf(i, j - 1, k)*(dy(j) + dy(j + 1))
                                dq%sf(i, j, k) = dq%sf(i, j, k)/ &
                                                 ((dy(j) + dy(j - 1))*(dy(j) + dy(j + 1)))
                            end if
                        end do
                    end do
                end do
            else
                ! Gradient in z dir.
                $:GPU_PARALLEL_LOOP(collapse=3)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            dq%sf(i, j, k) = 0._wp
                            if (q%sf(i, j, k) /= q%sf(i, j, k + 1) .or. q%sf(i, j, k) /= q%sf(i, j, k - 1) .or. &
                                q%sf(i, j, k - 1) /= q%sf(i, j, k + 1)) then

                                dq%sf(i, j, k) = q%sf(i, j, k)*(dz(k + 1) - dz(k - 1)) &
                                                 + q%sf(i, j, k + 1)*(dz(k) + dz(k - 1)) &
                                                 - q%sf(i, j, k - 1)*(dz(k) + dz(k + 1))
                                dq%sf(i, j, k) = dq%sf(i, j, k)/ &
                                                 ((dz(k) + dz(k - 1))*(dz(k) + dz(k + 1)))
                            end if
                        end do
                    end do
                end do
            end if
        end if

        ! call s_mpi_barrier()

    end subroutine s_gradient_dir

    !> Subroutine that writes on each time step the changes of the lagrangian bubbles.
        !!  @param q_time Current time
    impure subroutine s_write_lag_particles(qtime, replace)

        real(wp), intent(in) :: qtime
        logical, intent(in) :: replace
        integer :: k

        logical :: file_exist
        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'lag_bubble_evol_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (.not. file_exist .or. replace) then
            open (11, FILE=trim(file_loc), FORM='formatted', position='rewind')
            write (11, *) 'currentTime, particleID, x, y, z, ', &
                'coreVaporMass, coreVaporConcentration, radius, interfaceVelocity, ', &
                'corePressure'
        else
            open (11, FILE=trim(file_loc), FORM='formatted', position='append')
        end if

        if (lag_params%write_only_bub_id == dflt_int) then
            ! Cycle through list
            do k = 1, nBubs

                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    write (11, '(6X,f12.6,I12.6,13e24.8)') &
                        qtime, &
                        lag_id(k, 1), &
                        mtn_pos(k, 1, 1), &
                        mtn_pos(k, 2, 1), &
                        mtn_pos(k, 3, 1), &
                        gas_mv(k, 1), &
                        gas_mv(k, 1)/(gas_mv(k, 1) + gas_mg(k)), &
                        intfc_rad(k, 1), &
                        intfc_vel(k, 1), &
                        gas_p(k, 1), &
                        bub_interact(k), &
                        bub_qvis(k), &
                        bub_qth(k), &
                        mrmtnt_shell(k, 1), &
                        mrmtnt_Rrupt(k)
                end if
            end do

        elseif (nBubs > 0) then
            ! One specified bubble id only
            k = lag_params%write_only_bub_id

            if (k == lag_id(k, 1)) then
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    write (11, '(6X,f12.6,I12.6,12e24.8)') &
                        qtime, &
                        lag_id(k, 1), &
                        mtn_pos(k, 1, 1), &
                        mtn_pos(k, 2, 1), &
                        mtn_pos(k, 3, 1), &
                        gas_mv(k, 1), &
                        gas_mv(k, 1)/(gas_mv(k, 1) + gas_mg(k)), &
                        intfc_rad(k, 1), &
                        intfc_vel(k, 1), &
                        gas_p(k, 1), &
                        bub_interact(k), &
                        bub_qvis(k), &
                        bub_qth(k), &
                        mrmtnt_shell(k, 1)
                end if

            end if
        end if

        close (11)

    end subroutine s_write_lag_particles

    !>  Subroutine that writes some useful statistics related to the volume fraction
            !!       of the particles (void fraction) in the computatioational domain
            !!       on each time step.
            !!  @param q_time Current time
    impure subroutine s_write_void_evol(qtime, replace)

        real(wp), intent(in) :: qtime
        logical, intent(in) :: replace
        real(wp) :: volcell, voltot
        real(wp) :: lag_void_max, lag_void_avg, lag_vol
        real(wp) :: void_max_glb, void_avg_glb, vol_glb
        real(wp) :: aux_glb, nBubs_all

        integer :: i, j, k

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        if (proc_rank == 0) then
            write (file_loc, '(A)') 'voidfraction.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (.not. file_exist .or. replace) then
                open (12, FILE=trim(file_loc), FORM='formatted', position='rewind')
                !write (12, *) 'currentTime, averageVoidFraction, ', &
                !    'maximumVoidFraction, totalParticlesVolume', 'maxRadius', 'minRadius'
                !write (12, *) 'The averageVoidFraction value does ', &
                !    'not reflect the real void fraction in the cloud since the ', &
                !    'cells which do not have bubbles are not accounted'
            else
                open (12, FILE=trim(file_loc), FORM='formatted', position='append')
            end if
        end if

        lag_void_max = 0._wp
        lag_void_avg = 0._wp
        lag_vol = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, reduction='[[lag_vol, lag_void_avg], &
            & [lag_void_max]]', reductionOp='[+, MAX]', &
            & copy='[lag_vol, lag_void_avg, lag_void_max]')
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    lag_void_max = max(lag_void_max, 1._wp - q_beta%vf(1)%sf(i, j, k))
                    call s_get_char_vol(i, j, k, volcell)
                    if ((1._wp - q_beta%vf(1)%sf(i, j, k)) > 5.0d-11) then
                        lag_void_avg = lag_void_avg + (1._wp - q_beta%vf(1)%sf(i, j, k))*volcell
                        lag_vol = lag_vol + volcell
                    end if
                end do
            end do
        end do
        nBubs_all = real(nBubs, wp)

$:GPU_UPDATE(host='[Rmax_glb, Rmin_glb, Rmean_glb]')
        
#ifdef MFC_MPI
        if (num_procs > 1) then
            call s_mpi_allreduce_max(lag_void_max, void_max_glb)
            lag_void_max = void_max_glb
            call s_mpi_allreduce_sum(lag_vol, vol_glb)
            lag_vol = vol_glb
            call s_mpi_allreduce_sum(lag_void_avg, void_avg_glb)
            lag_void_avg = void_avg_glb
            call s_mpi_allreduce_max(Rmax_glb, aux_glb)
            Rmax_glb = aux_glb
            call s_mpi_allreduce_min(Rmin_glb, aux_glb)
            Rmin_glb = aux_glb
            call s_mpi_allreduce_sum(Rmean_glb, aux_glb)
            Rmean_glb = aux_glb
            call s_mpi_allreduce_sum(nBubs_all, aux_glb)
            nBubs_all = aux_glb
        end if
#endif
        voltot = lag_void_avg

        ! This voidavg value does not reflect the real void fraction in the cloud
        ! since the cell which does not have bubbles are not accounted
        if (lag_vol > 0._wp) lag_void_avg = lag_void_avg/lag_vol

        if (proc_rank == 0) then

            if (hifu_params%moments) then
            write (12, '(6X,8e24.8)') &
                qtime, &
                lag_void_avg, &
                lag_void_max, &
                voltot, &
                Rmean_glb/nBubs_all, &
                nBubs_all, &
                Rmax_glb, &
                Rmin_glb
            else
                write (12, '(6X,4e24.8)') &
                qtime, &
                lag_void_avg, &
                lag_void_max, &
                voltot
            end if
            close (12)
        end if

    end subroutine s_write_void_evol

    !>  Subroutine that writes the restarting files for the particles in the lagrangian solver.
        !!  @param t_step Current time step
    impure subroutine s_write_restart_lag_bubbles(t_step)

        ! Generic string used to store the address of a particular file
        integer, intent(in) :: t_step

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist
        integer :: bub_id, tot_part, tot_part_wrtn, npart_wrtn
        integer :: i, k

#ifdef MFC_MPI
        ! For Parallel I/O
        integer :: ifile, ierr
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer :: view
        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer, dimension(num_procs) :: part_order, part_ord_mpi
        integer :: varsExtra = 7

        bub_id = 0._wp
        if (nBubs /= 0) then
            do k = 1, nBubs
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    bub_id = bub_id + 1
                end if
            end do
        end if

        if (.not. parallel_io) return

        ! Total number of particles
        call MPI_ALLREDUCE(bub_id, tot_part, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        ! Total number of particles written so far
        call MPI_ALLREDUCE(npart_wrtn, tot_part_wrtn, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        lsizes(1) = max(1, bub_id)
        lsizes(2) = 21 + varsExtra

        ! if the particle number is zero, put 1 since MPI cannot deal with writing
        ! zero particle
        part_order(:) = 1
        part_order(proc_rank + 1) = max(1, bub_id)

        call MPI_ALLREDUCE(part_order, part_ord_mpi, num_procs, MPI_integer, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

        gsizes(1) = sum(part_ord_mpi(1:num_procs))
        gsizes(2) = 21 + varsExtra

        start_idx_part(1) = sum(part_ord_mpi(1:proc_rank + 1)) - part_ord_mpi(proc_rank + 1)
        start_idx_part(2) = 0

        write (file_loc, '(A,I0,A)') 'lag_bubbles_mpi_io_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if

        ! Writing down the total number of particles
        if (proc_rank == 0) then
            open (9, FILE=trim(file_loc), FORM='unformatted', STATUS='unknown')
            write (9) gsizes(1), mytime, dt
            close (9)
        end if

        call MPI_type_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        allocate (MPI_IO_DATA_lag_bubbles(1:max(1, bub_id), 1:(21 + varsExtra)))

        ! Open the file to write all flow variables
        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                           mpi_info_int, ifile, ierr)

        disp = 0._wp

        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, &
                               'native', mpi_info_null, ierr)

        ! Cycle through list
        i = 1

        if (bub_id == 0) then
            MPI_IO_DATA_lag_bubbles(1, 1:(21 + varsExtra)) = 0._wp
        else

            do k = 1, nBubs

                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then

                    MPI_IO_DATA_lag_bubbles(i, 1) = real(lag_id(k, 1))
                    MPI_IO_DATA_lag_bubbles(i, 2:4) = mtn_pos(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 5:7) = mtn_posPrev(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 8:10) = mtn_vel(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 11) = intfc_rad(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 12) = intfc_vel(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 13) = bub_R0(k)
                    MPI_IO_DATA_lag_bubbles(i, 14) = Rmax_stats(k)
                    MPI_IO_DATA_lag_bubbles(i, 15) = Rmin_stats(k)
                    MPI_IO_DATA_lag_bubbles(i, 16) = bub_dphidt(k)
                    if (.not. polytropic) then
                        MPI_IO_DATA_lag_bubbles(i, 17) = gas_p(k, 1)
                    else
                        MPI_IO_DATA_lag_bubbles(i, 17) = gas_p(k, 2)
                    end if
                    MPI_IO_DATA_lag_bubbles(i, 18) = gas_mv(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 19) = gas_mg(k)
                    MPI_IO_DATA_lag_bubbles(i, 20) = gas_betaT(k)
                    MPI_IO_DATA_lag_bubbles(i, 21) = gas_betaC(k)
                    ! Marmotant
                    MPI_IO_DATA_lag_bubbles(i, 22) = mrmtnt_shell(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 23) = mrmtnt_Rbuck(k)
                    MPI_IO_DATA_lag_bubbles(i, 24) = mrmtnt_Rrupt(k)
                    ! hifu
                    MPI_IO_DATA_lag_bubbles(i, 25) = bub_qvis(k)
                    MPI_IO_DATA_lag_bubbles(i, 26) = bub_qth(k)
                    MPI_IO_DATA_lag_bubbles(i, 27) = intfc_ac(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 28) = bub_hifu_rad(k)

                    i = i + 1

                end if

            end do

        end if

        call MPI_FILE_write_ALL(ifile, MPI_IO_DATA_lag_bubbles, (21 + varsExtra)*max(1, bub_id), &
                                mpi_p, status, ierr)

        call MPI_FILE_CLOSE(ifile, ierr)

        deallocate (MPI_IO_DATA_lag_bubbles)

#endif

    end subroutine s_write_restart_lag_bubbles

    !>  This procedure calculates the maximum and minimum radius of each bubble.
    subroutine s_calculate_lag_bubble_stats()

        integer :: k

        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)
        Rmean_glb = 0._wp
        $:GPU_UPDATE(device='[Rmax_glb, Rmin_glb, Rmean_glb]')

        $:GPU_PARALLEL_LOOP(reduction='[[Rmax_glb], [Rmin_glb], [Rmean_glb]]', &
            & reductionOp='[MAX, MIN, +]', copy='[Rmax_glb,Rmin_glb,Rmean_glb]')
        do k = 1, nBubs
            Rmax_glb = max(Rmax_glb, intfc_rad(k, 1))
            Rmin_glb = min(Rmin_glb, intfc_rad(k, 1))
            Rmean_glb = Rmean_glb + intfc_rad(k, 1)
            Rmax_stats(k) = max(Rmax_stats(k), intfc_rad(k, 1)/bub_R0(k))
            Rmin_stats(k) = min(Rmin_stats(k), intfc_rad(k, 1)/bub_R0(k))
        end do

    end subroutine s_calculate_lag_bubble_stats

    !>  Subroutine that writes the maximum and minimum radius of each bubble.
    impure subroutine s_write_lag_bubble_stats()

        integer :: k
        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'stats_lag_bubbles_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        $:GPU_UPDATE(host='[Rmax_glb,Rmin_glb,Rmean_glb]')

        open (13, FILE=trim(file_loc), FORM='formatted', position='rewind')
        write (13, *) 'proc_rank, Rmax_glb, Rmin_glb, Rmean_glb'
        write (13, '(6X,I24.8,3e24.8)') &
            proc_rank, &
            Rmax_glb, &
            Rmin_glb, &
            Rmean_glb/nBubs
        write (13, *) ' '
        write (13, *) 'particleID, x, y, z, Rmax, Rmin'

        do k = 1, nBubs
            write (13, '(6X,I24.8,7e24.8)') &
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                Rmax_stats(k), &
                Rmin_stats(k), &
                gas_betaT(k), &
                gas_betaC(k)
        end do

        close (13)

    end subroutine s_write_lag_bubble_stats

    !> The purpose of this subroutine is to remove one specific particle if dt is too small.
          !! @param bub_id Particle id
    impure subroutine s_remove_lag_bubble(bub_id)

        integer, intent(in) :: bub_id

        integer :: i

        $:GPU_LOOP(parallelism='[seq]')
        do i = bub_id, nBubs - 1
            if (i == bub_id) print *, 'In loop remove bub:', i
            lag_id(i, 1) = lag_id(i + 1, 1)
            bub_R0(i) = bub_R0(i + 1)
            Rmax_stats(i) = Rmax_stats(i + 1)
            Rmin_stats(i) = Rmin_stats(i + 1)
            gas_mg(i) = gas_mg(i + 1)
            gas_betaT(i) = gas_betaT(i + 1)
            gas_betaC(i) = gas_betaC(i + 1)
            bub_dphidt(i) = bub_dphidt(i + 1)
            gas_p(i, 1:2) = gas_p(i + 1, 1:2)
            gas_mv(i, 1:2) = gas_mv(i + 1, 1:2)
            intfc_rad(i, 1:2) = intfc_rad(i + 1, 1:2)
            intfc_vel(i, 1:2) = intfc_vel(i + 1, 1:2)
            intfc_ac(i, 1:2) = intfc_ac(i + 1, 1:2)
            mtn_pos(i, 1:3, 1:2) = mtn_pos(i + 1, 1:3, 1:2)
            mtn_posPrev(i, 1:3, 1:2) = mtn_posPrev(i + 1, 1:3, 1:2)
            mtn_vel(i, 1:3, 1:2) = mtn_vel(i + 1, 1:3, 1:2)
            mtn_s(i, 1:3, 1:2) = mtn_s(i + 1, 1:3, 1:2)
            intfc_draddt(i, 1:lag_num_ts) = intfc_draddt(i + 1, 1:lag_num_ts)
            intfc_dveldt(i, 1:lag_num_ts) = intfc_dveldt(i + 1, 1:lag_num_ts)
            gas_dpdt(i, 1:lag_num_ts) = gas_dpdt(i + 1, 1:lag_num_ts)
            gas_dmvdt(i, 1:lag_num_ts) = gas_dmvdt(i + 1, 1:lag_num_ts)
            ! mtn_dposdt(i, 1:3, 1:lag_num_ts) = mtn_dposdt(i + 1, 1:3, 1:lag_num_ts)
            ! mtn_dveldt(i, 1:3, 1:lag_num_ts) = mtn_dveldt(i + 1, 1:3, 1:lag_num_ts)
            mrmtnt_shell(i, 1:2) = mrmtnt_shell(i + 1, 1:2)
            mrmtnt_Rbuck(i) = mrmtnt_Rbuck(i + 1)
            mrmtnt_Rrupt(i) = mrmtnt_Rrupt(i + 1)
            bub_qvis(i) = bub_qvis(i + 1)
            bub_qth(i) = bub_qth(i + 1)
            bub_hifu_rad(i) = bub_hifu_rad(i + 1)
        end do

        nBubs = nBubs - 1
        dt = 5._wp*dt
        $:GPU_UPDATE(device='[nBubs, dt]')

        print *, 'Bubble removed, nBubs now: in processor: ', nBubs, proc_rank

    end subroutine s_remove_lag_bubble

    subroutine s_free_memory_stg3()

        @:DEALLOCATE(Rmax_stats)
        @:DEALLOCATE(Rmin_stats)
        @:DEALLOCATE(gas_mg)
        @:DEALLOCATE(gas_betaT)
        @:DEALLOCATE(gas_betaC)
        @:DEALLOCATE(bub_dphidt)
        @:DEALLOCATE(gas_p)
        @:DEALLOCATE(gas_mv)
        @:DEALLOCATE(intfc_ac)
        @:DEALLOCATE(mtn_vel)
        @:DEALLOCATE(intfc_draddt)
        @:DEALLOCATE(intfc_dveldt)
        @:DEALLOCATE(gas_dpdt)
        @:DEALLOCATE(gas_dmvdt)
        ! @:DEALLOCATE(mtn_dposdt)
        ! @:DEALLOCATE(mtn_dveldt)
        ! Marmotant model
        @:DEALLOCATE(mrmtnt_shell)
        @:DEALLOCATE(mrmtnt_Rbuck)
        @:DEALLOCATE(mrmtnt_Rrupt)
        ! bubble interaction
        @:DEALLOCATE(bub_interact)
        if (lag_params%pressure_corrector .and. any(lag_params%interaction_model == (/2, 3/))) then
            @:DEALLOCATE(bub_int_ids)
        end if
        !@:DEALLOCATE(bub_lambda_c)
        if (hifu_params%moments) then
            @:DEALLOCATE(moments_bubs)
        end if
        if (hifu_params%power_balance) then
            @:DEALLOCATE(acPw_bubs)
        end if

    end subroutine s_free_memory_stg3

    !> The purpose of this subroutine is to deallocate variables
    impure subroutine s_finalize_lagrangian_solver()

        integer :: i

        do i = 1, q_beta_idx
            @:DEALLOCATE(q_beta%vf(i)%sf)
        end do
        @:DEALLOCATE(q_beta%vf)

        !Deallocating space
        @:DEALLOCATE(lag_id)
        @:DEALLOCATE(bub_R0)
        @:DEALLOCATE(intfc_rad)
        @:DEALLOCATE(intfc_vel)
        @:DEALLOCATE(mtn_pos)
        @:DEALLOCATE(mtn_posPrev)
        @:DEALLOCATE(mtn_s)
        ! hifu
        @:DEALLOCATE(bub_qvis)
        @:DEALLOCATE(bub_qth)
        @:DEALLOCATE(bub_hifu_rad)

        if (.not. hifu_params%heatSolver) then
            @:DEALLOCATE(Rmax_stats)
            @:DEALLOCATE(Rmin_stats)
            @:DEALLOCATE(gas_mg)
            @:DEALLOCATE(gas_betaT)
            @:DEALLOCATE(gas_betaC)
            @:DEALLOCATE(bub_dphidt)
            @:DEALLOCATE(gas_p)
            @:DEALLOCATE(gas_mv)
            @:DEALLOCATE(intfc_ac)
            @:DEALLOCATE(mtn_vel)
            @:DEALLOCATE(intfc_draddt)
            @:DEALLOCATE(intfc_dveldt)
            @:DEALLOCATE(gas_dpdt)
            @:DEALLOCATE(gas_dmvdt)
            ! @:DEALLOCATE(mtn_dposdt)
            ! @:DEALLOCATE(mtn_dveldt)
            ! Marmotant model
            @:DEALLOCATE(mrmtnt_shell)
            @:DEALLOCATE(mrmtnt_Rbuck)
            @:DEALLOCATE(mrmtnt_Rrupt)
            ! bubble interaction
            @:DEALLOCATE(bub_interact)
            if (lag_params%pressure_corrector .and. any(lag_params%interaction_model == (/2, 3/))) then
                @:DEALLOCATE(bub_int_ids)
            end if
            !@:DEALLOCATE(bub_lambda_c)
            if (hifu_params%moments) then
                @:DEALLOCATE(moments_bubs)
            end if
            if (hifu_params%power_balance) then
                @:DEALLOCATE(acPw_bubs)
            end if
        end if

    end subroutine s_finalize_lagrangian_solver

end module m_bubbles_EL
