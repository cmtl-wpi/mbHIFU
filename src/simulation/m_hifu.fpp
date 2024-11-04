!>
!! @file m_hifu.f90
!! @brief Contains module m_viscous

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to study HIFU
module m_hifu

    ! Dependencies ============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    !use m_viscous

    ! ==========================================================================
    implicit none
    private; public ::  s_update_hifu_vars_stg2, &
                        s_update_Pmax, &
                        s_update_intensity_HIFU, &
                        s_initialize_HIFU, &
                        s_rhs_heatEqn, &
                        s_populate_HIFU_variables_buffers, &
                        s_open_run_time_information_samplingHIFU, &
                        s_close_run_time_information_samplingHIFU

contains

    ! ==========================================================================
    ! Important variables and where they are stored:
    ! q_cons_hifu(1)%sf(j,k,l): Temperature distribution                (T_hifu_idx)
    ! q_cons_hifu(2)%sf(j,k,l): RHS value from heat transfer eqn        (T_hifu_idx+1)
    ! q_cons_hifu(3)%sf(j,k,l): Maximum pressure at each cell           (P_hifu_idx)
    ! q_cons_hifu(4)%sf(j,k,l): Total sampling time                     (tt_hifu_idx)
    ! q_cons_hifu(5)%sf(j,k,l): Acoustic damping (Sum over time)        (qus_hifu_idx / qus_prms_hifu_idx)
    ! q_cons_hifu(6)%sf(j,k,l): Viscous damping  (Sum over time)        (qvis_hifu_idx)
    ! q_cons_hifu(7)%sf(j,k,l): Any extra variable / Analit. Sol.       (dmb_hifu_idx)
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ==========================================================================

    ! ==========================================================================
    ! Initializes the HIFU solver using the MFC arquitecture
    subroutine s_initialize_HIFU(q_cons_hifu)

        type(scalar_field), dimension(sys_size_hifu), intent(inout) :: q_cons_hifu
        type(int_bounds_info) :: ix_t, iy_t, iz_t !<
        integer :: i, j !< Generic loop iterators

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

        !Zeroing all the hifu variables
        do j = 1, sys_size_hifu
            q_cons_hifu(j)%sf(ix_t%beg:ix_t%end, &
                              iy_t%beg:iy_t%end, &
                        iz_t%beg:iz_t%end) = 0.0d0
        end do

        !Initial Temperature
        q_cons_hifu(T_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
                                   iz_t%beg:iz_t%end) = hifu_Tref !User input

        !Initialize Pmax
        q_cons_hifu(P_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
                                   iz_t%beg:iz_t%end) =  min(dflt_real,-dflt_real)

        !Initialize Pmin
        q_cons_hifu(P_hifu_idx+1)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
                                   iz_t%beg:iz_t%end) =  max(dflt_real,-dflt_real)

        call s_open_run_time_information_samplingHIFU

    end subroutine s_initialize_HIFU ! =========================================

    subroutine s_update_HIFU_vars_stg2(q_cons_sf, q_prim_sf, q_cons_hifu, t_step, hdid)

        type(scalar_field), dimension(sys_size_hifu), intent(in) :: q_cons_sf
        type(scalar_field), dimension(sys_size_hifu), intent(in) :: q_prim_sf
        type(scalar_field), dimension(sys_size_hifu), intent(inout) :: q_cons_hifu
        integer, intent(in) :: t_step
        real(kind(0.d0)), intent(in) :: hdid

        !> update Pmax
        call s_update_Pmax(q_cons_sf, q_cons_hifu, t_step)

        !> Update heat deposition source terms
        call s_update_intensity_HIFU(q_cons_sf, q_prim_sf, q_cons_hifu, t_step, hdid)

    end subroutine s_update_HIFU_vars_stg2 
    

    ! Update the maximum and minimum pressure through time along the 
    ! axisymmetric and radial axes
    subroutine s_update_Pmax(q_cons_vf,q_cons_hifu,t_step)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        type(scalar_field), dimension(sys_size_hifu), intent(INOUT) :: q_cons_hifu
        integer, intent(IN) :: t_step
        
        logical :: axialCondition, radialCondition, condition
        integer :: i, j, k, l, s, q  !< generic loop variables
        integer :: unitFile

	    real(kind(0d0)) :: rho
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: pres
	    real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: qv
        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)) :: G
        real(kind(0d0)) :: rhoYks(1:num_species)

        ! Zeroing out flow variables for all processors      
        rho = 0d0
        do s = 1, num_dims
           vel(s) = 0d0
        end do
        pres = 0d0
        gamma = 0d0
        pi_inf = 0d0

        unitFile = 100+proc_rank

        if (cyl_coord .and. p==0) then
            l = 0
            ! if (mod(t_step,t_step_save)==0) then
            !     ! write (file_path, '(A,I0,A)') '/D/Pmax_', proc_rank, '.dat'
            !     ! file_path = trim(case_dir)//trim(file_path)        
            !     ! open (unitFile, FILE=trim(file_path), FORM='formatted', STATUS='unknown')
            !     ! write (unitFile, *) 'timeStep, x_cc, y_cc, Pmax, Pmin, j, k, l'
            ! end if

            do  j = 0, m
                do k = 0, n
                    call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                            rho, gamma, pi_inf, qv, &
                                                            Re, G, fluid_pp(:)%G)
                    do s = 1, num_dims
                        vel(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho
                    end do

                    call s_compute_pressure(q_cons_vf(E_idx)%sf(j, k, l), &
                                0d0, 0.5d0*rho*dot_product(vel, vel), pi_inf, gamma, rho, qv, rhoYks, pres)

                    q_cons_hifu(P_hifu_idx)%sf(j,k,l)   = max(q_cons_hifu(P_hifu_idx)%sf(j,k,l),pres)
                    q_cons_hifu(P_hifu_idx+1)%sf(j,k,l) = min(q_cons_hifu(P_hifu_idx+1)%sf(j,k,l),pres)
                    
                    ! Specify enough conditions for axial and radial probe lines
                    if (mod(t_step,t_step_save)==0) then
                        axialCondition= (dy(k)>y_cc(k) .and. y_cc(k)>0)
                        radialCondition= (x_cb(j-1)<focLength_bc .and. focLength_bc<x_cb(j))
                        condition= (axialCondition .or. radialCondition)
                        if (condition) then
                            write (100, '(6x,I24,4E24.8)') &
                                t_step, &
                                x_cc(j), &
                                y_cc(k), &
                                q_cons_hifu(P_hifu_idx)%sf(j,k,l), &
                                q_cons_hifu(P_hifu_idx+1)%sf(j,k,l)
                        end if
                        
                    end if
                end do
            end do

            ! if (mod(t_step,t_step_save)==0) then
            !     close(unitFile)
            ! end if

        end if
   
    end subroutine s_update_Pmax ! =============================================

    ! ==========================================================================
    !> Heat deposition in HIFU.
    !> Obtain the generated heat source "q_us_ac". Heating from the primary ultrasound source
    !> and, if mb exist, from acoustic emission from bubble oscillations.
    subroutine s_update_intensity_HIFU(q_cons_vf, q_prim_vf, q_cons_hifu, t_step, hdid)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size_hifu), intent(inout) :: q_cons_hifu
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: t_step
        real(kind(0.d0)), intent(in) :: hdid

        real(kind(0d0)) :: nondim_time !< Non-dimensional time	
        real(kind(0d0)) :: rho
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: pres
        real(kind(0d0)) :: pres_old
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)) :: G
        real(kind(0d0)) :: qv
        real(kind(0d0)) :: c
        real(kind(0d0)), dimension(num_fluids) :: alpha
        real(kind(0d0)) :: rhoYks(1:num_species)

        real(kind(0d0)) :: shearVisc, bulkVisc
        real(kind(0d0)) :: absCoef, spdsound 
        real(kind(0d0)) :: angFreq
        real(kind(0d0)) :: varA, varB
        real(kind(0d0)) :: duxdx, duxdr, durdx, durdr
        real(kind(0d0)) :: ep11, ep22, ep33, ep13
        real(kind(0d0)) :: intensity_ac, sumIntensity_ac, tmp, focalIntensity_ac, intensity_ac_prms
        real(kind(0d0)) :: focalIntensity_th, sumIntensity_th
        real(kind(0d0)) :: intensity_vis, sumIntensity_vis, focalIntensity_vis, focalIntensity_ac_prms
        real(kind(0d0)) :: focal_u, focal_v

        real(kind(0d0)), dimension(num_fluids) :: myalpha_rho, myalpha
        real(kind(0d0)) :: n_tait, B_tait, myRho, lamda

        integer :: ii, i, j, k, l, s, q  !< generic loop variables

        logical :: file_exist, printFlag1, printFlag2
        logical :: axialCondition, radialCondition, condition
        real(kind(0d0)) :: val, dist
        
        ! Zeroing out flow variables for all processors      
        rho = 0d0
        do s = 1, num_dims
           vel(s) = 0d0
        end do
        pres = 0d0
        pres_old = 0d0
        qv = 0d0
        c = 0d0
        gamma = 0d0
        pi_inf = 0d0
        focalIntensity_ac = 0d0
        focalIntensity_ac_prms = 0d0
        focalIntensity_vis = 0d0
        focalIntensity_th = 0d0
        sumIntensity_ac = 0d0
        sumIntensity_vis = 0d0
        sumIntensity_th = 0d0
        focal_u = 0.0d0
        focal_v = 0.0d0

        ! if (cfl_dt) then
        !     save_count_start = n_start
        ! else
        !     save_count_start = t_step_start
        ! end if


        if (cyl_coord .and. p==0) then  !Axysimetric		
            l = 0

            ! if (proc_rank==0 .and. t_step==save_count_start ) then
            !     write (file_path, '(A,I0,A)') '/D/sumIntensity-HIFU.dat'
            !     file_path = trim(case_dir)//trim(file_path)
            !     open (99, FILE=trim(file_path), FORM='formatted', POSITION='append', STATUS='unknown')
            !     write (99, *) 'timeStep, totalSamplingTime, acousticFocalIntensity, acousticFocalIntensityPRMS, viscousFocalIntensity, ', &
            !                   'thermalFocalIntensity, sumAcousticIntensity, sumViscousIntensity, sumThermalIntensity, focalxVel, focalyVel'
            ! end if

            do  j = 0, m
                do k = 0, n

                    !Get viscosities (and absorption coeff.) which are user inputs
                    shearVisc = 0.0d0
                    bulkVisc  = 0.0d0
                    absCoef   = 0.0d0
                    do i = 1, num_fluids
                        shearVisc = shearVisc + q_prim_vf(E_idx+i)%sf(j,k,l)*fluid_pp(i)%Re(1)
                        bulkVisc  = bulkVisc  + q_prim_vf(E_idx+i)%sf(j,k,l)*fluid_pp(i)%Re(2)
                        absCoef   = absCoef   + q_prim_vf(E_idx+i)%sf(j,k,l)*fluid_pp(i)%absCoef
                        alpha(i)  = q_prim_vf(E_idx+i)%sf(j,k,l)
                    end do
                    shearVisc = 1/shearVisc
                    bulkVisc  = 1/bulkVisc

                    if (absCoef<=0.0) call s_mpi_abort('Check absCoef values!')

                    !>> Get the strain rate tensor (using central finite difference)
                    varA = 0.0d0
                    varB = 0.0d0

                    ! Only for axysimmetric assumption
                    duxdx = (q_prim_vf(mom_idx%beg)%sf(j+1, k, 0) - q_prim_vf(mom_idx%beg)%sf(j-1, k, 0)) / (x_cc(j+1) - x_cc(j-1))
                    duxdr = (q_prim_vf(mom_idx%beg)%sf(j, k+1, 0) - q_prim_vf(mom_idx%beg)%sf(j, k-1, 0)) / (y_cc(k+1) - y_cc(k-1))

                    durdx = (q_prim_vf(mom_idx%beg+1)%sf(j+1, k, 0) - q_prim_vf(mom_idx%beg+1)%sf(j-1, k, 0)) / (x_cc(j+1) - x_cc(j-1))
                    durdr = (q_prim_vf(mom_idx%beg+1)%sf(j, k+1, 0) - q_prim_vf(mom_idx%beg+1)%sf(j, k-1, 0)) / (y_cc(k+1) - y_cc(k-1))

                    !>> Get pressure, density and speed of sound
                    call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                            rho, gamma, pi_inf, qv, &
                                                            Re, G, fluid_pp(:)%G)
                    do s = 1, num_dims
                        vel(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho
                    end do

                    call s_compute_pressure(q_cons_vf(E_idx)%sf(j, k, l), &
                                0d0, 0.5d0*rho*dot_product(vel, vel), pi_inf, gamma, rho, qv, rhoYks, pres)

                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                      ((gamma + 1d0)*pres + pi_inf)/rho, alpha, 0d0, c)

                    !>> Compute intensity form acoustic damping
                    
                    intensity_ac_prms = 0d0
                    intensity_ac_prms = intensity_ac_prms + absCoef*(q_cons_hifu(P_hifu_idx)%sf(j,k,l)-hifu_atmPres)**2/(rho*c)

                    intensity_ac = 0d0
                    ep11 = durdr
                    ep22 = vel(2) / y_cc(k)
                    ep33 = duxdx
                    ep13 = 0.5d0*(durdx + duxdr)
                    varA = ep11**2.0 + ep22**2.0 + ep33**2.0
                    varB = (8.0d0/3.0d0)*varA - (4.0d0/3.0d0)*(ep11*ep22 + ep11*ep33 + ep22*ep33) + 6.0d0*(ep13**2.0)
                    intensity_ac = intensity_ac + bulkVisc*varA + 2.0d0*shearVisc*varB !intensity is "q_us_ac"

                    q_cons_hifu(tt_hifu_idx)%sf(j,k,l) = q_cons_hifu(tt_hifu_idx)%sf(j,k,l) + hdid ! Update total sampling time
                    q_cons_hifu(qus_hifu_idx)%sf(j,k,l) = q_cons_hifu(qus_hifu_idx)%sf(j,k,l) + intensity_ac * hdid ! Sampling acoustic intensity 
                    q_cons_hifu(qus_prms_hifu_idx)%sf(j,k,l) = q_cons_hifu(qus_prms_hifu_idx)%sf(j,k,l) + intensity_ac_prms * hdid ! Sampling acoustic intensity (prms)

                    !Update average velocities for streaming
                    q_cons_hifu(u_hifu_idx)%sf(j,k,l) =  q_cons_hifu(u_hifu_idx)%sf(j,k,l) + vel(1) * hdid ! Sampling x-vel
                    q_cons_hifu(v_hifu_idx)%sf(j,k,l) =  q_cons_hifu(v_hifu_idx)%sf(j,k,l) + vel(2) * hdid ! Sampling y-vel

                    !Get focal intensity and velocities
                    axialCondition= (dy(k)>y_cc(k) .and. y_cc(k)>0.0)
                    !radialCondition= (x_cb(j-1)<=focalPoint_x .and. focalPoint_x<=x_cb(j))
                    radialCondition= (x_cb(j-1)<focLength_bc .and. focLength_bc<x_cb(j))
                    condition= (axialCondition .and. radialCondition)
                    if (condition) then
                            focalIntensity_ac = q_cons_hifu(qus_hifu_idx)%sf(j,k,l)
                            focalIntensity_ac_prms = q_cons_hifu(qus_prms_hifu_idx)%sf(j,k,l)
                            focalIntensity_vis = q_cons_hifu(qvis_hifu_idx)%sf(j,k,l)
                            focalIntensity_th = q_cons_hifu(qth_hifu_idx)%sf(j,k,l)
                            focal_u = q_cons_hifu(u_hifu_idx)%sf(j,k,l)
                            focal_v = q_cons_hifu(v_hifu_idx)%sf(j,k,l)
                    end if

                    !Intensity summation through the domain, avoid acoustic source influence (0.8*Focal length)
                    sumIntensity_ac = sumIntensity_ac + q_cons_hifu(qus_hifu_idx)%sf(j,k,l)
                    sumIntensity_vis = sumIntensity_vis + q_cons_hifu(qvis_hifu_idx)%sf(j,k,l)
                    sumIntensity_th = sumIntensity_th + q_cons_hifu(qth_hifu_idx)%sf(j,k,l)

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

            if (proc_rank==0) write (99, '(6x,I24.8,f24.8,9e24.8)') &
                t_step, q_cons_hifu(tt_hifu_idx)%sf(0,0,0), focalIntensity_ac, &
                focalIntensity_ac_prms, focalIntensity_vis, focalIntensity_th, &
                sumIntensity_ac, sumIntensity_vis, sumIntensity_th, focal_u, focal_v

            !if (t_step == t_step_stop .and. proc_rank==0) close (99)
  
        end if

    end subroutine s_update_intensity_HIFU ! ===================================

    ! ==========================================================================
    !Calculate the rhs value from heat transfer eqn discretized with finite volumes.
    subroutine s_rhs_heatEqn(q_cons_hifu, q_cons_vf, t_step)

        type(scalar_field), dimension(sys_size_hifu), intent(inout) :: q_cons_hifu
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer :: i, j, k, l !< Generic loop iterators
	real(kind(0d0)) :: dTdx_L, dTdx_R, dTdr_L, dTdr_R, val1, val2
        real(kind(0d0)) :: Tx_L, Tx_R, Tr_L, Tr_R
        real(kind(0d0)) :: Ux_L, Ux_R, Ur_L, Ur_R
        real(kind(0d0)) :: Ux_max_L, Ux_max_R, Ur_max_L, Ur_max_R
        real(kind(0d0)) :: Ux_min_L, Ux_min_R, Ur_min_L, Ur_min_R
        integer, intent(in) :: t_step
        integer :: qus_hifu_idx_ht
        real(kind(0d0)) :: absCoef, rho_cp, tdiff
        real(kind(0d0)) :: alpha !volume fraction of each species

        if (hifu_intPrms) then
            qus_hifu_idx_ht = qus_prms_hifu_idx
        else
            qus_hifu_idx_ht = qus_hifu_idx
        end if

        call s_populate_HIFU_variables_buffers(q_cons_hifu)

        l = 0
        do  j = 0, m
            do k = 0, n

                q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) = 0.0d0 !Zeroing RHS_heat 

                !Find temperature derivatives at the faces of the cell
                dTdx_L = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) - q_cons_hifu(T_hifu_idx)%sf(j-1,k,l)) / (x_cc(j)-x_cc(j-1))
                dTdx_R = (q_cons_hifu(T_hifu_idx)%sf(j+1,k,l) - q_cons_hifu(T_hifu_idx)%sf(j,k,l)) / (x_cc(j+1)-x_cc(j))
                dTdr_L = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) - q_cons_hifu(T_hifu_idx)%sf(j,k-1,l)) / (y_cc(k)-y_cc(k-1))
                dTdr_R = (q_cons_hifu(T_hifu_idx)%sf(j,k+1,l) - q_cons_hifu(T_hifu_idx)%sf(j,k,l)) / (y_cc(k+1)-y_cc(k))

                !Find temperature and streaming velocities at the faces of the cell
                Tx_L = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) + q_cons_hifu(T_hifu_idx)%sf(j-1,k,l)) / 2.0d0
                Ux_L = (q_cons_hifu(u_hifu_idx)%sf(j,k,l) + q_cons_hifu(u_hifu_idx)%sf(j-1,k,l)) / 2.0d0

                Tx_R = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) + q_cons_hifu(T_hifu_idx)%sf(j+1,k,l)) / 2.0d0
                Ux_R = (q_cons_hifu(u_hifu_idx)%sf(j,k,l) + q_cons_hifu(u_hifu_idx)%sf(j+1,k,l)) / 2.0d0

                Tr_L = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) + q_cons_hifu(T_hifu_idx)%sf(j,k-1,l)) / 2.0d0
                Ur_L = (q_cons_hifu(v_hifu_idx)%sf(j,k,l) + q_cons_hifu(v_hifu_idx)%sf(j,k-1,l)) / 2.0d0

                Tr_R = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) + q_cons_hifu(T_hifu_idx)%sf(j,k+1,l)) / 2.0d0
                Ur_R = (q_cons_hifu(v_hifu_idx)%sf(j,k,l) + q_cons_hifu(v_hifu_idx)%sf(j,k+1,l)) / 2.0d0

                !Get thermal properties
                alpha  = 0.0d0
                rho_cp = 0.0d0
                tdiff  = 0.0d0
                do i = 1, num_fluids
                    alpha  = q_cons_vf(advxb+i-1)%sf(j,k,l) !q_prim_vf(E_idx+i)%sf(j,k,l)
                    rho_cp = rho_cp + alpha*fluid_pp(i)%rho_cp
                    tdiff  = tdiff  + alpha*fluid_pp(i)%tdiff
                end do

                if ((rho_cp<=0.0) .or. (tdiff<=0.0)) then
                    print*, 'alpha, rho_cp, tdiff', alpha, rho_cp, tdiff
                    call s_mpi_abort('Check thermal properties HIFU!')
                end if

                !Obtain rhs (see notes)
                q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) + &
                       tdiff * (1 / dx(j)) * (dTdx_R - dTdx_L) + &
                       tdiff * (1 / (2*y_cc(k)*dy(k))) * ( (2*y_cc(k)+dy(k))*dTdr_R - (2*y_cc(k)-dy(k))*dTdr_L)
                                
                if ( (q_cons_hifu(tt_hifu_idx)%sf(j,k,l) > 0.0d0) .and. (t_step< hifu_t_step_stopSource)) then !Adding the heat source terms
                    q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) + &
                          (1/(rho_cp)) * (1 / q_cons_hifu(tt_hifu_idx)%sf(j,k,l)) * q_cons_hifu(qus_hifu_idx_ht)%sf(j,k,l) + &  !Acoustic intensity
                          (1/(rho_cp)) * (1 / q_cons_hifu(tt_hifu_idx)%sf(j,k,l)) * q_cons_hifu(qvis_hifu_idx)%sf(j,k,l) + &    !Viscous intensity
                          (1/(rho_cp)) * (1 / q_cons_hifu(tt_hifu_idx)%sf(j,k,l)) * q_cons_hifu(qth_hifu_idx)%sf(j,k,l)         !Thermal intensity
                
                    if (hifu_streaming) then !Convected heat flux
                        q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) - &
                                    (1 / q_cons_hifu(tt_hifu_idx)%sf(j,k,l)) * ( & !Double check this
                                        (1 / dx(j)) * (Ux_R*Tx_R - Ux_L*Tx_L) + &
                                        (1 / (2*y_cc(k)*dy(k))) * ( (2*y_cc(k)+dy(k))*Ur_R*Tr_R - (2*y_cc(k)-dy(k))*Ur_L*Tr_L) )
                    end if

                end if

            end do
        end do
 
    end subroutine s_rhs_heatEqn ! =============================================

    ! ==========================================================================
    !> The purpose of this procedure is to populate the buffers
    !! of the conservative variables, depending on the selected
    !! boundary conditions.
    subroutine s_populate_HIFU_variables_buffers(q_cons_hifu) ! ---------------

        integer :: i, j, k, l, r, q !< Generic loop iterators
        type(scalar_field), dimension(sys_size_hifu), intent(inout) :: q_cons_hifu

        ! Population of Buffers in x-direction =============================

        if (bc_x%beg <= -3) then         ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size_hifu
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_hifu(i)%sf(-j, k, l) = &
                                q_cons_hifu(i)%sf(0, k, l)
                        end do
                    end do
                end do
            end do

        elseif (bc_x%beg == -2) then     ! Symmetry BC at beginning

            do l = 0, p
                do k = 0, n
                    do j = 1, buff_size
                        do i = 1, sys_size_hifu
                            q_cons_hifu(i)%sf(-j, k, l) = &
                                q_cons_hifu(i)%sf(j - 1, k, l)
                        end do
                    end do
                end do
            end do

        elseif (bc_x%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size_hifu
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_hifu(i)%sf(-j, k, l) = &
                                q_cons_hifu(i)%sf(m - (j - 1), k, l)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_variables_buffers( &
                               q_cons_hifu, mpi_dir=1, pbc_loc=-1)
            !call s_mpi_sendrecv_variables_buffers( &
            !    q_cons_hifu, pb, mv, 1, -1)

        end if

        if (bc_x%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size_hifu
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_hifu(i)%sf(m + j, k, l) = &
                                q_cons_hifu(i)%sf(m, k, l)
                        end do
                    end do
                end do
            end do

        elseif (bc_x%end == -2) then     ! Symmetry BC at end

            do l = 0, p
                do k = 0, n
                    do j = 1, buff_size
                        do i = 1, sys_size_hifu
                            q_cons_hifu(i)%sf(m + j, k, l) = &
                                q_cons_hifu(i)%sf(m - (j - 1), k, l)
                        end do
                    end do
                end do
            end do

        elseif (bc_x%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size_hifu
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_hifu(i)%sf(m + j, k, l) = &
                                q_cons_hifu(i)%sf(j - 1, k, l)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_variables_buffers( &
                               q_cons_hifu, mpi_dir=1, pbc_loc=1)
            !call s_mpi_sendrecv_variables_buffers( &
            !    q_cons_hifu, pb, mv,  1, 1)

        end if

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        if (n == 0) then

            return

        elseif (bc_y%beg <= -3 .and. bc_y%beg /= -13) then     ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size_hifu
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(l, -j, k) = &
                                q_cons_hifu(i)%sf(l, 0, k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%beg == -13) then    ! Axis BC at beginning

            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        if (z_cc(k) < pi) then
                            do i = 1, sys_size_hifu
                                q_cons_hifu(i)%sf(l, -j, k) = &
                                    q_cons_hifu(i)%sf(l, j - 1, k + ((p + 1)/2))
                            end do
                        else
                            do i = 1, sys_size_hifu
                                q_cons_hifu(i)%sf(l, -j, k) = &
                                    q_cons_hifu(i)%sf(l, j - 1, k - ((p + 1)/2))
                            end do
                        end if
                    end do
                end do
            end do

        elseif (bc_y%beg == -2) then     ! Symmetry BC at beginning

            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        do i = 1, sys_size_hifu
                            q_cons_hifu(i)%sf(l, -j, k) = &
                                q_cons_hifu(i)%sf(l, j - 1, k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%beg == -1) then     ! Periodic BC at beginning
            do i = 1, sys_size_hifu
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(l, -j, k) = &
                                q_cons_hifu(i)%sf(l, n - (j - 1), k)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_variables_buffers( &
                              q_cons_hifu, mpi_dir=2, pbc_loc=-1)
            !call s_mpi_sendrecv_variables_buffers( &
            !    q_cons_hifu, pb, mv,  2, -1)

        end if

        if (bc_y%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size_hifu
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(l, n + j, k) = &
                                q_cons_hifu(i)%sf(l, n, k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%end == -2) then     ! Symmetry BC at end

            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        do i = 1, sys_size_hifu
                            q_cons_hifu(i)%sf(l, n + j, k) = &
                                q_cons_hifu(i)%sf(l, n - (j - 1), k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size_hifu
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(l, n + j, k) = &
                                q_cons_hifu(i)%sf(l, j - 1, k)
                        end do
                    end do
                end do
            end do      

        else                            ! Processor BC at end

            call s_mpi_sendrecv_variables_buffers( &
                                q_cons_hifu, mpi_dir=2, pbc_loc=1)
            !call s_mpi_sendrecv_variables_buffers( &
            !    q_cons_hifu, pb, mv,  2, 1)

        end if

        ! END: Population of Buffers in y-direction ========================

        ! Population of Buffers in z-direction =============================

        if (p == 0) then

            return

        elseif (bc_z%beg <= -3) then     ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size_hifu
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(k, l, -j) = &
                                q_cons_hifu(i)%sf(k, l, 0)
                        end do
                    end do
                end do
            end do

        elseif (bc_z%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        do i = 1, sys_size_hifu
                            q_cons_hifu(i)%sf(k, l, -j) = &
                                q_cons_hifu(i)%sf(k, l, j - 1)
                        end do
                    end do
                end do
            end do

        elseif (bc_z%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size_hifu
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(k, l, -j) = &
                                q_cons_hifu(i)%sf(k, l, p - (j - 1))
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_variables_buffers( &
                               q_cons_hifu, mpi_dir=3, pbc_loc=-1)
            !call s_mpi_sendrecv_variables_buffers( &
            !    q_cons_hifu, pb, mv,  3, -1)

        end if

        if (bc_z%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size_hifu
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(k, l, p + j) = &
                                q_cons_hifu(i)%sf(k, l, p)
                        end do
                    end do
                end do
            end do

        elseif (bc_z%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        do i = 1, sys_size_hifu
                            q_cons_hifu(i)%sf(k, l, p + j) = &
                                q_cons_hifu(i)%sf(k, l, p - (j - 1))
                        end do
                    end do
                end do
            end do

        elseif (bc_z%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size_hifu
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_hifu(i)%sf(k, l, p + j) = &
                                q_cons_hifu(i)%sf(k, l, j - 1)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at end
        
            call s_mpi_sendrecv_variables_buffers( &
                                q_cons_hifu, mpi_dir=3, pbc_loc=1)
            !call s_mpi_sendrecv_variables_buffers( &
            !    q_cons_hifu, pb, mv,  3, 1)

        end if

        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_HIFU_variables_buffers ! -------------

    subroutine s_open_run_time_information_samplingHIFU()

        character(LEN=path_len + 3*name_len) :: file_path
        integer :: unitFile

        !Open files to save Pmax data at the axial and radial axes
        write (file_path, '(A,I0,A)') '/D/Pmax_', proc_rank, '.dat'
        file_path = trim(case_dir)//trim(file_path)        
        open (100, FILE=trim(file_path), FORM='formatted', STATUS='unknown')
        write (100, *) 'timeStep, x_cc, y_cc, Pmax, Pmin, j, k, l'

        if (proc_rank==0) then

            print*, 'HIFU simulation >>>> Stage 2: Obtaining time-averaged heat source terms, Pmax and Pmin'

            !Open files to save intensity sampling information at focus
            write (file_path, '(A)') '/D/sumIntensity-HIFU.dat'
            file_path = trim(case_dir)//trim(file_path)
            open (99, FILE=trim(file_path), FORM='formatted', POSITION='append', STATUS='unknown')
            write (99, *)   'timeStep, totalSamplingTime, acousticFocalIntensity, acousticFocalIntensityPRMS, ', &
                            'viscousFocalIntensity, thermalFocalIntensity, sumAcousticIntensity, ', &
                            'sumViscousIntensity, sumThermalIntensity, focalxVel, focalyVel'
        

            !Open files to save viscous and thermal intensity sampling information for a single bubble
            write (file_path, '(A,I0,A)') '/D/viscous_thermal_kernel-HIFU_', proc_rank, '.dat'
            file_path = trim(case_dir)//trim(file_path)
            open (98, FILE=trim(file_path), FORM='formatted', POSITION='append', STATUS='unknown')
            write (98, *) 'Recommended to use only with one particle to test and compare the performance of the smootheing function'
            write (98, *) 'Requieres to uncomment some command lines in s_update_RK (m_particle.fpp)'
            write (98, *)   'dt_did, totalSamplingTime, viscousIntensity_beforeKernel, viscousIntensity_afterKernel, ', &
                            'thermalIntensity_beforeKernel, thermalIntensity_afterKernel, radius, velocity'

        end if
        

    end subroutine s_open_run_time_information_samplingHIFU

    subroutine s_close_run_time_information_samplingHIFU()

        integer :: unitFile

        !Close files to save Pmax data at the axial and radial axes
        close(100)

        !Close file to save intensity sampling information at focus
        if (proc_rank==0) close (99)

        !Close file to save viscous and thermal intensity sampling information for a single bubble
        close(98)

    end subroutine s_close_run_time_information_samplingHIFU

end module
