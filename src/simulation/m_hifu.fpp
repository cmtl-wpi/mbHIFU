!>
!! @file m_hifu.f90
!! @brief Contains module m_viscous

#:include 'macros.fpp'
#:include 'inline_conversions.fpp'

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
    private; public ::  s_restart_Pmax, &
                        s_update_Pmax, &
                        s_update_intensity_HIFU, &
                        s_initialize_HIFU, &
                        s_rhs_heatEqn, &
                        s_populate_HIFU_variables_buffers

contains

    ! ==========================================================================
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! READ ME: !!!!!!!!
    ! Important variables and where they are stored:
    ! q_cons_hifu(1)%sf(j,k,l): Heat intensity (generated). Summation over time (finite number of samples).
    ! q_cons_hifu(2)%sf(j,k,l): Number of finite samples to compute q_cons_hifu(1).
    ! q_cons_hifu(3)%sf(j,k,l): Temperature.
    ! q_cons_hifu(4)%sf(j,k,l): RHS value from heat transfer eqn (finite volume discretization).
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! NEW VERSION
    ! Important variables and where they are stored:
    ! q_cons_hifu(1)%sf(j,k,l): Temperature distribution                (T_hifu_idx)
    ! q_cons_hifu(2)%sf(j,k,l): RHS value from heat transfer eqn        (T_hifu_idx+1)
    ! q_cons_hifu(3)%sf(j,k,l): Maximum pressure at each cell           (P_hifu_idx)
    ! q_cons_hifu(4)%sf(j,k,l): Number of samples to avg intensities    (N_hifu_idx)
    ! q_cons_hifu(5)%sf(j,k,l): Acoustic damping (Sum over time)        (qus_hifu_idx)
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

        !do i = ix_t%beg, ix_t%end
        !    do j = iy_t%beg, iy_t%end
        !        q_cons_hifu(dmb_hifu_idx)%sf(i,j,0) = (i+j)*1.0d0
        !    end do
        !end do

        !Initial Temperature
        !if (proc_rank==0) print*,T_hifu_idx, size(q_cons_hifu)
        q_cons_hifu(T_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
                                   iz_t%beg:iz_t%end) = hifu_Tref !User input

        !Initialize Pmax
        q_cons_hifu(P_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
                                   iz_t%beg:iz_t%end) =  min(dflt_real,-dflt_real)

        !Initialize umax, vmax
        !q_cons_hifu(umax_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
        !                           iz_t%beg:iz_t%end) =  min(dflt_real,-dflt_real)

        !q_cons_hifu(vmax_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
        !                           iz_t%beg:iz_t%end) =  min(dflt_real,-dflt_real)

        !Initialize umin, vmin
        !q_cons_hifu(umin_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
        !                           iz_t%beg:iz_t%end) =  max(dflt_real,-dflt_real)

        !q_cons_hifu(vmin_hifu_idx)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, &
        !                           iz_t%beg:iz_t%end) =  max(dflt_real,-dflt_real)

    end subroutine s_initialize_HIFU ! =========================================

    ! ==========================================================================
    ! Restarting or initializing routine to get the maximum and minimum pressure 
    ! through time along the axisymmetric and radial axes
    subroutine s_restart_Pmax()

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist
        real(kind(0d0)) :: trsh1, trsh2, trsh3, presMax, presMin
        integer :: i,j, k
        integer :: unitFile
        logical :: printFlag = .true.

        ! HIFU Open Pmax file if it exists
        write (file_loc, '(A,I0,A)') '/D/Pmax_', proc_rank, '.dat'
        file_loc = trim(case_dir)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (proc_rank==0 .and. printFlag) then
             printFlag = .false.
             print*, 'dx, dy, dt', dx(1), dy(1), dt
        end if

        unitFile = 100+proc_rank

        if (file_exist) then
             !open (unitFile, FILE=trim(file_loc), FORM='formatted', STATUS='old', &
             !                                   access='sequential', action='read')
             !334 read(unitFile,'(6x,f12.6,f24.8,f24.8,f24.8,f24.8,I24,I24,I24)',end=335) &
             !        trsh1, &
             !        trsh2, &
             !        trsh3, &
             !        presMax, &
             !        presMin, &
             !        i, &
             !        j, &
             !        k
                
             !Pmax(i,j)=presMax
             !Pmin(i,j)=presMin
             !goto 334
             !335 continue
             !close(unitFile)
        else
            !print*, 'Are we here?', dflt_real, -dflt_real
            !do i=0,m
            !    do j=0,n
                    !Pmax(i,j)= min(dflt_real,-dflt_real)
                    !Pmin(i,j)= max(dflt_real,-dflt_real)
            !    end do
            !end do
            !print*, 'Done generating initial Pmax and Pmin'

        end if

    end subroutine s_restart_Pmax ! ============================================

    ! ==========================================================================
    ! Update the maximum and minimum pressure through time along the 
    ! axisymmetric and radial axes
    subroutine s_update_Pmax(q_cons_vf,q_cons_hifu,t_step)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        type(scalar_field), dimension(sys_size_hifu), intent(INOUT) :: q_cons_hifu
        integer, intent(IN) :: t_step
        
        real(kind(0d0)) :: nondim_time !< Non-dimensional time	
        character(LEN=path_len + 3*name_len) :: file_path !<
        logical :: axialCondition, radialCondition, condition

        integer :: i, j, k, l, s, q  !< generic loop variables
        integer :: unitFile

	real(kind(0d0)) :: rho
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: pres
        real(kind(0d0)) :: presMax_old, presMin_old
	real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: qv
	real(kind(0d0)), dimension(2) :: Re
	real(kind(0d0)) :: G

        ! Zeroing out flow variables for all processors      
        rho = 0d0
        do s = 1, num_dims
           vel(s) = 0d0
        end do
        pres = 0d0
        presMax_old = 0d0
        presMin_old =0d0
        gamma = 0d0
        pi_inf = 0d0

        unitFile = 100+proc_rank

        if (cyl_coord .and. p==0) then
            l = 0
            if (mod(t_step,t_step_save)==0) then
                write (file_path, '(A,I0,A)') '/D/Pmax_', proc_rank, '.dat'
                file_path = trim(case_dir)//trim(file_path)        
                open (unitFile, FILE=trim(file_path), FORM='formatted', STATUS='unknown')
            end if

            do  j = 0, m
                do k = 0, n
                    call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                            rho, gamma, pi_inf, qv, &
                                                            Re, G, fluid_pp(:)%G)
                    do s = 1, num_dims
                        vel(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho
                    end do

                    call s_compute_pressure(q_cons_vf(E_idx)%sf(j, k, l), &
                                0d0, 0.5d0*rho*dot_product(vel, vel), pi_inf, gamma, rho, qv, pres)

                    !presMax_old = Pmax(j,k)
                    !presMin_old = Pmin(j,k)
                    !Pmax(j,k) = max(pres,presMax_old)
                    !Pmin(j,k) = min(pres,presMin_old)

                    q_cons_hifu(P_hifu_idx)%sf(j,k,l)=max(q_cons_hifu(P_hifu_idx)%sf(j,k,l),pres)
                                      
                    ! Specify enough conditions for axial and radial probe lines
                    axialCondition= (dy(k)>y_cc(k) .and. y_cc(k)>0)
                    radialCondition= (x_cb(j-1)<mono(1)%foc_length .and. mono(1)%foc_length<x_cb(j))
                    condition= (axialCondition .or. radialCondition)

                    if (mod(t_step,t_step_save)==0 .and. condition) then
                        write (unitFile, '(6x,I24,f24.8,f24.8,f24.8,f24.8,I24,I24,I24)') &
                            t_step, &
                            x_cc(j), &
                            y_cc(k), &
                            q_cons_hifu(P_hifu_idx)%sf(j,k,l), &
                            pres, &
                            j, &
                            k, &
                            l
                    end if

                end do
            end do

            if (mod(t_step,t_step_save)==0) then
                close(unitFile)
            end if

        end if
   
    end subroutine s_update_Pmax ! =============================================

    @:s_compute_speed_of_sound()

    ! ==========================================================================
    !> Heat deposition in HIFU.
    !> Obtain the generated heat source "q_us_ac". Heating from the primary ultrasound source
    !> and, if mb exist, from acoustic emission from bubble oscillations.
    subroutine s_update_intensity_HIFU(q_cons_vf, q_prim_vf, t_step, q_cons_hifu)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        type(scalar_field), dimension(sys_size_hifu), intent(INOUT) :: q_cons_hifu
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        integer, intent(IN) :: t_step

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

	real(kind(0d0)) :: shearVisc, bulkVisc
	real(kind(0d0)) :: absCoef, spdsound 
	real(kind(0d0)) :: angFreq
	real(kind(0d0)) :: varA, varB
	real(kind(0d0)) :: duxdx, duxdr, durdx, durdr
	real(kind(0d0)) :: ep11, ep22, ep33, ep13
	real(kind(0d0)) :: intensity_ac, sumIntensity_ac, tmp, focalIntensity_ac
        real(kind(0d0)) :: intensity_vis, sumIntensity_vis, focalIntensity_vis
        real(kind(0d0)) :: focal_u, focal_v

	real(kind(0d0)), dimension(num_fluids) :: myalpha_rho, myalpha
        real(kind(0d0)) :: n_tait, B_tait, myRho, lamda

        integer :: ii, i, j, k, l, s, q  !< generic loop variables

        character(LEN=path_len + 3*name_len) :: file_path !<
	logical :: file_exist, printFlag1, printFlag2
        logical :: axialCondition, radialCondition, condition
        real(kind(0d0)) :: focalPoint_x, val, dist
        
        
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
        focalIntensity_vis = 0d0
        focalPoint_x = mono(1)%foc_length + mono(1)%loc(1)
        sumIntensity_ac = 0d0
        sumIntensity_vis = 0d0
        focal_u = 0.0d0
        focal_v = 0.0d0

        if (t_step==t_step_start) then
            stepsPerWave = int(mono(1)%length / dt) !Fix me if bubbles present
            if (proc_rank==0) print*, 'Every wave is discretized in ', mono(1)%length/dt, 'time steps (', stepsPerWave,')'
        end if

        if (cyl_coord .and. p==0) then  !Axysimetric		
            l = 0

            if (proc_rank==0 .and. t_step_old==dflt_int) then
                write (file_path, '(A,I0,A)') '/D/sumIntensity-HIFU.dat'
                file_path = trim(case_dir)//trim(file_path)
                open (99, FILE=trim(file_path), FORM='formatted', POSITION='append', STATUS='unknown')
            end if

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
                                0d0, 0.5d0*rho*dot_product(vel, vel), pi_inf, gamma, rho, qv, pres)

                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                      ((gamma + 1d0)*pres + pi_inf)/rho, alpha, 0d0, c)

                    !>> Compute intensity form acoustic damping
                    intensity_ac = 0d0
                    
                    if (hifu_intPrms) then
                        if (absCoef<=0.0) call s_mpi_abort('Prms: check absCoef values!')
                        intensity_ac = intensity_ac + absCoef*(q_cons_hifu(P_hifu_idx)%sf(j,k,l)-hifu_atmPres)**2/(rho*c)
                    else
                        ep11 = durdr
                        ep22 = vel(2) / y_cc(k)
                        ep33 = duxdx
                        ep13 = 0.5d0*(durdx + duxdr)
                        varA = ep11**2.0 + ep22**2.0 + ep33**2.0
                        varB = (8.0d0/3.0d0)*varA - (4.0d0/3.0d0)*(ep11*ep22 + ep11*ep33 + ep22*ep33) + 6.0d0*(ep13**2.0)
                        intensity_ac = intensity_ac + bulkVisc*varA + 2.0d0*shearVisc*varB !intensity is "q_us_ac"
                        
                    end if

                    q_cons_hifu(qus_hifu_idx)%sf(j,k,l) = q_cons_hifu(qus_hifu_idx)%sf(j,k,l) + intensity_ac ! sum intensity to old value
                    q_cons_hifu(N_hifu_idx)%sf(j,k,l) = q_cons_hifu(N_hifu_idx)%sf(j,k,l) + 1.0d0 !number of samples taken

                    !Update average velocities for streaming
                    !if (stepStreaming == 1) then
                    !    q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) = vel(1)
                    !    q_cons_hifu(u_hifu_idx+2)%sf(j,k,l) = vel(1)
                    !    q_cons_hifu(v_hifu_idx+1)%sf(j,k,l) = vel(2)
                    !    q_cons_hifu(v_hifu_idx+2)%sf(j,k,l) = vel(2)
                    !else if (stepStreaming == stepsPerWave) then
                    !    q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) + vel(1)
                    !    q_cons_hifu(u_hifu_idx)%sf(j,k,l)   = q_cons_hifu(u_hifu_idx)%sf(j,k,l) + ( q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) - &
                    !                                          q_cons_hifu(u_hifu_idx+2)%sf(j,k,l) ) * (1.0d0 / stepsPerWave)
                    !    q_cons_hifu(v_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(v_hifu_idx+1)%sf(j,k,l) + vel(2)
                    !    q_cons_hifu(v_hifu_idx)%sf(j,k,l)   = q_cons_hifu(v_hifu_idx)%sf(j,k,l) + ( q_cons_hifu(v_hifu_idx+1)%sf(j,k,l) - &
                    !                                          q_cons_hifu(v_hifu_idx+2)%sf(j,k,l) ) * (1.0d0 / stepsPerWave)
                    !else
                    !    q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) + vel(1)
                    !    q_cons_hifu(v_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(v_hifu_idx+1)%sf(j,k,l) + vel(2)
                    !end if
                    if (hifu_streaming) then
                         q_cons_hifu(u_hifu_idx)%sf(j,k,l) =  q_cons_hifu(u_hifu_idx)%sf(j,k,l) + vel(1)
                         q_cons_hifu(v_hifu_idx)%sf(j,k,l) =  q_cons_hifu(v_hifu_idx)%sf(j,k,l) + vel(2)
                    !    if (q_cons_hifu(N_hifu_idx)%sf(j,k,l) == 1.0d0) then
                    !        q_cons_hifu(u_hifu_idx)%sf(j,k,l) = vel(1)          !Avg velocity 
                    !        q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) = vel(1)        !Old slope
                    !        q_cons_hifu(u_hifu_idx+2)%sf(j,k,l) = vel(1)        !Vel Peak on avg velocity plot
                    !        q_cons_hifu(u_hifu_idx+3)%sf(j,k,l) = vel(1)        !Vel Peak on avg velocity plot
                    !    else
                    !        q_cons_hifu(u_hifu_idx)%sf(j,k,l) = q_cons_hifu(u_hifu_idx)%sf(j,k,l) * &
                    !                                            (q_cons_hifu(N_hifu_idx)%sf(j,k,l)-1.0d0) + vel(1)
                    !        q_cons_hifu(u_hifu_idx)%sf(j,k,l) = q_cons_hifu(u_hifu_idx)%sf(j,k,l) / q_cons_hifu(N_hifu_idx)%sf(j,k,l)
                    !        q_cons_hifu(u_hifu_idx+1)%sf(j,k,l) = vel(1)
                    !    end if
                    end if


                    !Get focal intensity and velocities
                    axialCondition= (dy(k)>y_cc(k) .and. y_cc(k)>0.0)
                    radialCondition= (x_cb(j-1)<=focalPoint_x .and. focalPoint_x<=x_cb(j))
                    condition= (axialCondition .and. radialCondition)
                    if (condition) then
                            focalIntensity_ac = q_cons_hifu(qus_hifu_idx)%sf(j,k,l)
                            focalIntensity_vis = q_cons_hifu(qvis_hifu_idx)%sf(j,k,l)
                            focal_u = q_cons_hifu(u_hifu_idx)%sf(j,k,l)
                            focal_v = q_cons_hifu(v_hifu_idx)%sf(j,k,l)
                            !vel_u = vel(1)
                            !vel_v = vel(2)
                    end if

                    !Intensity summation through the domain, avoid monopole influence (0.8*Focal length)
                    sumIntensity_ac = sumIntensity_ac + q_cons_hifu(qus_hifu_idx)%sf(j,k,l)
                    sumIntensity_vis = sumIntensity_vis + q_cons_hifu(qvis_hifu_idx)%sf(j,k,l)

                end do
            end do

            stepStreaming = stepStreaming + 1

            if (stepStreaming > stepsPerWave) stepStreaming = 1

            tmp = sumIntensity_ac
            call s_mpi_allreduce_sum(tmp, sumIntensity_ac)
            tmp = sumIntensity_vis
            call s_mpi_allreduce_sum(tmp, sumIntensity_vis)

            tmp = focalIntensity_vis
            call s_mpi_allreduce_max(tmp, focalIntensity_vis)
            tmp = focalIntensity_ac
            call s_mpi_allreduce_max(tmp, focalIntensity_ac)

            tmp = focal_u
            call s_mpi_allreduce_sum(tmp, focal_u)
            tmp = focal_v
            call s_mpi_allreduce_sum(tmp, focal_v)

            !tmp = vel_u
            !call s_mpi_allreduce_sum(tmp, vel_u)
            !tmp = vel_v
            !call s_mpi_allreduce_sum(tmp, vel_v)

            if (proc_rank==0) write (99, '(6x,I24.8,f24.8,e24.8,e24.8,e24.8,e24.8,e24.8,e24.8)') &
                t_step, q_cons_hifu(N_hifu_idx)%sf(0,0,0), focalIntensity_ac, focalIntensity_vis, &
                        sumIntensity_ac, sumIntensity_vis, focal_u, focal_v

            if (t_step == t_step_stop .and. proc_rank==0) close (99)
  
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
        real(kind(0d0)) :: absCoef, rho_cp, tdiff
        real(kind(0d0)) :: alpha !volume fraction of each species


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
                !Ux_max_L = (q_cons_hifu(umax_hifu_idx)%sf(j,k,l) + q_cons_hifu(umax_hifu_idx)%sf(j-1,k,l)) / 2.0d0
                !Ux_min_L = (q_cons_hifu(umin_hifu_idx)%sf(j,k,l) + q_cons_hifu(umin_hifu_idx)%sf(j-1,k,l)) / 2.0d0
                !Ux_L = Ux_max_L-0.5d0*(Ux_max_L-Ux_min_L)

                Tx_R = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) + q_cons_hifu(T_hifu_idx)%sf(j+1,k,l)) / 2.0d0
                Ux_R = (q_cons_hifu(u_hifu_idx)%sf(j,k,l) + q_cons_hifu(u_hifu_idx)%sf(j+1,k,l)) / 2.0d0
                !Ux_max_R = (q_cons_hifu(umax_hifu_idx)%sf(j,k,l) + q_cons_hifu(umax_hifu_idx)%sf(j+1,k,l)) / 2.0d0
                !Ux_min_R = (q_cons_hifu(umin_hifu_idx)%sf(j,k,l) + q_cons_hifu(umin_hifu_idx)%sf(j+1,k,l)) / 2.0d0
                !Ux_R = Ux_max_R-0.5d0*(Ux_max_R-Ux_min_R)

                Tr_L = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) + q_cons_hifu(T_hifu_idx)%sf(j,k-1,l)) / 2.0d0
                Ur_L = (q_cons_hifu(v_hifu_idx)%sf(j,k,l) + q_cons_hifu(v_hifu_idx)%sf(j,k-1,l)) / 2.0d0
                !Ur_max_L = (q_cons_hifu(vmax_hifu_idx)%sf(j,k,l) + q_cons_hifu(vmax_hifu_idx)%sf(j,k-1,l)) / 2.0d0
                !Ur_min_L = (q_cons_hifu(vmin_hifu_idx)%sf(j,k,l) + q_cons_hifu(vmin_hifu_idx)%sf(j,k-1,l)) / 2.0d0
                !Ur_L = Ur_max_L-0.5d0*(Ur_max_L-Ur_min_L)

                Tr_R = (q_cons_hifu(T_hifu_idx)%sf(j,k,l) + q_cons_hifu(T_hifu_idx)%sf(j,k+1,l)) / 2.0d0
                Ur_R = (q_cons_hifu(v_hifu_idx)%sf(j,k,l) + q_cons_hifu(v_hifu_idx)%sf(j,k+1,l)) / 2.0d0
                !Ur_max_R = (q_cons_hifu(vmax_hifu_idx)%sf(j,k,l) + q_cons_hifu(vmax_hifu_idx)%sf(j,k+1,l)) / 2.0d0
                !Ur_min_R = (q_cons_hifu(vmin_hifu_idx)%sf(j,k,l) + q_cons_hifu(vmin_hifu_idx)%sf(j,k+1,l)) / 2.0d0
                !Ur_R = Ur_max_R-0.5d0*(Ur_max_R-Ur_min_R)

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
                                
                if ( (q_cons_hifu(N_hifu_idx)%sf(j,k,l) > 0.0) .and. (t_step< hifu_t_step_stopSource)) then !Adding the heat source terms
                    q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) + &
                          (1/(rho_cp)) * (1 / q_cons_hifu(N_hifu_idx)%sf(j,k,l)) * q_cons_hifu(qus_hifu_idx)%sf(j,k,l) + &
                          (1/(rho_cp)) * (1 / q_cons_hifu(N_hifu_idx)%sf(j,k,l)) * q_cons_hifu(qvis_hifu_idx)%sf(j,k,l)
                
                    if (hifu_streaming) then !Convected heat flux
                        q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) = q_cons_hifu(T_hifu_idx+1)%sf(j,k,l) - &
                                  ( (1 / dx(j)) * (Ux_R*Tx_R - Ux_L*Tx_L) + &
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

            call s_mpi_sendrecv_conservative_variables_buffers( &
                               q_cons_hifu, mpi_dir=1, pbc_loc=-1)
            !call s_mpi_sendrecv_conservative_variables_buffers( &
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

            call s_mpi_sendrecv_conservative_variables_buffers( &
                               q_cons_hifu, mpi_dir=1, pbc_loc=1)
            !call s_mpi_sendrecv_conservative_variables_buffers( &
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

            call s_mpi_sendrecv_conservative_variables_buffers( &
                              q_cons_hifu, mpi_dir=2, pbc_loc=-1)
            !call s_mpi_sendrecv_conservative_variables_buffers( &
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

            call s_mpi_sendrecv_conservative_variables_buffers( &
                                q_cons_hifu, mpi_dir=2, pbc_loc=1)
            !call s_mpi_sendrecv_conservative_variables_buffers( &
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

            call s_mpi_sendrecv_conservative_variables_buffers( &
                               q_cons_hifu, mpi_dir=3, pbc_loc=-1)
            !call s_mpi_sendrecv_conservative_variables_buffers( &
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
        
            call s_mpi_sendrecv_conservative_variables_buffers( &
                                q_cons_hifu, mpi_dir=3, pbc_loc=1)
            !call s_mpi_sendrecv_conservative_variables_buffers( &
            !    q_cons_hifu, pb, mv,  3, 1)

        end if

        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_HIFU_variables_buffers ! -------------

end module
