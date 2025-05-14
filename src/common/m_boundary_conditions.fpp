!>
!! @file m_boundary_conditions.fpp
!! @brief Contains module m_boundary_conditions

!> @brief The purpose of the module is to apply noncharacteristic and processor
!! boundary condiitons
module m_boundary_conditions

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy

    use m_constants

    implicit none

#ifdef MFC_SIMULATION
    private; public :: s_populate_variables_buffers, s_populate_capillary_buffers
#else
    private; public :: s_populate_variables_buffers
#endif

contains

    !>  The purpose of this procedure is to populate the buffers
        !!      of the primitive variables, depending on the selected
        !!      boundary conditions.
    subroutine s_populate_variables_buffers(q_prim_vf, pb, mv)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf

        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv

        integer :: bc_loc, bc_dir

        ! Population of Buffers in x-direction

        select case (bc_x%beg)
        case (-13:-3) ! Ghost-cell extrap. BC at beginning
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, -1)
        case (-2)     ! Symmetry BC at beginning
            call s_symmetry(q_prim_vf, pb, mv, 1, -1)
        case (-1)     ! Periodic BC at beginning
            call s_periodic(q_prim_vf, pb, mv, 1, -1)
        case (-22)    ! Periodic BC at beginning (translational)
            call s_periodic_rotational(q_prim_vf, pb, mv, 1, -1)
        case (-15)    ! Slip wall BC at beginning
            call s_slip_wall(q_prim_vf, pb, mv, 1, -1)
        case (-16)    ! No-slip wall BC at beginning
            call s_no_slip_wall(q_prim_vf, pb, mv, 1, -1)
        case (-20)    ! Sinusoudal pressure (acoustic transducer)
            call s_acoustic_bc(q_prim_vf, pb, mv, 1, -1)
        case default ! Processor BC at beginning
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 1, -1)
        end select

        select case (bc_x%end)
        case (-13:-3) ! Ghost-cell extrap. BC at end
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, 1)
        case (-2)     ! Symmetry BC at end
            call s_symmetry(q_prim_vf, pb, mv, 1, 1)
        case (-1)     ! Periodic BC at end
            call s_periodic(q_prim_vf, pb, mv, 1, 1)
        case (-22)    ! Periodic BC at beginning (translational)
            call s_periodic_rotational(q_prim_vf, pb, mv, 1, 1)
        case (-15)    ! Slip wall BC at end
            call s_slip_wall(q_prim_vf, pb, mv, 1, 1)
        case (-16)    ! No-slip wall bc at end
            call s_no_slip_wall(q_prim_vf, pb, mv, 1, 1)
        case default ! Processor BC at end
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 1, 1)
        end select

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            select case (bc_x%beg)
            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                call s_qbmm_extrapolation(pb, mv, 1, -1)
            case (-15)    ! Slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 1, -1)
            case (-16)    ! No-slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 1, -1)
            end select

            select case (bc_x%end)
            case (-13:-3) ! Ghost-cell extrap. BC at end
                call s_qbmm_extrapolation(pb, mv, 1, 1)
            case (-15)    ! Slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 1, 1)
            case (-16)    ! No-slip wall bc at end
                call s_qbmm_extrapolation(pb, mv, 1, 1)
            end select
        end if
#endif

        ! Population of Buffers in y-direction

        if (n == 0) return

        select case (bc_y%beg)
        case (-13:-3) ! Ghost-cell extrap. BC at beginning
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, -1)
        case (-14)    ! Axis BC at beginning
            call s_axis(q_prim_vf, pb, mv, 2, -1)
        case (-2)     ! Symmetry BC at beginning
            call s_symmetry(q_prim_vf, pb, mv, 2, -1)
        case (-1)     ! Periodic BC at beginning
            call s_periodic(q_prim_vf, pb, mv, 2, -1)
        case (-22)    ! Periodic BC at beginning (translational)
            call s_periodic_rotational(q_prim_vf, pb, mv, 2, -1)
        case (-15)    ! Slip wall BC at beginning
            call s_slip_wall(q_prim_vf, pb, mv, 2, -1)
        case (-16)    ! No-slip wall BC at beginning
            call s_no_slip_wall(q_prim_vf, pb, mv, 2, -1)
        case (-20)    ! Sinusoudal pressure (acoustic transducer)
            call s_acoustic_bc(q_prim_vf, pb, mv, 2, -1)
        case (-21)    ! Axis BC in a cylindrical sector HIFU
            call s_axis_cylindrical_sector_hifu(q_prim_vf, pb, mv, 2, -1)
        case default ! Processor BC at beginning
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 2, -1)
        end select

        select case (bc_y%end)
        case (-13:-3) ! Ghost-cell extrap. BC at end
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, 1)
        case (-2)     ! Symmetry BC at end
            call s_symmetry(q_prim_vf, pb, mv, 2, 1)
        case (-1)     ! Periodic BC at end
            call s_periodic(q_prim_vf, pb, mv, 2, 1)
        case (-22)    ! Periodic BC at beginning (translational)
            call s_periodic_rotational(q_prim_vf, pb, mv, 2, 1)
        case (-15)    ! Slip wall BC at end
            call s_slip_wall(q_prim_vf, pb, mv, 2, 1)
        case (-16)    ! No-slip wall BC at end
            call s_no_slip_wall(q_prim_vf, pb, mv, 2, 1)
        case default ! Processor BC at end
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 2, 1)
        end select

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            select case (bc_y%beg)
            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                call s_qbmm_extrapolation(pb, mv, 2, -1)
            case (-15)    ! Slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 2, -1)
            case (-16)    ! No-slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 2, -1)
            end select

            select case (bc_y%end)
            case (-13:-3) ! Ghost-cell extrap. BC at end
                call s_qbmm_extrapolation(pb, mv, 2, 1)
            case (-15)    ! Slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 2, 1)
            case (-16)    ! No-slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 2, 1)
            end select
        end if
#endif

        ! Population of Buffers in z-direction

        if (p == 0) return

        select case (bc_z%beg)
        case (-13:-3) ! Ghost-cell extrap. BC at beginning
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, -1)
        case (-2)     ! Symmetry BC at beginning
            call s_symmetry(q_prim_vf, pb, mv, 3, -1)
        case (-1)     ! Periodic BC at beginning
            call s_periodic(q_prim_vf, pb, mv, 3, -1)
        case (-22)    ! Periodic BC at beginning (translational)
            call s_periodic_rotational(q_prim_vf, pb, mv, 3, -1)
        case (-15)    ! Slip wall BC at beginning
            call s_slip_wall(q_prim_vf, pb, mv, 3, -1)
        case (-16)    ! No-slip wall BC at beginning
            call s_no_slip_wall(q_prim_vf, pb, mv, 3, -1)
        case (-20)    ! Sinusoudal pressure (acoustic transducer)
            call s_acoustic_bc(q_prim_vf, pb, mv, 3, -1)
        case default ! Processor BC at beginning
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 3, -1)
        end select

        select case (bc_z%end)
        case (-13:-3) ! Ghost-cell extrap. BC at end
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, 1)
        case (-2)     ! Symmetry BC at end
            call s_symmetry(q_prim_vf, pb, mv, 3, 1)
        case (-1)     ! Periodic BC at end
            call s_periodic(q_prim_vf, pb, mv, 3, 1)
        case (-22)    ! Periodic BC at beginning (translational)
            call s_periodic_rotational(q_prim_vf, pb, mv, 3, 1)
        case (-15)    ! Slip wall BC at end
            call s_slip_wall(q_prim_vf, pb, mv, 3, 1)
        case (-16)    ! No-slip wall BC at end
            call s_no_slip_wall(q_prim_vf, pb, mv, 3, 1)
        case default ! Processor BC at end
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 3, 1)
        end select

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            select case (bc_z%beg)
            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                call s_qbmm_extrapolation(pb, mv, 3, -1)
            case (-15)    ! Slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 3, -1)
            case (-16)    ! No-slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 3, -1)
            end select

            select case (bc_z%end)
            case (-13:-3) ! Ghost-cell extrap. BC at end
                call s_qbmm_extrapolation(pb, mv, 3, 1)
            case (-15)    ! Slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 3, 1)
            case (-16)    ! No-slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 3, 1)
            end select
        end if
#endif
        ! END: Population of Buffers in z-direction

    end subroutine s_populate_variables_buffers

    subroutine s_acoustic_bc(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer :: j, k, l, q, i
        real(wp) :: tau, gFun, rc, rbeta, radial_cc

        !< x-direction =========================================================
        if (bc_dir == 1) then !< x-direction

            if (bc_loc == -1) then !bc_x%beg

                !$acc parallel loop collapse(3) gang vector default(present) copyin(mytime)
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            tau = 0._wp

                            ! Velocities (# dim), and fluids' volume fraction
                            !$acc loop seq
                            do i = 1, sys_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(0, k, l)
                            end do

                            if (acoustic_bc_params%iwave == 1) then
                                ! Pressure : Planar wave
                                tau = mytime

                                if (tau < (acoustic_bc_params%ncycles/acoustic_bc_params%freq)) then
                                    q_prim_vf(momxe + 1)%sf(-j, k, l) = acoustic_bc_params%Pbase + &
                                                                        acoustic_bc_params%Pamp* &
                                                                        sin(2._wp*pi*acoustic_bc_params%freq*tau)
                                    q_prim_vf(1)%sf(-j, k, l) = acoustic_bc_params%rho
                                else
                                    q_prim_vf(momxe + 1)%sf(-j, k, l) = acoustic_bc_params%Pbase
                                end if

                            elseif (acoustic_bc_params%iwave == 2) then
                                ! Single hemispherical transdurer
                                rc = (0.5_wp*acoustic_bc_params%apert)/sqrt(1._wp - ((0.5_wp*acoustic_bc_params%apert)/ &
                                                                                     (acoustic_bc_params%focLen + acoustic_bc_params%focCal))**2._wp)
                                rbeta = sqrt(1._wp + ((0.5_wp*acoustic_bc_params%apert)/ &
                                                      (acoustic_bc_params%focLen + acoustic_bc_params%focCal))**2._wp)
                                radial_cc = y_cc(k)
                                if (p > 0) radial_cc = sqrt(y_cc(k)**2._wp + z_cc(l)**2._wp)

                                if (radial_cc < rc) then
                                    tau = mytime + radial_cc**2._wp/(2._wp*acoustic_bc_params%cson* &
                                                                     (acoustic_bc_params%focLen + acoustic_bc_params%focCal))
                                    gFun = (1._wp/rbeta)

                                    if (tau < (acoustic_bc_params%ncycles/acoustic_bc_params%freq)) then
                                        q_prim_vf(momxe + 1)%sf(-j, k, l) = acoustic_bc_params%Pbase + &
                                                                            acoustic_bc_params%Pamp* &
                                                                            sin(2._wp*pi*acoustic_bc_params%freq*tau)
                                        q_prim_vf(1)%sf(-j, k, l) = acoustic_bc_params%rho
                                    else
                                        q_prim_vf(momxe + 1)%sf(-j, k, l) = acoustic_bc_params%Pbase
                                    end if
                                else
                                    q_prim_vf(momxe + 1)%sf(-j, k, l) = acoustic_bc_params%Pbase
                                end if

                            else
                                stop "acoustic_bc_params%iwave incorrect value (1: planar wave, 2: axisymmetric spherical transducer)."
                            end if

                        end do
                    end do
                end do

            else !< bc_x%end

                call s_mpi_abort('Transucer not available for end x-dir')

            end if

            !< y-direction =========================================================
        elseif (bc_dir == 2) then !< y-direction

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(3) gang vector default(present) copyin(mytime)
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            tau = 0._wp

                            ! Velocities (# dim), and fluids' volume fraction
                            !$acc loop seq
                            do i = 1, sys_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, 0, k)
                            end do

                            if (acoustic_bc_params%iwave == 1) then
                                ! Pressure : Planar wave
                                tau = mytime

                                if (tau < (acoustic_bc_params%ncycles/acoustic_bc_params%freq)) then
                                    q_prim_vf(momxe + 1)%sf(l, -j, k) = acoustic_bc_params%Pbase + &
                                                                        acoustic_bc_params%Pamp* &
                                                                        sin(2._wp*pi*acoustic_bc_params%freq*tau)
                                    q_prim_vf(1)%sf(l, -j, k) = acoustic_bc_params%rho
                                else
                                    q_prim_vf(momxe + 1)%sf(l, -j, k) = acoustic_bc_params%Pbase
                                end if

                            elseif (acoustic_bc_params%iwave == 2) then
                                ! Single hemispherical transdurer
                                stop "Axisymmetric spherical transducer only valid with bc_x%beg."
                            else
                                stop "acoustic_bc_params%iwave incorrect value (1: planar wave, 2: axisymmetric spherical transducer)."
                            end if

                        end do
                    end do
                end do

            else !< bc_y%end

                call s_mpi_abort('Transucer not available for y-dir')

            end if

            !< z-direction =========================================================
        elseif (bc_dir == 3) then !< z-direction

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(3) gang vector default(present) copyin(mytime)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            tau = 0._wp

                            ! Velocities (# dim), and fluids' volume fraction
                            !$acc loop seq
                            do i = 1, sys_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, 0)
                            end do

                            if (acoustic_bc_params%iwave == 1) then
                                ! Pressure : Planar wave
                                tau = mytime

                                if (tau < (acoustic_bc_params%ncycles/acoustic_bc_params%freq)) then
                                    q_prim_vf(momxe + 1)%sf(k, l, -j) = acoustic_bc_params%Pbase + &
                                                                        acoustic_bc_params%Pamp* &
                                                                        sin(2._wp*pi*acoustic_bc_params%freq*tau)
                                    q_prim_vf(1)%sf(k, l, -j) = acoustic_bc_params%rho
                                else
                                    q_prim_vf(momxe + 1)%sf(k, l, -j) = acoustic_bc_params%Pbase
                                end if

                            elseif (acoustic_bc_params%iwave == 2) then
                                ! Single hemispherical transdurer
                                rc = (0.5_wp*acoustic_bc_params%apert)/sqrt(1._wp - ((0.5_wp*acoustic_bc_params%apert)/ &
                                                                                     (acoustic_bc_params%focLen + acoustic_bc_params%focCal))**2._wp)
                                rbeta = sqrt(1._wp + ((0.5_wp*acoustic_bc_params%apert)/ &
                                                      (acoustic_bc_params%focLen + acoustic_bc_params%focCal))**2._wp)
                                radial_cc = sqrt(y_cc(k)**2._wp + x_cc(l)**2._wp)

                                if (radial_cc < rc) then
                                    tau = mytime + radial_cc**2._wp/(2._wp*acoustic_bc_params%cson* &
                                                                     (acoustic_bc_params%focLen + acoustic_bc_params%focCal))
                                    gFun = (1._wp/rbeta)

                                    if (tau < (acoustic_bc_params%ncycles/acoustic_bc_params%freq)) then
                                        q_prim_vf(momxe + 1)%sf(k, l, -j) = acoustic_bc_params%Pbase + &
                                                                            acoustic_bc_params%Pamp* &
                                                                            sin(2._wp*pi*acoustic_bc_params%freq*tau)
                                        q_prim_vf(1)%sf(k, l, -j) = acoustic_bc_params%rho
                                    else
                                        q_prim_vf(momxe + 1)%sf(k, l, -j) = acoustic_bc_params%Pbase
                                    end if
                                else
                                    q_prim_vf(momxe + 1)%sf(k, l, -j) = acoustic_bc_params%Pbase
                                end if

                            else
                                stop "acoustic_bc_params%iwave incorrect value (1: planar wave, 2: axisymmetric spherical transducer)."
                            end if

                        end do
                    end do
                end do

            else !< bc_z%end

                call s_mpi_abort('Transucer not available for z dir')

            end if

        end if
        !< =====================================================================

    end subroutine s_acoustic_bc

    subroutine s_ghost_cell_extrapolation(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer :: j, k, l, q, i

        !< x-direction
        if (bc_dir == 1) then !< x-direction

            if (bc_loc == -1) then !bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction
        elseif (bc_dir == 2) then !< y-direction

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, 0, k)
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, n, k)
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction
        elseif (bc_dir == 3) then !< z-direction

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, 0)
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, p)
                            end do
                        end do
                    end do
                end do

            end if

        end if

    end subroutine s_ghost_cell_extrapolation

    subroutine s_symmetry(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            !$acc loop seq
                            do i = 1, contxe
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(j - 1, k, l)
                            end do

                            q_prim_vf(momxb)%sf(-j, k, l) = &
                                -q_prim_vf(momxb)%sf(j - 1, k, l)

                            !$acc loop seq
                            do i = momxb + 1, sys_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(j - 1, k, l)
                            end do

                            if (hyperelasticity) then
                                q_prim_vf(xibeg)%sf(-j, k, l) = &
                                    -q_prim_vf(xibeg)%sf(j - 1, k, l)
                            end if

                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(-j, k, l, q, i) = &
                                            pb(j - 1, k, l, q, i)
                                        mv(-j, k, l, q, i) = &
                                            mv(j - 1, k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif

            else !< bc_x%end

                !$acc parallel loop collapse(3) default(present)
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size

                            !$acc loop seq
                            do i = 1, contxe
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(m - (j - 1), k, l)
                            end do

                            q_prim_vf(momxb)%sf(m + j, k, l) = &
                                -q_prim_vf(momxb)%sf(m - (j - 1), k, l)

                            !$acc loop seq
                            do i = momxb + 1, sys_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(m - (j - 1), k, l)
                            end do

                            if (hyperelasticity) then
                                q_prim_vf(xibeg)%sf(m + j, k, l) = &
                                    -q_prim_vf(xibeg)%sf(m - (j - 1), k, l)
                            end if

                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(m + j, k, l, q, i) = &
                                            pb(m - (j - 1), k, l, q, i)
                                        mv(m + j, k, l, q, i) = &
                                            mv(m - (j - 1), k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            end if

            !< y-direction
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(3) gang vector default(present)
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k)
                            end do

                            q_prim_vf(momxb + 1)%sf(l, -j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, j - 1, k)

                            !$acc loop seq
                            do i = momxb + 2, sys_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k)
                            end do

                            if (hyperelasticity) then
                                q_prim_vf(xibeg + 1)%sf(l, -j, k) = &
                                    -q_prim_vf(xibeg + 1)%sf(l, j - 1, k)
                            end if
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, -j, k, q, i) = &
                                            pb(l, j - 1, k, q, i)
                                        mv(l, -j, k, q, i) = &
                                            mv(l, j - 1, k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            else !< bc_y%end

                !$acc parallel loop collapse(3) gang vector default(present)
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, n - (j - 1), k)
                            end do

                            q_prim_vf(momxb + 1)%sf(l, n + j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, n - (j - 1), k)

                            !$acc loop seq
                            do i = momxb + 2, sys_size
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, n - (j - 1), k)
                            end do

                            if (hyperelasticity) then
                                q_prim_vf(xibeg + 1)%sf(l, n + j, k) = &
                                    -q_prim_vf(xibeg + 1)%sf(l, n - (j - 1), k)
                            end if
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, n + j, k, q, i) = &
                                            pb(l, n - (j - 1), k, q, i)
                                        mv(l, n + j, k, q, i) = &
                                            mv(l, n - (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            end if

            !< z-direction
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            !$acc loop seq
                            do i = 1, momxb + 1
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, j - 1)
                            end do

                            q_prim_vf(momxe)%sf(k, l, -j) = &
                                -q_prim_vf(momxe)%sf(k, l, j - 1)

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, j - 1)
                            end do

                            if (hyperelasticity) then
                                q_prim_vf(xiend)%sf(k, l, -j) = &
                                    -q_prim_vf(xiend)%sf(k, l, j - 1)
                            end if
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, -j, q, i) = &
                                            pb(k, l, j - 1, q, i)
                                        mv(k, l, -j, q, i) = &
                                            mv(k, l, j - 1, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            else !< bc_z%end

                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            !$acc loop seq
                            do i = 1, momxb + 1
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, p - (j - 1))
                            end do

                            q_prim_vf(momxe)%sf(k, l, p + j) = &
                                -q_prim_vf(momxe)%sf(k, l, p - (j - 1))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, p - (j - 1))
                            end do

                            if (hyperelasticity) then
                                q_prim_vf(xiend)%sf(k, l, p + j) = &
                                    -q_prim_vf(xiend)%sf(k, l, p - (j - 1))
                            end if
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, p + j, q, i) = &
                                            pb(k, l, p - (j - 1), q, i)
                                        mv(k, l, p + j, q, i) = &
                                            mv(k, l, p - (j - 1), q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            end if

        end if

    end subroutine s_symmetry

    subroutine s_periodic(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(m - (j - 1), k, l)
                            end do
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(-j, k, l, q, i) = &
                                            pb(m - (j - 1), k, l, q, i)
                                        mv(-j, k, l, q, i) = &
                                            mv(m - (j - 1), k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(j - 1, k, l)
                            end do
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(m + j, k, l, q, i) = &
                                            pb(j - 1, k, l, q, i)
                                        mv(m + j, k, l, q, i) = &
                                            mv(j - 1, k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            end if

            !< y-direction
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, n - (j - 1), k)
                            end do
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, -j, k, q, i) = &
                                            pb(l, n - (j - 1), k, q, i)
                                        mv(l, -j, k, q, i) = &
                                            mv(l, n - (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k)
                            end do
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, n + j, k, q, i) = &
                                            pb(l, (j - 1), k, q, i)
                                        mv(l, n + j, k, q, i) = &
                                            mv(l, (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            end if

            !< z-direction
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, p - (j - 1))
                            end do
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, -j, q, i) = &
                                            pb(k, l, p - (j - 1), q, i)
                                        mv(k, l, -j, q, i) = &
                                            mv(k, l, p - (j - 1), q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, j - 1)
                            end do
                        end do
                    end do
                end do
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, p + j, q, i) = &
                                            pb(k, l, j - 1, q, i)
                                        mv(k, l, p + j, q, i) = &
                                            mv(k, l, j - 1, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#endif
            end if

        end if

    end subroutine s_periodic

    subroutine s_periodic_pair(bc_dir, bc_dir_pair, bc_loc_pair)

        integer, intent(in) :: bc_dir
        integer, intent(out) :: bc_dir_pair, bc_loc_pair

        if (bc_dir == 1) then
            if (bc_y%beg == -22) then; bc_dir_pair = 2; bc_loc_pair = -1; end if
            if (bc_y%end == -22) then; bc_dir_pair = 2; bc_loc_pair = 1; end if
            if (bc_z%beg == -22) then; bc_dir_pair = 3; bc_loc_pair = -1; end if
            if (bc_z%end == -22) then; bc_dir_pair = 3; bc_loc_pair = 1; end if
        elseif (bc_dir == 2) then
            if (bc_x%beg == -22) then; bc_dir_pair = 1; bc_loc_pair = -1; end if
            if (bc_x%end == -22) then; bc_dir_pair = 1; bc_loc_pair = 1; end if
            if (bc_z%beg == -22) then; bc_dir_pair = 3; bc_loc_pair = -1; end if
            if (bc_z%end == -22) then; bc_dir_pair = 3; bc_loc_pair = 1; end if
        else
            if (bc_x%beg == -22) then; bc_dir_pair = 1; bc_loc_pair = -1; end if
            if (bc_x%end == -22) then; bc_dir_pair = 1; bc_loc_pair = 1; end if
            if (bc_y%beg == -22) then; bc_dir_pair = 2; bc_loc_pair = -1; end if
            if (bc_y%end == -22) then; bc_dir_pair = 2; bc_loc_pair = 1; end if
        end if

    end subroutine s_periodic_pair

    ! Rotate the velocity vector usinf Rodrigues' formula
    subroutine s_rotate_velocity(vel_vect, axis_rot, theta, vel_rot)
        !$acc routine seq
        real(wp), dimension(3), intent(in) :: vel_vect, axis_rot
        real(wp), intent(in) :: theta
        real(wp), dimension(3), intent(out) :: vel_rot
        real(wp) :: cos_theta, sin_theta, dot_kv
        real(wp), dimension(3) :: cross_kv

        cos_theta = cos(theta)
        sin_theta = sin(theta)

        if (f_approx_equal(cos_theta, 0._wp)) cos_theta = 0._wp
        if (f_approx_equal(sin_theta, 0._wp)) sin_theta = 0._wp
        if (f_approx_equal(cos_theta, 1._wp)) cos_theta = 1._wp
        if (f_approx_equal(sin_theta, 1._wp)) sin_theta = 1._wp

        ! Compute cross product k × v
        cross_kv(1) = axis_rot(2)*vel_vect(3) - axis_rot(3)*vel_vect(2)
        cross_kv(2) = axis_rot(3)*vel_vect(1) - axis_rot(1)*vel_vect(3)
        cross_kv(3) = axis_rot(1)*vel_vect(2) - axis_rot(2)*vel_vect(1)

        ! Compute dot product k · v
        dot_kv = axis_rot(1)*vel_vect(1) + axis_rot(2)*vel_vect(2) + &
                 axis_rot(3)*vel_vect(3)

        ! Apply Rodrigues' formula
        vel_rot(1) = vel_vect(1)*cos_theta + cross_kv(1)*sin_theta + &
                     axis_rot(1)*dot_kv*(1._wp - cos_theta)
        vel_rot(2) = vel_vect(2)*cos_theta + cross_kv(2)*sin_theta + &
                     axis_rot(2)*dot_kv*(1._wp - cos_theta)
        vel_rot(3) = vel_vect(3)*cos_theta + cross_kv(3)*sin_theta + &
                     axis_rot(3)*dot_kv*(1._wp - cos_theta)

    end subroutine s_rotate_velocity

    function f_rotation_angle(pos, pos_rot)
        !$acc routine seq
        real(wp), dimension(3), intent(in) :: pos, pos_rot
        real(wp) :: dot_pr, norm_pos, norm_rot, cos_theta
        real(wp) :: f_rotation_angle

        ! Compute dot product
        dot_pr = pos(1)*pos_rot(1) + pos(2)*pos_rot(2) + pos(3)*pos_rot(3)

        ! Compute magnitudes
        norm_pos = sqrt(pos(1)**2._wp + pos(2)**2._wp + pos(3)**2._wp)
        norm_rot = sqrt(pos_rot(1)**2._wp + pos_rot(2)**2._wp + pos_rot(3)**2._wp)

        ! Prevent divide by zero
        if (f_approx_equal(norm_pos, 0._wp) .or. f_approx_equal(norm_rot, 0._wp)) then
            f_rotation_angle = 0._wp  ! or return -1 to indicate invalid input
            print *, 'Division by zero in f_rotation_angle: Rotationally periodic BC'
            return
        end if

        ! Compute cosine of angle
        cos_theta = dot_pr/(norm_pos*norm_rot)

        ! Clamp the value to [-1, 1] to avoid domain error in acos
        ! cos_theta = max(-1._wp, min(1._wp, cos_theta))

        ! Compute angle in radians
        f_rotation_angle = acos(cos_theta)

    end function f_rotation_angle

    subroutine s_periodic_rotational(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: bc_dir_pair, bc_loc_pair
        integer :: j, k, l, q, i
        real(wp) :: theta
        real(wp), dimension(3) :: vel_vect, axis_rot, vel_rot, pos, pos_rot

        call s_periodic_pair(bc_dir, bc_dir_pair, bc_loc_pair)

        !< x-direction
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                if (bc_dir_pair == 2) then

                    if (bc_loc_pair == -1) then !< bc_y%beg -> pair

                        !call s_mpi_abort('xbeg <- ybeg need to be implemented (periodic rotational)')
                        !$acc parallel loop collapse(3) gang vector default(present) &
                        !$acc private(pos, pos_rot, vel_vect, axis_rot, vel_rot)
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size

                                    !$acc loop seq
                                    do i = 1, contxe
                                        q_prim_vf(i)%sf(-j, k, l) = &
                                            q_prim_vf(i)%sf(k, j - 1, l)
                                    end do

                                    !> Unit vector of the rotation axis
                                    axis_rot(1) = 0._wp
                                    axis_rot(2) = 0._wp
                                    axis_rot(3) = 1._wp

                                    !> Rotation angle
                                    pos(1:3) = 0._wp; pos_rot(1:3) = 0._wp
                                    pos(1) = x_cc(-j); pos_rot(1) = x_cc(k)
                                    pos(2) = y_cc(k); pos_rot(2) = y_cc(j - 1)
                                    if (p /= 0) pos(3) = z_cc(l)
                                    if (p /= 0) pos_rot(3) = z_cc(l)
                                    theta = f_rotation_angle(pos, pos_rot)

                                    vel_vect(1:3) = 0._wp
                                    vel_vect(1) = q_prim_vf(momxb)%sf(k, j - 1, l)
                                    vel_vect(2) = q_prim_vf(momxb + 1)%sf(k, j - 1, l)
                                    if (p /= 0) vel_vect(3) = q_prim_vf(momxb + 2)%sf(k, j - 1, l)

                                    call s_rotate_velocity(vel_vect, axis_rot, theta, vel_rot)

                                    !> Rotate velocty
                                    !$acc loop seq
                                    do i = momxb, momxe
                                        q_prim_vf(i)%sf(-j, k, l) = vel_rot(i - momxb + 1)
                                    end do

                                    !$acc loop seq
                                    do i = E_idx, sys_size
                                        q_prim_vf(i)%sf(-j, k, l) = &
                                            q_prim_vf(i)%sf(k, j - 1, l)
                                    end do

                                end do
                            end do
                        end do

                    else !< bc_y%end -> pair

                        call s_mpi_abort('xbeg <- yend need to be implemented (periodic rotational)')

                    end if

                else

                    if (bc_loc_pair == -1) then !< bc_z%beg -> pair

                        call s_mpi_abort('xbeg <- zbeg need to be implemented (periodic rotational)')

                    else !< bc_z%end -> pair

                        call s_mpi_abort('xbeg <- zend need to be implemented (periodic rotational)')

                    end if

                end if

            else !< bc_x%end

                if (bc_dir_pair == 2) then

                    if (bc_loc_pair == -1) then !< bc_y%beg -> pair

                        call s_mpi_abort('xend <- ybeg need to be implemented (periodic rotational)')

                    else !< bc_y%end -> pair

                        call s_mpi_abort('xend <- yend need to be implemented (periodic rotational)')

                    end if

                else

                    if (bc_loc_pair == -1) then !< bc_z%beg -> pair

                        call s_mpi_abort('xend <- zbeg need to be implemented (periodic rotational)')

                    else !< bc_z%end -> pair

                        call s_mpi_abort('xend <- zend need to be implemented (periodic rotational)')

                    end if

                end if

            end if

            ! y-direction
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                if (bc_dir_pair == 1) then

                    if (bc_loc_pair == -1) then !< bc_x%beg -> pair

                        ! call s_mpi_abort('ybeg <- xbeg need to be implemented (periodic rotational)')
                        !$acc parallel loop collapse(3) gang vector default(present) &
                        !$acc private(pos, pos_rot, vel_vect, axis_rot, vel_rot)
                        do l = 0, p
                            do j = 1, buff_size
                                do k = 0, m

                                    !$acc loop seq
                                    do i = 1, contxe
                                        q_prim_vf(i)%sf(k, -j, l) = &
                                            q_prim_vf(i)%sf(j - 1, k, l)
                                    end do

                                    !> Unit vector of the rotation axis
                                    axis_rot(1) = 0._wp
                                    axis_rot(2) = 0._wp
                                    axis_rot(3) = -1._wp

                                    !> Rotation angle
                                    pos(1:3) = 0._wp; pos_rot(1:3) = 0._wp
                                    pos(1) = x_cc(k); pos_rot(1) = x_cc(j - 1)
                                    pos(2) = y_cc(-j); pos_rot(2) = y_cc(k)
                                    if (p /= 0) pos(3) = z_cc(l)
                                    if (p /= 0) pos_rot(3) = z_cc(l)
                                    theta = abs(f_rotation_angle(pos, pos_rot))

                                    vel_vect(1:3) = 0._wp
                                    vel_vect(1) = q_prim_vf(momxb)%sf(j - 1, k, l)
                                    vel_vect(2) = q_prim_vf(momxb + 1)%sf(j - 1, k, l)
                                    if (p /= 0) vel_vect(3) = q_prim_vf(momxb + 2)%sf(j - 1, k, l)

                                    call s_rotate_velocity(vel_vect, axis_rot, theta, vel_rot)

                                    !> Rotate velocty
                                    !$acc loop seq
                                    do i = momxb, momxe
                                        q_prim_vf(i)%sf(k, -j, l) = vel_rot(i - momxb + 1)
                                    end do

                                    !$acc loop seq
                                    do i = E_idx, sys_size
                                        q_prim_vf(i)%sf(k, -j, l) = &
                                            q_prim_vf(i)%sf(j - 1, k, l)
                                    end do

                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(3) gang vector default(present) &
                        !$acc private(pos, pos_rot, vel_vect, axis_rot, vel_rot)
                        do l = 0, p
                            do j = 1, buff_size
                                do k = -buff_size, -1

                                    !$acc loop seq
                                    do i = 1, contxe
                                        q_prim_vf(i)%sf(k, -j, l) = &
                                            q_prim_vf(i)%sf(j - 1, k, l)
                                    end do

                                    !> Unit vector of the rotation axis
                                    axis_rot(1) = 0._wp
                                    axis_rot(2) = 0._wp
                                    axis_rot(3) = -1._wp

                                    !> Rotation angle
                                    pos(1:3) = 0._wp; pos_rot(1:3) = 0._wp
                                    pos(1) = x_cc(k); pos_rot(1) = x_cc(j - 1)
                                    pos(2) = y_cc(-j); pos_rot(2) = y_cc(k)
                                    if (p /= 0) pos(3) = z_cc(l)
                                    if (p /= 0) pos_rot(3) = z_cc(l)
                                    theta = abs(f_rotation_angle(pos, pos_rot))

                                    vel_vect(1:3) = 0._wp
                                    vel_vect(1) = q_prim_vf(momxb)%sf(j - 1, k, l)
                                    vel_vect(2) = q_prim_vf(momxb + 1)%sf(j - 1, k, l)
                                    if (p /= 0) vel_vect(3) = q_prim_vf(momxb + 2)%sf(j - 1, k, l)

                                    call s_rotate_velocity(vel_vect, axis_rot, theta, vel_rot)

                                    !> Rotate velocty
                                    !$acc loop seq
                                    do i = momxb, momxe
                                        q_prim_vf(i)%sf(k, -j, l) = vel_rot(i - momxb + 1)
                                    end do

                                    !$acc loop seq
                                    do i = E_idx, sys_size
                                        q_prim_vf(i)%sf(k, -j, l) = &
                                            q_prim_vf(i)%sf(j - 1, k, l)
                                    end do

                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(4) gang vector default(present)
                        do i = 1, sys_size
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = m + 1, m + buff_size
                                        q_prim_vf(i)%sf(l, -j, k) = &
                                            q_prim_vf(i)%sf(m, -j, k)
                                    end do
                                end do
                            end do
                        end do

                    else !< bc_x%end -> pair

                        call s_mpi_abort('ybeg <- xend need to be implemented (periodic rotational)')

                    end if

                else

                    if (bc_loc_pair == -1) then !< bc_z%beg -> pair

                        call s_mpi_abort('ybeg <- zbeg need to be implemented (periodic rotational)')
                        ! !$acc parallel loop collapse(3) gang vector default(present)
                        ! do k = 0, p
                        !     do j = 1, buff_size
                        !         do l = -buff_size, m + buff_size

                        !             !$acc loop seq
                        !             do i = 1, contxe
                        !                 q_prim_vf(i)%sf(l, -j, k) = &
                        !                         q_prim_vf(i)%sf(l, k, j - 1)
                        !             end do

                        !             !> Unit vector of the rotation axis
                        !             axis_rot(1) = 1._wp
                        !             axis_rot(2) = 0._wp
                        !             axis_rot(3) = 0._wp

                        !             !> Rotation angle
                        !             pos(1:3) = 0._wp; pos_rot(1:3) = 0._wp
                        !             pos(1) = x_cc(l); pos_rot(1) = x_cc(l)
                        !             pos(2) = y_cc(-j); pos_rot(2) = y_cc(k)
                        !             if (p /= 0 ) pos(3) = z_cc(k)
                        !             if (p /= 0 ) pos_rot(3) = z_cc(j - 1)
                        !             theta = abs(f_rotation_angle(pos, pos_rot))

                        !             vel_vect(1:3) = 0._wp
                        !             vel_vect(1) = q_prim_vf(1)%sf(l, k, j - 1)
                        !             vel_vect(2) = q_prim_vf(2)%sf(l, k, j - 1)
                        !             if (p /= 0 ) vel_vect(3) = q_prim_vf(3)%sf(l, k, j - 1)

                        !             call s_rotate_velocity(vel_vect, axis_rot, theta, vel_rot)

                        !             !> Rotate velocty
                        !             !$acc loop seq
                        !             do i = momxb, momxe
                        !                 q_prim_vf(i)%sf(l, -j, k) = vel_rot(i - momxb + 1)
                        !             end do

                        !             !$acc loop seq
                        !             do i = E_idx, sys_size
                        !                 q_prim_vf(i)%sf(l, -j, k) = &
                        !                         q_prim_vf(i)%sf(l, k, j - 1)
                        !             end do

                        !         end do
                        !     end do
                        ! end do

                    else !< bc_z%end -> pair

                        call s_mpi_abort('ybeg <- zend need to be implemented (periodic rotational)')
                        ! !$acc parallel loop collapse(3) gang vector default(present)
                        ! do k = 0, p
                        !     do j = 1, buff_size
                        !         do l = -buff_size, m + buff_size

                        !             !$acc loop seq
                        !             do i = 1, contxe
                        !                 q_prim_vf(i)%sf(l, -j, k) = &
                        !                         q_prim_vf(i)%sf(l, k, p - (j - 1))
                        !             end do

                        !             !> Unit vector of the rotation axis and rotation angle
                        !             a1 = -1._wp
                        !             a2 = 0._wp
                        !             a3 = 0._wp
                        !             theta = pi/4._wp
                        !             vel_vect(1) = q_prim_vf(1)%sf(l, k, p - (j - 1))
                        !             vel_vect(2) = q_prim_vf(2)%sf(l, k, p - (j - 1))
                        !             if (p /= 0 ) vel_vect(3) = q_prim_vf(3)%sf(l, k, p - (j - 1))

                        !             !> Rotate velocty
                        !             !$acc loop seq
                        !             do i = momxb, momxe
                        !                 q_prim_vf(i)%sf(l, -j, k) = f_rotate_velocity(vel_vect, a1, a2, a3, theta, i)
                        !             end do

                        !             !$acc loop seq
                        !             do i = E_idx, sys_size
                        !                 q_prim_vf(i)%sf(l, -j, k) = &
                        !                         q_prim_vf(i)%sf(l, k, p - (j - 1))
                        !             end do

                        !         end do
                        !     end do
                        ! end do

                    end if

                end if

            else !< bc_y%end

                if (bc_dir_pair == 1) then

                    if (bc_loc_pair == -1) then !< bc_x%beg -> pair

                        call s_mpi_abort('yend <- xbeg need to be implemented (periodic rotational)')

                    else !< bc_x%end -> pair

                        call s_mpi_abort('yend <- xend need to be implemented (periodic rotational)')

                    end if

                else

                    if (bc_loc_pair == -1) then !< bc_z%beg -> pair

                        call s_mpi_abort('yend <- zbeg need to be implemented (periodic rotational)')

                    else !< bc_z%end -> pair

                        call s_mpi_abort('yend <- zend need to be implemented (periodic rotational)')

                    end if

                end if

            end if

            ! z-direction
        else

            if (bc_loc == -1) then !< bc_z%beg

                if (bc_dir_pair == 1) then

                    if (bc_loc_pair == -1) then !< bc_x%beg -> pair

                        call s_mpi_abort('zbeg <- xbeg need to be implemented (periodic rotational)')

                    else !< bc_x%end -> pair

                        call s_mpi_abort('zbeg <- xend need to be implemented (periodic rotational)')

                    end if

                else

                    if (bc_loc_pair == -1) then !< bc_y%beg -> pair

                        call s_mpi_abort('zbeg <- ybeg need to be implemented (periodic rotational)') !!!NEEEDEEED!!
                        ! !$acc parallel loop collapse(3) gang vector default(present)
                        ! do j = 1, buff_size
                        !     do l = -buff_size, n + buff_size
                        !         do k = -buff_size, m + buff_size
                        !             !$acc loop seq
                        !             do i = 1, contxe
                        !                 q_prim_vf(i)%sf(k, l, -j) = &
                        !                         q_prim_vf(i)%sf(k, j - 1, l)
                        !             end do

                        !             !> Unit vector of the rotation axis
                        !             axis_rot(1) = -1._wp
                        !             axis_rot(2) = 0._wp
                        !             axis_rot(3) = 0._wp

                        !             !> Rotation angle
                        !             pos(1:3) = 0._wp; pos_rot(1:3) = 0._wp
                        !             pos(1) = x_cc(k); pos_rot(1) = x_cc(k)
                        !             pos(2) = y_cc(l); pos_rot(2) = y_cc(j - 1)
                        !             if (p /= 0 ) pos(3) = z_cc(-j)
                        !             if (p /= 0 ) pos_rot(3) = z_cc(l)
                        !             theta = abs(f_rotation_angle(pos, pos_rot))

                        !             vel_vect(1:3) = 0._wp
                        !             vel_vect(1) = q_prim_vf(1)%sf(k, j - 1, l)
                        !             vel_vect(2) = q_prim_vf(2)%sf(k, j - 1, l)
                        !             if (p /= 0 ) vel_vect(3) = q_prim_vf(3)%sf(k, j - 1, l)

                        !             call s_rotate_velocity(vel_vect, axis_rot, theta, vel_rot)

                        !             !> Rotate velocty
                        !             !$acc loop seq
                        !             do i = momxb, momxe
                        !                 q_prim_vf(i)%sf(k, l, -j) = vel_rot(i - momxb + 1)
                        !             end do

                        !             !$acc loop seq
                        !             do i = E_idx, sys_size
                        !                 q_prim_vf(i)%sf(k, l, -j) = &
                        !                         q_prim_vf(i)%sf(k, j - 1, l)
                        !             end do

                        !         end do
                        !     end do
                        ! end do

                    else !< bc_z%end -> pair

                        call s_mpi_abort('zbeg <- yend need to be implemented (periodic rotational)')

                    end if

                end if

            else !< bc_z%end

                if (bc_dir_pair == 1) then

                    if (bc_loc_pair == -1) then !< bc_x%beg -> pair

                        call s_mpi_abort('zend <- xbeg need to be implemented (periodic rotational)')

                    else !< bc_x%end -> pair

                        call s_mpi_abort('zend <- xend need to be implemented (periodic rotational)')

                    end if

                else

                    if (bc_loc_pair == -1) then !< bc_y%beg -> pair

                        call s_mpi_abort('zend <- ybeg need to be implemented (periodic rotational)')

                    else !< bc_z%end -> pair

                        call s_mpi_abort('zend <- yend need to be implemented (periodic rotational)')

                    end if

                end if

            end if

        end if

    end subroutine s_periodic_rotational

    subroutine s_axis(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 0, p
            do j = 1, buff_size
                do l = -buff_size, m + buff_size

                    if (z_cc(k) < pi) then

                        if (hifu_params%heatSolver) then
                            !Temperature only
                            q_prim_vf(1)%sf(l, -j, k) = &
                                q_prim_vf(1)%sf(l, j - 1, k + ((p + 1)/2))
                        else

                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k + ((p + 1)/2))
                            end do

                            q_prim_vf(momxb + 1)%sf(l, -j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, j - 1, k + ((p + 1)/2))

                            q_prim_vf(momxe)%sf(l, -j, k) = &
                                -q_prim_vf(momxe)%sf(l, j - 1, k + ((p + 1)/2))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k + ((p + 1)/2))
                            end do

                        end if
                    else

                        if (hifu_params%heatSolver) then
                            !Temperature only
                            q_prim_vf(1)%sf(l, -j, k) = &
                                q_prim_vf(1)%sf(l, j - 1, k - ((p + 1)/2))
                        else

                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k - ((p + 1)/2))
                            end do

                            q_prim_vf(momxb + 1)%sf(l, -j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, j - 1, k - ((p + 1)/2))

                            q_prim_vf(momxe)%sf(l, -j, k) = &
                                -q_prim_vf(momxe)%sf(l, j - 1, k - ((p + 1)/2))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k - ((p + 1)/2))
                            end do

                        end if
                    end if
                end do
            end do
        end do
#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do q = 1, nnode
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                pb(l, -j, k, q, i) = &
                                    pb(l, j - 1, k - ((p + 1)/2), q, i)
                                mv(l, -j, k, q, i) = &
                                    mv(l, j - 1, k - ((p + 1)/2), q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if
#endif
    end subroutine s_axis

    subroutine s_axis_cylindrical_sector_hifu(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 0, p
            do j = 1, buff_size
                do l = -buff_size, m + buff_size
                    q_prim_vf(hifu_params%T_idx)%sf(l, -j, k) = &
                        q_prim_vf(hifu_params%T_idx)%sf(l, j - 1, k)
                end do
            end do
        end do

    end subroutine s_axis_cylindrical_sector_hifu

    subroutine s_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                                else
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        q_prim_vf(i)%sf(0, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve1
                                else
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        q_prim_vf(i)%sf(m, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb + 1) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb2
                                else
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        q_prim_vf(i)%sf(l, 0, k)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb + 1) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve2
                                else
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        q_prim_vf(i)%sf(l, n, k)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxe) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb3
                                else
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        q_prim_vf(i)%sf(k, l, 0)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxe) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve3
                                else
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        q_prim_vf(i)%sf(k, l, p)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

        end if

    end subroutine s_slip_wall

    subroutine s_no_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb3
                                else
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        q_prim_vf(i)%sf(0, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve3
                                else
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        q_prim_vf(i)%sf(m, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb3
                                else
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        q_prim_vf(i)%sf(l, 0, k)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve3
                                else
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        q_prim_vf(i)%sf(l, n, k)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb3
                                else
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        q_prim_vf(i)%sf(k, l, 0)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve3
                                else
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        q_prim_vf(i)%sf(k, l, p)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

        end if

    end subroutine s_no_slip_wall

    subroutine s_qbmm_extrapolation(pb, mv, bc_dir, bc_loc)

        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    pb(-j, k, l, q, i) = &
                                        pb(0, k, l, q, i)
                                    mv(-j, k, l, q, i) = &
                                        mv(0, k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    pb(m + j, k, l, q, i) = &
                                        pb(m, k, l, q, i)
                                    mv(m + j, k, l, q, i) = &
                                        mv(m, k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    pb(l, -j, k, q, i) = &
                                        pb(l, 0, k, q, i)
                                    mv(l, -j, k, q, i) = &
                                        mv(l, 0, k, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    pb(l, n + j, k, q, i) = &
                                        pb(l, n, k, q, i)
                                    mv(l, n + j, k, q, i) = &
                                        mv(l, n, k, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    pb(k, l, -j, q, i) = &
                                        pb(k, l, 0, q, i)
                                    mv(k, l, -j, q, i) = &
                                        mv(k, l, 0, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    pb(k, l, p + j, q, i) = &
                                        pb(k, l, p, q, i)
                                    mv(k, l, p + j, q, i) = &
                                        mv(k, l, p, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            end if

        end if

    end subroutine s_qbmm_extrapolation

#ifdef MFC_SIMULATION
    subroutine s_populate_capillary_buffers(c_divs)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer :: i, j, k, l

        ! x - direction
        if (bc_x%beg <= -3) then !< ghost cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs(i)%sf(-j, k, l) = &
                                c_divs(i)%sf(0, k, l)
                        end do
                    end do
                end do
            end do
        elseif (bc_x%beg == -2) then !< slip wall or reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            if (i == 1) then
                                c_divs(i)%sf(-j, k, l) = &
                                    -c_divs(i)%sf(j - 1, k, l)
                            else
                                c_divs(i)%sf(-j, k, l) = &
                                    c_divs(i)%sf(j - 1, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        elseif (bc_x%beg == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs(i)%sf(-j, k, l) = &
                                c_divs(i)%sf(m - (j - 1), k, l)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, -1)
        end if

        if (bc_x%end <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs(i)%sf(m + j, k, l) = &
                                c_divs(i)%sf(m, k, l)
                        end do
                    end do
                end do
            end do
        elseif (bc_x%end == -2) then
            !$acc parallel loop collapse(4) default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            if (i == 1) then
                                c_divs(i)%sf(m + j, k, l) = &
                                    -c_divs(i)%sf(m - (j - 1), k, l)
                            else
                                c_divs(i)%sf(m + j, k, l) = &
                                    c_divs(i)%sf(m - (j - 1), k, l)
                            end if
                        end do
                    end do
                end do
            end do
        else if (bc_x%end == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs(i)%sf(m + j, k, l) = &
                                c_divs(i)%sf(j - 1, k, l)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, 1)
        end if

        if (n == 0) then
            return
        elseif (bc_y%beg <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs(i)%sf(l, -j, k) = &
                                c_divs(i)%sf(l, 0, k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%beg == -2) then !< slip wall or reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            if (i == 2) then
                                c_divs(i)%sf(l, -j, k) = &
                                    -c_divs(i)%sf(l, j - 1, k)
                            else
                                c_divs(i)%sf(l, -j, k) = &
                                    c_divs(i)%sf(l, j - 1, k)
                            end if
                        end do
                    end do
                end do
            end do
        elseif (bc_y%beg == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs(i)%sf(l, -j, k) = &
                                c_divs(i)%sf(l, n - (j - 1), k)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 2, -1)
        end if

        if (bc_y%end <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs(i)%sf(l, n + j, k) = &
                                c_divs(i)%sf(l, n, k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%end == -2) then !< slip wall or reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            if (i == 2) then
                                c_divs(i)%sf(l, n + j, k) = &
                                    -c_divs(i)%sf(l, n - (j - 1), k)
                            else
                                c_divs(i)%sf(l, n + j, k) = &
                                    c_divs(i)%sf(l, n - (j - 1), k)
                            end if
                        end do
                    end do
                end do
            end do
        elseif (bc_y%end == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs(i)%sf(l, n + j, k) = &
                                c_divs(i)%sf(l, j - 1, k)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 2, 1)
        end if

        if (p == 0) then
            return
        elseif (bc_z%beg <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs(i)%sf(k, l, -j) = &
                                c_divs(i)%sf(k, l, 0)
                        end do
                    end do
                end do
            end do
        elseif (bc_z%beg == -2) then !< symmetry
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            if (i == 3) then
                                c_divs(i)%sf(k, l, -j) = &
                                    -c_divs(i)%sf(k, l, j - 1)
                            else
                                c_divs(i)%sf(k, l, -j) = &
                                    c_divs(i)%sf(k, l, j - 1)
                            end if
                        end do
                    end do
                end do
            end do
        elseif (bc_z%beg == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs(i)%sf(k, l, -j) = &
                                c_divs(i)%sf(k, l, p - (j - 1))
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 3, -1)
        end if

        if (bc_z%end <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs(i)%sf(k, l, p + j) = &
                                c_divs(i)%sf(k, l, p)
                        end do
                    end do
                end do
            end do
        elseif (bc_z%end == -2) then !< symmetry
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            if (i == 3) then
                                c_divs(i)%sf(k, l, p + j) = &
                                    -c_divs(i)%sf(k, l, p - (j - 1))
                            else
                                c_divs(i)%sf(k, l, p + j) = &
                                    c_divs(i)%sf(k, l, p - (j - 1))
                            end if
                        end do
                    end do
                end do
            end do
        elseif (bc_z%end == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs(i)%sf(k, l, p + j) = &
                                c_divs(i)%sf(k, l, j - 1)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 3, 1)
        end if

    end subroutine s_populate_capillary_buffers
#endif

end module m_boundary_conditions
