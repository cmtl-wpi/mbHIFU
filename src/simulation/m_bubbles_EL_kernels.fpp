!>
!! @file m_bubbles_EL_kernels.f90
!! @brief Contains module m_bubbles_EL_kernels

#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian bubbles
!!        in the Eulerian framework.
module m_bubbles_EL_kernels

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    implicit none

    integer :: bcxb, bcxe, bcyb, bcye, bczb, bcze
    !$acc declare create(bcxb, bcxe, bcyb, bcye, bczb, bcze)

contains

    subroutine s_initialize_bubbles_EL_kernels()

        bcxb = bc_x%beg
        bcxe = bc_x%end
        !$acc update device(bcxb, bcxe)

        if (n > 0) then
            bcyb = bc_y%beg
            bcye = bc_y%end
            !$acc update device(bcyb, bcye)
        end if

        if (p > 0) then
            bczb = bc_z%beg
            bcze = bc_z%end
            !$acc update device(bczb, bcze)
        end if

    end subroutine s_initialize_bubbles_EL_kernels

    !> The purpose of this subroutine is to smear the strength of the lagrangian
            !!      bubbles into the Eulerian framework using different approaches.
            !! @param nBubs Number of lagrangian bubbles in the current domain
            !! @param lbk_rad Radius of the bubbles
            !! @param lbk_vel Interface velocity of the bubbles
            !! @param lbk_s Computational coordinates of the bubbles
            !! @param lbk_pos Spatial coordinates of the bubbles
            !! @param updatedvar Eulerian variable to be updated
    subroutine s_smoothfunction(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(vector_field), intent(inout) :: updatedvar
        real(wp), dimension(1:lag_params%nBubs_glb), intent(in), optional :: lbk_qvis, lbk_qth

        if (hifu_params%heatSolver) then
            if (hifu_params%cartesian .or. (p>0 .and. .not. cyl_coord)) then
                !call s_deltafunc_hifu(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)
                call s_gaussian_hifu(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)
            else
                call s_deltafunc(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)
            end if
        else
            smoothfunc:select case(lag_params%smooth_type)
            case (1)
            call s_gaussian(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)
            case (2)
            call s_deltafunc(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)
            end select smoothfunc
        end if

    end subroutine s_smoothfunction

    !> The purpose of this procedure contains the algorithm to use the delta kernel function to map the effect of the bubbles.
            !!      The effect of the bubbles only affects the cell where the bubble is located.
    subroutine s_deltafunc(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(vector_field), intent(inout) :: updatedvar
        real(wp), dimension(1:lag_params%nBubs_glb), intent(in), optional :: lbk_qvis, lbk_qth

        integer, dimension(3) :: cell
        real(wp) :: strength_vel, strength_vol

        real(wp) :: addFun1, addFun2, addFun3
        real(wp) :: volpart, Vol
        real(wp), dimension(3) :: s_coord
        integer :: l

        integer :: thetaCell
        real(wp) :: thetaPos
        logical :: bubble_in_hifu_domain

        !$acc parallel loop gang vector default(present) private(l, s_coord, cell)
        do l = 1, nBubs

            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)

            if (hifu_params%heatSolver .and. num_dims == 3) then
                ! Find azimuthal cell location
                ! -pi to pi -> 0 to 2*pi
                thetaPos = lbk_pos(l, 3, 1)
                if (thetaPos < 0._wp) thetaPos = 2._wp*pi + lbk_pos(l, 3, 1)

                bubble_in_hifu_domain = .true.
                thetaCell = -buff_size
                do while (thetaPos < z_cb(thetaCell - 1) .and. bubble_in_hifu_domain)
                    thetaCell = thetaCell - 1
                    if (thetaCell < -1 - buff_size) then
                        bubble_in_hifu_domain = .false.
                        thetaCell = thetaCell + 1
                    end if
                end do
                do while (thetaPos > z_cb(thetaCell) .and. bubble_in_hifu_domain)
                    thetaCell = thetaCell + 1
                    if (thetaCell > p + buff_size) then
                        bubble_in_hifu_domain = .false.
                        thetaCell = thetaCell - 1
                    end if
                end do
                if (bubble_in_hifu_domain) s_coord(3) = thetaCell + (thetaPos - z_cb(thetaCell - 1))/dz(thetaCell)
            end if

            if (.not. bubble_in_hifu_domain) cycle

            call s_get_cell(s_coord, cell)

            strength_vol = volpart
            strength_vel = 4._wp*pi*lbk_rad(l, 2)**2._wp*lbk_vel(l, 2)

            if (num_dims == 2) then
                Vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                if (cyl_coord) Vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                Vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            end if

            if (hifu_params%heatSolver) then
                ! Smear the viscous and thermal intensities in the 3D domain
                !Update heat source field (qvis)
                addFun1 = lbk_qvis(l)/Vol
                !$acc atomic update
                updatedvar%vf(hifu_params%qvis_idx)%sf(cell(1), cell(2), cell(3)) = &
                    updatedvar%vf(hifu_params%qvis_idx)%sf(cell(1), cell(2), cell(3)) + &
                    addFun1

                !Update heat source field (qth)
                addFun2 = lbk_qth(l)/Vol
                !$acc atomic update
                updatedvar%vf(hifu_params%qth_idx)%sf(cell(1), cell(2), cell(3)) = &
                    updatedvar%vf(hifu_params%qth_idx)%sf(cell(1), cell(2), cell(3)) + &
                    addFun2
            else
                !Update void fraction field
                addFun1 = strength_vol/Vol
                !$acc atomic update
                updatedvar%vf(1)%sf(cell(1), cell(2), cell(3)) = updatedvar%vf(1)%sf(cell(1), cell(2), cell(3)) + addFun1

                !Update time derivative of void fraction
                addFun2 = strength_vel/Vol
                !$acc atomic update
                updatedvar%vf(2)%sf(cell(1), cell(2), cell(3)) = updatedvar%vf(2)%sf(cell(1), cell(2), cell(3)) + addFun2

                !Product of two smeared functions
                !Update void fraction * time derivative of void fraction
                if (p == 0) then
                    addFun3 = (strength_vol*strength_vel)/Vol
                    !$acc atomic update
                    updatedvar%vf(5)%sf(cell(1), cell(2), cell(3)) = updatedvar%vf(5)%sf(cell(1), cell(2), cell(3)) + addFun3
                end if
            end if
        end do

    end subroutine s_deltafunc

    !> The purpose of this procedure contains the algorithm to use the gaussian kernel function to map the effect of the bubbles.
            !!      The effect of the bubbles affects the 3X3x3 cells that surround the bubble.
    subroutine s_gaussian(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(vector_field), intent(inout) :: updatedvar
        real(wp), dimension(1:lag_params%nBubs_glb), intent(in), optional :: lbk_qvis, lbk_qth

        real(wp), dimension(3) :: center
        integer, dimension(3) :: cell
        real(wp) :: stddsv
        real(wp) :: strength_vel, strength_vol, gpu_sum

        real(wp), dimension(3) :: nodecoord
        real(wp) :: addFun1, addFun2, addFun3
        real(wp) :: func, func2, volpart
        integer, dimension(3) :: cellaux
        real(wp), dimension(3) :: s_coord
        integer :: l, i, j, k
        logical :: celloutside
        integer :: smearGrid, smearGridz

        smearGrid = mapCells - (-mapCells) ! Include the cell that contains the bubble (3+1+3)
        smearGridz = smearGrid
        if (p == 0) smearGridz = 0

        gpu_sum = 0._wp
        !$acc parallel loop collapse(4) gang vector default(present) private(l, s_coord, cell, center, cellaux, nodecoord) &
        !$acc copyin(smearGrid, smearGridz) reduction(+:gpu_sum) copy(gpu_sum)
        do l = 1, nBubs
            do i = 0, smearGrid
                do j = 0, smearGrid
                    do k = 0, smearGridz

                        nodecoord(1:3) = 0
                        center(1:3) = 0._wp
                        volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
                        s_coord(1:3) = lbk_s(l, 1:3, 2)
                        center(1:2) = lbk_pos(l, 1:2, 2)
                        if (p > 0) center(3) = lbk_pos(l, 3, 2)
                        call s_get_cell(s_coord, cell)
                        call s_compute_stddsv(cell, volpart, stddsv)
                        strength_vol = volpart
                        strength_vel = 4._wp*pi*lbk_rad(l, 2)**2._wp*lbk_vel(l, 2)

                        cellaux(1) = cell(1) + i - mapCells
                        cellaux(2) = cell(2) + j - mapCells
                        cellaux(3) = cell(3) + k - mapCells
                        if (p == 0) cellaux(3) = 0

                        ! if (i==0 .and. j==0 .and. k==0) print*, cell(1), cell(2), cell(3), l, proc_rank, m, n, p

                        !Check if the cells intended to smear the bubbles in are in the computational domain
                        !and redefine the cells for symmetric boundary
                        call s_check_celloutside(center, cellaux, nodecoord, celloutside)

                        if (.not. celloutside) then

                            nodecoord(1) = x_cc(cellaux(1))
                            nodecoord(2) = y_cc(cellaux(2))
                            if (p > 0) nodecoord(3) = z_cc(cellaux(3))
                            call s_applygaussian(center, cellaux, nodecoord, stddsv, 0._wp, func)
                            if (p == 0 .and. .not. lag_params%newModel_2D) then
                                call s_applygaussian(center, cellaux, nodecoord, stddsv, 1._wp, func2)
                            end if

                            ! Relocate cells for bubbles intersecting symmetric boundaries
                            if (any((/bcxb, bcxe, bcyb, bcye, bczb, bcze/) == BC_REFLECTIVE)) then
                                call s_shift_cell_symmetric_bc(cellaux, cell)
                            end if
                        else
                            func = 0._wp
                            func2 = 0._wp
                            cellaux(1) = cell(1)
                            cellaux(2) = cell(2)
                            cellaux(3) = cell(3)
                            if (p == 0) cellaux(3) = 0
                        end if

                        !Update void fraction field
                        addFun1 = func*strength_vol
                        !$acc atomic update
                        updatedvar%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                            updatedvar%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                            + addFun1

                        !Update time derivative of void fraction
                        addFun2 = func*strength_vel
                        !$acc atomic update
                        updatedvar%vf(2)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                            updatedvar%vf(2)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                            + addFun2

                        !Product of two smeared functions
                        !Update void fraction * time derivative of void fraction
                        if (p == 0 .and. .not. lag_params%newModel_2D) then
                            addFun3 = func2*strength_vol*strength_vel
                            !$acc atomic update
                            updatedvar%vf(5)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar%vf(5)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + addFun3
                        end if

                        if (cellaux(1) >= 0 .and. cellaux(1) <= m .and. &
                            cellaux(2) >= 0 .and. cellaux(2) <= n) then
                            addFun1 = dx(cellaux(1))*dy(cellaux(2))*y_cc(cellaux(2))*2._wp*pi
                        else
                            addFun1 = 0._wp
                        end if

                        gpu_sum = gpu_sum + func*strength_vol*addFun1

                        ! if (i==3 .and. j==3) then ! shows error in the standard deviation
                        !     print*, l, addFun1, strength_vol, stddsv, func, strength_vel
                        ! end if

                    end do
                end do
            end do
        end do

        ! Populate symmetric boundaries
        if (any((/bcxb, bcxe, bcyb, bcye, bczb, bcze/) == BC_REFLECTIVE)) then
            call s_populate_symmetric_bc(updatedvar)
        end if

    end subroutine s_gaussian

    !> The purpose of this subroutine is to apply the gaussian kernel function for each bubble (Maeda and Colonius, 2018)).
    subroutine s_applygaussian(center, cellaux, nodecoord, stddsv, strength_idx, func)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_applygaussian
#else
        !$acc routine seq
#endif
        real(wp), dimension(3), intent(in) :: center
        integer, dimension(3), intent(in) :: cellaux
        real(wp), dimension(3), intent(in) :: nodecoord
        real(wp), intent(in) :: stddsv
        real(wp), intent(in) :: strength_idx
        real(wp), intent(out) :: func

        real(wp) :: distance
        real(wp) :: theta, dtheta, L2, dzp, Lz2
        real(wp) :: Nr, Nr_count

        distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + (center(3) - nodecoord(3))**2._wp)

        if (num_dims == 3 .or. lag_params%newModel_2D) then
            !< 3D gaussian function
            func = exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
        else
            if (cyl_coord) then
                !< 2D cylindrical function:
                ! We smear particles in the azimuthal direction for given r
                theta = 0._wp
                Nr = ceiling(2._wp*pi*nodecoord(2)/(y_cb(cellaux(2)) - y_cb(cellaux(2) - 1)))
                dtheta = 2._wp*pi/Nr
                L2 = center(2)**2._wp + nodecoord(2)**2._wp - 2._wp*center(2)*nodecoord(2)*cos(theta)
                distance = sqrt((center(1) - nodecoord(1))**2._wp + L2)
                ! Factor 2._wp is for symmetry (upper half of the 2D field (+r) is considered)
                func = dtheta/2._wp/pi*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
                Nr_count = 0._wp
                do while (Nr_count < Nr - 1._wp)
                    Nr_count = Nr_count + 1._wp
                    theta = Nr_count*dtheta
                    ! trigonometric relation
                    L2 = center(2)**2._wp + nodecoord(2)**2._wp - 2._wp*center(2)*nodecoord(2)*cos(theta)
                    distance = sqrt((center(1) - nodecoord(1))**2._wp + L2)
                    ! nodecoord(2)*dtheta is the azimuthal width of the cell
                    func = func + &
                           dtheta/2._wp/pi*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**(3._wp*(strength_idx + 1._wp))
                end do
            else

                !< 2D cartesian function:
                ! We smear particles considering a virtual depth (lag_params%charwidth)
                theta = 0._wp
                Nr = ceiling(lag_params%charwidth/(y_cb(cellaux(2)) - y_cb(cellaux(2) - 1)))
                Nr_count = 1._wp - mapCells*1._wp
                !dzp = y_cb(cellaux(2) + 1) - y_cb(cellaux(2))
                dzp = dy(cellaux(2))
                Lz2 = (center(3) - (dzp*(0.5_wp + Nr_count) - lag_params%charwidth/2._wp))**2._wp
                distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + Lz2)
                func = dzp/lag_params%charwidth*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
                do while (Nr_count < Nr - 1._wp + ((mapCells - 1)*1._wp))
                    Nr_count = Nr_count + 1._wp
                    Lz2 = (center(3) - (dzp*(0.5_wp + Nr_count) - lag_params%charwidth/2._wp))**2._wp
                    distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + Lz2)
                    func = func + &
                           dzp/lag_params%charwidth*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**(3._wp*(strength_idx + 1._wp))
                end do
            end if
        end if

    end subroutine s_applygaussian

    subroutine s_gaussian_hifu(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(vector_field), intent(inout) :: updatedvar
        real(wp), dimension(1:lag_params%nBubs_glb), intent(in), optional :: lbk_qvis, lbk_qth

        real(wp), dimension(3) :: center
        integer, dimension(3) :: cell
        real(wp) :: stddsv
        real(wp) :: strength_vel, strength_vol, gpu_sum

        real(wp), dimension(3) :: nodecoord, scoord
        real(wp) :: addFun1, addFun2, addFun3, sumFun
        real(wp) :: func, func2, volpart, normGaussSum
        integer, dimension(3) :: cellaux
        integer :: l, i, j, k
        logical :: celloutside, particle_in_domain
        integer :: smearGrid, smearGridz

        smearGrid = mapCells - (-mapCells) ! Include the cell that contains the bubble (3+1+3)

        if (hifu_params%cartesian) then

            !$acc parallel loop gang vector default(present) private(l, cell, scoord, center) &
            !$acc copyin(smearGrid)
            do l = 1, nBubs

                !> Is the particle in the domain?
                particle_in_domain = ((lbk_pos(l, 1, 1) < x_cb_hf(m_hf + buff_size)) .and. (lbk_pos(l, 1, 1) >= x_cb_hf(-1 - buff_size)) .and. &
                                    (lbk_pos(l, 2, 1) < y_cb_hf(n_hf + buff_size)) .and. (lbk_pos(l, 2, 1) >= y_cb_hf(-1 - buff_size)) .and. &
                                    (lbk_pos(l, 3, 1) < z_cb_hf(p_hf + buff_size)) .and. (lbk_pos(l, 3, 1) >= z_cb_hf(-1 - buff_size)))

                if (particle_in_domain) then

                    !> Find cell of the bubble in the heat solver domain
                    cell = -buff_size
                    !x
                    do while (lbk_pos(l, 1, 1) < x_cb_hf(cell(1) - 1))
                        cell(1) = cell(1) - 1
                    end do
                    do while (lbk_pos(l, 1, 1) > x_cb_hf(cell(1)))
                        cell(1) = cell(1) + 1
                    end do
                    !y
                    do while (lbk_pos(l, 2, 1) < y_cb_hf(cell(2) - 1))
                        cell(2) = cell(2) - 1
                    end do
                    do while (lbk_pos(l, 2, 1) > y_cb_hf(cell(2)))
                        cell(2) = cell(2) + 1
                    end do
                    !z
                    do while (lbk_pos(l, 3, 1) < z_cb_hf(cell(3) - 1))
                        cell(3) = cell(3) - 1
                    end do
                    do while (lbk_pos(l, 3, 1) > z_cb_hf(cell(3)))
                        cell(3) = cell(3) + 1
                    end do
                    !coordinates in computational space
                    scoord(1) = cell(1) + (lbk_pos(l, 1, 1) - x_cb_hf(cell(1) - 1))/dx_hf(cell(1))
                    scoord(2) = cell(2) + (lbk_pos(l, 2, 1) - y_cb_hf(cell(2) - 1))/dy_hf(cell(2))
                    scoord(3) = cell(3) + (lbk_pos(l, 3, 1) - z_cb_hf(cell(3) - 1))/dz_hf(cell(3))
                    !cell
                    cell(:) = int(scoord(:))

                    !> Gasussian parameters and cells to smear
                    volpart = 4._wp/3._wp*pi*lbk_rad(l, 1)**3._wp
                    call s_compute_stddsv(cell, volpart, stddsv)
                    center(1:3) = lbk_pos(l, 1:3, 2)

                    !> Smearing
                    normGaussSum = 0._wp
                    !$acc loop collapse(3) gang vector private(cellaux, nodecoord) reduction(+: normGaussSum)
                    do i = 0, smearGrid
                        do j = 0, smearGrid
                            do k = 0, smearGrid

                                cellaux(1) = cell(1) + i - mapCells
                                cellaux(2) = cell(2) + j - mapCells
                                cellaux(3) = cell(3) + k - mapCells

                                !> Check if the cells intended to smear the bubbles in are in the computational domain (heat solver)
                                celloutside = .false.
                                if ((cellaux(3) < -buff_size) .or. (cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                    celloutside = .true.
                                end if
                                if ((cellaux(3) > p_hf + buff_size) .or. (cellaux(2) > n_hf + buff_size) .or. (cellaux(1) > m_hf + buff_size)) then
                                    celloutside = .true.
                                end if

                                if (.not. celloutside) then
                                    nodecoord(1) = x_cc_hf(cellaux(1))
                                    nodecoord(2) = y_cc_hf(cellaux(2))
                                    nodecoord(3) = z_cc_hf(cellaux(3))
                                    call s_applygaussian(center, cellaux, nodecoord, stddsv, 0._wp, func)
                                    !func = func / sumFun !Adjusted intensity
                                    ! Relocate cells for bubbles intersecting symmetric boundaries
                                    !if (any((/bcxb, bcxe, bcyb, bcye, bczb, bcze/) == BC_REFLECTIVE)) then
                                    !    call s_shift_cell_symmetric_bc(cellaux, cell)
                                    !end if

                                else
                                    func = 0._wp
                                    cellaux(1) = cell(1)
                                    cellaux(2) = cell(2)
                                    cellaux(3) = cell(3)
                                end if

                                !Summation of the normalized gaussian function
                                normGaussSum = normGaussSum + func*(dx_hf(cellaux(1))*dy_hf(cellaux(2))*dz_hf(cellaux(3)))

                                !Update qvis field
                                addFun1 = func*lbk_qvis(l)
                                !$acc atomic update
                                updatedvar%vf(hifu_params%qvis_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar%vf(hifu_params%qvis_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + addFun1

                                !Update qth field
                                addFun2 = func*lbk_qth(l)
                                !$acc atomic update
                                updatedvar%vf(hifu_params%qth_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar%vf(hifu_params%qth_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + addFun2

                            end do
                        end do
                    end do

                    ! Summation of normal gaussian weights must be equal to one, except from the cells at the buffers since some surrounding cells can be outside the domain.
                    ! Tolerance of 0.01 defined
                    if ((cell(1) > 0 .and. cell(1) <= m_hf) .and. (cell(2) > 0 .and. cell(2) <= n_hf) .and. &
                        (cell(3) > 0 .and. cell(3) <= p_hf) .and. abs(normGaussSum - 1._wp) > 0.1_wp) then
                        print *, 'Smeared bubble out of tolerance:', l, normGaussSum, cell(1), cell(2), cell(3)
                    end if

                end if

            end do

        elseif (p>0 .and. .not. cyl_coord) then

            !$acc parallel loop gang vector default(present) private(l, cell, scoord, center) &
            !$acc copyin(smearGrid)
            do l = 1, nBubs

                !> Is the particle in the domain?
                particle_in_domain = ((lbk_pos(l, 1, 1) < x_cb(m + buff_size)) .and. (lbk_pos(l, 1, 1) >= x_cb(-1 - buff_size)) .and. &
                                    (lbk_pos(l, 2, 1) < y_cb(n + buff_size)) .and. (lbk_pos(l, 2, 1) >= y_cb(-1 - buff_size)) .and. &
                                    (lbk_pos(l, 3, 1) < z_cb(p + buff_size)) .and. (lbk_pos(l, 3, 1) >= z_cb(-1 - buff_size)))

                if (particle_in_domain) then

                    !> Find cell of the bubble in the heat solver domain
                    cell = -buff_size
                    !x
                    do while (lbk_pos(l, 1, 1) < x_cb(cell(1) - 1))
                        cell(1) = cell(1) - 1
                    end do
                    do while (lbk_pos(l, 1, 1) > x_cb(cell(1)))
                        cell(1) = cell(1) + 1
                    end do
                    !y
                    do while (lbk_pos(l, 2, 1) < y_cb(cell(2) - 1))
                        cell(2) = cell(2) - 1
                    end do
                    do while (lbk_pos(l, 2, 1) > y_cb(cell(2)))
                        cell(2) = cell(2) + 1
                    end do
                    !z
                    do while (lbk_pos(l, 3, 1) < z_cb(cell(3) - 1))
                        cell(3) = cell(3) - 1
                    end do
                    do while (lbk_pos(l, 3, 1) > z_cb(cell(3)))
                        cell(3) = cell(3) + 1
                    end do
                    !coordinates in computational space
                    scoord(1) = cell(1) + (lbk_pos(l, 1, 1) - x_cb(cell(1) - 1))/dx(cell(1))
                    scoord(2) = cell(2) + (lbk_pos(l, 2, 1) - y_cb(cell(2) - 1))/dy(cell(2))
                    scoord(3) = cell(3) + (lbk_pos(l, 3, 1) - z_cb(cell(3) - 1))/dz(cell(3))
                    !cell
                    cell(:) = int(scoord(:))

                    !> Gasussian parameters and cells to smear
                    volpart = 4._wp/3._wp*pi*lbk_rad(l, 1)**3._wp
                    call s_compute_stddsv(cell, volpart, stddsv)
                    center(1:3) = lbk_pos(l, 1:3, 2)
                    
                    if (l==1) print*, 'in kernel', lbk_qvis(l), lbk_qth(l)

                    !> Smearing
                    normGaussSum = 0._wp
                    !$acc loop collapse(3) gang vector private(cellaux, nodecoord) reduction(+: normGaussSum)
                    do i = 0, smearGrid
                        do j = 0, smearGrid
                            do k = 0, smearGrid

                                cellaux(1) = cell(1) + i - mapCells
                                cellaux(2) = cell(2) + j - mapCells
                                cellaux(3) = cell(3) + k - mapCells

                                !> Check if the cells intended to smear the bubbles in are in the computational domain (heat solver)
                                celloutside = .false.
                                if ((cellaux(3) < -buff_size) .or. (cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                    celloutside = .true.
                                end if
                                if ((cellaux(3) > p + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                                    celloutside = .true.
                                end if

                                if (.not. celloutside) then
                                    nodecoord(1) = x_cc(cellaux(1))
                                    nodecoord(2) = y_cc(cellaux(2))
                                    nodecoord(3) = z_cc(cellaux(3))
                                    call s_applygaussian(center, cellaux, nodecoord, stddsv, 0._wp, func)

                                    ! Relocate cells for bubbles intersecting symmetric boundaries
                                    if (any((/bcxb, bcxe, bcyb, bcye, bczb, bcze/) == BC_REFLECTIVE)) then
                                        call s_shift_cell_symmetric_bc(cellaux, cell)
                                    end if

                                else
                                    func = 0._wp
                                    cellaux(1) = cell(1)
                                    cellaux(2) = cell(2)
                                    cellaux(3) = cell(3)
                                end if

                                !Summation of the normalized gaussian function
                                normGaussSum = normGaussSum + func*(dx(cellaux(1))*dy(cellaux(2))*dz(cellaux(3)))

                                !Update qvis field
                                addFun1 = func*lbk_qvis(l)
                                !$acc atomic update
                                updatedvar%vf(hifu_params%qvis_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar%vf(hifu_params%qvis_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + addFun1

                                !Update qth field
                                addFun2 = func*lbk_qth(l)
                                !$acc atomic update
                                updatedvar%vf(hifu_params%qth_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar%vf(hifu_params%qth_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + addFun2

                            end do
                        end do
                    end do

                    ! Summation of normal gaussian weights must be equal to one, except from the cells at the buffers since some surrounding cells can be outside the domain.
                    ! Tolerance of 0.01 defined
                    if ((cell(1) > 0 .and. cell(1) <= m) .and. (cell(2) > 0 .and. cell(2) <= n) .and. &
                        (cell(3) > 0 .and. cell(3) <= p) .and. abs(normGaussSum - 1._wp) > 0.1_wp) then
                        print *, 'Smeared bubble out of tolerance:', l, normGaussSum, cell(1), cell(2), cell(3)
                    end if

                end if

            end do

            ! Populate symmetric boundaries
            if (any((/bcxb, bcxe, bcyb, bcye, bczb, bcze/) == BC_REFLECTIVE)) then
                call s_populate_symmetric_bc(updatedvar)
            end if

        end if

        ! Populate symmetric boundaries
        !if (any((/bcxb, bcxe, bcyb, bcye, bczb, bcze/) == BC_REFLECTIVE)) then
        !    call s_populate_symmetric_bc(updatedvar)
        !end if

        if (proc_rank == 0) print *, 'Bubble sources smeared with Gaussian kernel: qvis & qth'

    end subroutine s_gaussian_hifu

    !> The purpose of this procedure contains the algorithm to use the delta kernel function to map the effect of the bubbles.
            !!      The effect of the bubbles only affects the cell where the bubble is located.
    subroutine s_deltafunc_hifu(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, lbk_qvis, lbk_qth)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(vector_field), intent(inout) :: updatedvar
        real(wp), dimension(1:lag_params%nBubs_glb), intent(in), optional :: lbk_qvis, lbk_qth

        integer, dimension(3) :: cell
        real(wp), dimension(3) :: scoord
        real(wp) :: addFun1, addFun2
        real(wp) :: volCell
        integer :: l
        logical :: particle_in_domain

        !$acc parallel loop gang vector default(present) private(l, cell, scoord)
        do l = 1, nBubs

            !> Is the particle in the domain?
            particle_in_domain = ((lbk_pos(l, 1, 1) < x_cb_hf(m_hf + buff_size)) .and. (lbk_pos(l, 1, 1) >= x_cb_hf(-1 - buff_size)) .and. &
                                  (lbk_pos(l, 2, 1) < y_cb_hf(n_hf + buff_size)) .and. (lbk_pos(l, 2, 1) >= y_cb_hf(-1 - buff_size)) .and. &
                                  (lbk_pos(l, 3, 1) < z_cb_hf(p_hf + buff_size)) .and. (lbk_pos(l, 3, 1) >= z_cb_hf(-1 - buff_size)))

            if (particle_in_domain) then

                !> Find cell of the bubble in the heat solver domain
                cell = -buff_size
                !x
                do while (lbk_pos(l, 1, 1) < x_cb_hf(cell(1) - 1))
                    cell(1) = cell(1) - 1
                end do
                do while (lbk_pos(l, 1, 1) > x_cb_hf(cell(1)))
                    cell(1) = cell(1) + 1
                end do
                !y
                do while (lbk_pos(l, 2, 1) < y_cb_hf(cell(2) - 1))
                    cell(2) = cell(2) - 1
                end do
                do while (lbk_pos(l, 2, 1) > y_cb_hf(cell(2)))
                    cell(2) = cell(2) + 1
                end do
                !z
                do while (lbk_pos(l, 3, 1) < z_cb_hf(cell(3) - 1))
                    cell(3) = cell(3) - 1
                end do
                do while (lbk_pos(l, 3, 1) > z_cb_hf(cell(3)))
                    cell(3) = cell(3) + 1
                end do
                !coordinates in computational space
                scoord(1) = cell(1) + (lbk_pos(l, 1, 1) - x_cb_hf(cell(1) - 1))/dx_hf(cell(1))
                scoord(2) = cell(2) + (lbk_pos(l, 2, 1) - y_cb_hf(cell(2) - 1))/dy_hf(cell(2))
                scoord(3) = cell(3) + (lbk_pos(l, 3, 1) - z_cb_hf(cell(3) - 1))/dz_hf(cell(3))
                !cell
                cell(:) = int(scoord(:))

                !> Volume of the cell
                volCell = dx_hf(cell(1))*dy_hf(cell(2))*dz_hf(cell(3))

                !> Smearing
                !Update qvis field
                addFun1 = lbk_qvis(l)/volCell
                !$acc atomic update
                updatedvar%vf(hifu_params%qvis_idx)%sf(cell(1), cell(2), cell(3)) = &
                    updatedvar%vf(hifu_params%qvis_idx)%sf(cell(1), cell(2), cell(3)) &
                    + addFun1

                !Update qth field
                addFun2 = lbk_qth(l)/volCell
                !$acc atomic update
                updatedvar%vf(hifu_params%qth_idx)%sf(cell(1), cell(2), cell(3)) = &
                    updatedvar%vf(hifu_params%qth_idx)%sf(cell(1), cell(2), cell(3)) &
                    + addFun2

            end if
        end do

        if (proc_rank == 0) print *, 'Bubble sources smeared with delta kernel: qvis & qth'

    end subroutine s_deltafunc_hifu

    !> The purpose of this subroutine is to check if the current cell is outside the computational domain or not (including ghost cells).
            !! @param cellaux Tested cell to smear the bubble effect in.
            !! @param celloutside If true, then cellaux is outside the computational domain.
    subroutine s_check_celloutside(center, cellaux, nodecoord, celloutside)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_check_celloutside
#else
        !$acc routine seq
#endif
        real(wp), dimension(3), intent(in) :: center, nodecoord
        integer, dimension(3), intent(inout) :: cellaux
        logical, intent(out) :: celloutside

        real(wp) :: distance, vrtDist, chardist

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

        if (lag_params%newModel_2D .and. .not. celloutside) then
            distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + (center(3) - nodecoord(3))**2._wp)
            vrtDist = 0.5_wp*(dx(cellaux(1)) + dy(cellaux(2)))
            chardist = (dx(cellaux(1))*dy(cellaux(2))*vrtDist)**(1._wp/3._wp)
            if (distance >= 5._wp*chardist) celloutside = .true.
            return
        end if

    end subroutine s_check_celloutside

    !> This subroutine relocates the current cell, if it intersects a symmetric boundary.
            !! @param cell Cell of the current bubble
            !! @param cellaux Cell to map the bubble effect in.
    subroutine s_shift_cell_symmetric_bc(cellaux, cell)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_shift_cell_symmetric_bc
#else
        !$acc routine seq
#endif
        integer, dimension(3), intent(inout) :: cellaux
        integer, dimension(3), intent(in) :: cell

        ! x-dir
        if (bcxb == BC_REFLECTIVE .and. (cell(1) <= mapCells - 1)) then
            cellaux(1) = abs(cellaux(1)) - 1
        end if
        if (bcxe == BC_REFLECTIVE .and. (cell(1) >= m + 1 - mapCells)) then
            cellaux(1) = cellaux(1) - (2*(cellaux(1) - m) - 1)
        end if

        !y-dir
        if (bcyb == BC_REFLECTIVE .and. (cell(2) <= mapCells - 1)) then
            cellaux(2) = abs(cellaux(2)) - 1
        end if
        if (bcye == BC_REFLECTIVE .and. (cell(2) >= n + 1 - mapCells)) then
            cellaux(2) = cellaux(2) - (2*(cellaux(2) - n) - 1)
        end if

        if (p > 0) then
            !z-dir
            if (bczb == BC_REFLECTIVE .and. (cell(3) <= mapCells - 1)) then
                cellaux(3) = abs(cellaux(3)) - 1
            end if
            if (bcze == BC_REFLECTIVE .and. (cell(3) >= p + 1 - mapCells)) then
                cellaux(3) = cellaux(3) - (2*(cellaux(3) - p) - 1)
            end if
        end if

    end subroutine s_shift_cell_symmetric_bc

    subroutine s_populate_symmetric_bc(updatedvar)

        type(vector_field), intent(inout) :: updatedvar

        integer :: j, k, l

        ! x-dir
        if (bcxb == BC_REFLECTIVE) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 1, buff_size
                        updatedvar%vf(1)%sf(-j, k, l) = updatedvar%vf(1)%sf(j - 1, k, l)
                        updatedvar%vf(2)%sf(-j, k, l) = updatedvar%vf(2)%sf(j - 1, k, l)
                        if (p == 0) then
                            updatedvar%vf(5)%sf(-j, k, l) = updatedvar%vf(5)%sf(j - 1, k, l)
                        end if
                    end do
                end do
            end do
        end if
        if (bcxe == BC_REFLECTIVE) then
            !$acc parallel loop collapse(3) default(present)
            do l = 0, p
                do k = 0, n
                    do j = 1, buff_size
                        updatedvar%vf(1)%sf(m + j, k, l) = updatedvar%vf(1)%sf(m - (j - 1), k, l)
                        updatedvar%vf(2)%sf(m + j, k, l) = updatedvar%vf(2)%sf(m - (j - 1), k, l)
                        if (p == 0) then
                            updatedvar%vf(5)%sf(m + j, k, l) = updatedvar%vf(5)%sf(m - (j - 1), k, l)
                        end if
                    end do
                end do
            end do
        end if

        !y-dir
        if (bcyb == BC_REFLECTIVE) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        updatedvar%vf(1)%sf(l, -j, k) = updatedvar%vf(1)%sf(l, j - 1, k)
                        updatedvar%vf(2)%sf(l, -j, k) = updatedvar%vf(2)%sf(l, j - 1, k)
                        if (p == 0) then
                            updatedvar%vf(5)%sf(l, -j, k) = updatedvar%vf(5)%sf(l, j - 1, k)
                        end if
                    end do
                end do
            end do
        end if
        if (bcye == BC_REFLECTIVE) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        updatedvar%vf(1)%sf(l, n + j, k) = updatedvar%vf(1)%sf(l, n - (j - 1), k)
                        updatedvar%vf(2)%sf(l, n + j, k) = updatedvar%vf(2)%sf(l, n - (j - 1), k)
                        if (p == 0) then
                            updatedvar%vf(5)%sf(l, n + j, k) = updatedvar%vf(5)%sf(l, n - (j - 1), k)
                        end if
                    end do
                end do
            end do
        end if

        if (p > 0) then
            !z-dir
            if (bczb == BC_REFLECTIVE) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            updatedvar%vf(1)%sf(k, l, -j) = updatedvar%vf(1)%sf(k, l, j - 1)
                            updatedvar%vf(2)%sf(k, l, -j) = updatedvar%vf(2)%sf(k, l, j - 1)
                            if (p == 0) then
                                updatedvar%vf(5)%sf(k, l, -j) = updatedvar%vf(5)%sf(k, l, j - 1)
                            end if
                        end do
                    end do
                end do
            end if
            if (bcze == BC_REFLECTIVE) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            updatedvar%vf(1)%sf(k, l, p + j) = updatedvar%vf(1)%sf(k, l, p - (j - 1))
                            updatedvar%vf(2)%sf(k, l, p + j) = updatedvar%vf(2)%sf(k, l, p - (j - 1))
                            if (p == 0) then
                                updatedvar%vf(5)%sf(k, l, p + j) = updatedvar%vf(5)%sf(k, l, p - (j - 1))
                            end if
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_populate_symmetric_bc

    !> Calculates the standard deviation of the bubble being smeared in the Eulerian framework.
            !! @param cell Cell where the bubble is located
            !! @param volpart Volume of the bubble
            !! @param stddsv Standard deviaton
    subroutine s_compute_stddsv(cell, volpart, stddsv)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_stddsv
#else
        !$acc routine seq
#endif
        integer, dimension(3), intent(in) :: cell
        real(wp), intent(in) :: volpart
        real(wp), intent(out) :: stddsv

        real(wp) :: chardist, charvol, vrtDist
        real(wp) :: rad

        if (hifu_params%cartesian .and. hifu_params%heatSolver) then

            !< Compute characteristic distance
            chardist = (dx_hf(cell(1))*dy_hf(cell(2))*dz_hf(cell(3)))**(1._wp/3._wp)

            !< Compute Standard deviaton
            rad = (3._wp*volpart/(4._wp*pi))**(1._wp/3._wp)
            stddsv = max(chardist, rad)

        else

            !< Compute characteristic distance
            chardist = sqrt(dx(cell(1))*dy(cell(2)))
            if (p > 0) chardist = (dx(cell(1))*dy(cell(2))*dz(cell(3)))**(1._wp/3._wp)

            !< Compute characteristic volume
            if (p > 0) then
                charvol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            else
                if (cyl_coord) then
                    charvol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
                else
                    charvol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                end if
            end if

            if (lag_params%newModel_2D) then
                vrtDist = 0.5_wp*(dx(cell(1)) + dy(cell(2)))
                chardist = (dx(cell(1))*dy(cell(2))*vrtDist)**(1._wp/3._wp)
                charvol = dx(cell(1))*dy(cell(2))*vrtDist
            end if

            !< Compute Standard deviaton
            if (((volpart/charvol) > 0.5_wp*lag_params%valmaxvoid) .or. (lag_params%smooth_type == 1)) then
                rad = (3._wp*volpart/(4._wp*pi))**(1._wp/3._wp)
                stddsv = 1._wp*lag_params%epsilonb*max(chardist, rad)
                !print*, rad, dx(cell(1)), dy(cell(2)), y_cc(cell(2)), cell(1), cell(2)
            else
                stddsv = 0._wp
            end if

        end if

    end subroutine s_compute_stddsv

    !> The purpose of this procedure is to calculate the characteristic cell volume
            !! @param cell Computational coordinates (x, y, z)
            !! @param Charvol Characteristic volume
    subroutine s_get_char_vol(cellx, celly, cellz, Charvol)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_get_char_vol
#else
        !$acc routine seq
#endif
        integer, intent(in) :: cellx, celly, cellz
        real(wp), intent(out) :: Charvol

        if (p > 0) then
            Charvol = dx(cellx)*dy(celly)*dz(cellz)
        else
            if (cyl_coord) then
                Charvol = dx(cellx)*dy(celly)*y_cc(celly)*2._wp*pi
            else
                Charvol = dx(cellx)*dy(celly)*lag_params%charwidth
            end if
        end if

        !print*, cellx, celly, cellz, Charvol

    end subroutine s_get_char_vol

    !> This subroutine transforms the computational coordinates of the bubble from
            !!      real type into integer.
            !! @param s Computational coordinates of the bubble, real type
            !! @param get_cell Computational coordinates of the bubble, integer type
    subroutine s_get_cell(s_cell, get_cell)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_get_cell
#else
        !$acc routine seq
#endif
        real(wp), dimension(3), intent(in) :: s_cell
        integer, dimension(3), intent(out) :: get_cell
        integer :: i

        get_cell(:) = int(s_cell(:))
        do i = 1, num_dims
            if (s_cell(i) < 0._wp) get_cell(i) = get_cell(i) - 1
        end do

    end subroutine s_get_cell

end module m_bubbles_EL_kernels
