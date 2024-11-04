!>
!! @file m_kernel_functions.f90
!! @brief Contains module m_kernel_functions

#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian bubbles
!!        in the Eulerian framework.
module m_kernel_functions

    ! Dependencies =============================================================

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! ==========================================================================

    implicit none

contains

    !> The purpose of this procedure is to smear the strength of the lagrangian
        !!      bubbles into the Eulerian framework.
        !! @param updatedvar Variable to be updated
        !! @param center Coordinates of the lagragian bubble
        !! @param cell  Computational coordinates of the cell that contains the bubble.
        !! @param strength Variable to be smeared
        !! @param kernelused Kernel smoothening function
        !! @param stddsv Standard deviation of the bubble contribution in the Eulerian space
        !! @param strength2 Define the product of two smeared functions (strength*strength2)
    subroutine s_smoothfunction(updatedvar, center, cell, strength, kernelused, stddsv, strength2)

        type(scalar_field) :: updatedvar
        real(kind(0.d0)), dimension(3) :: center
        integer, dimension(3) :: cell
        real(kind(0.d0)) :: strength, stddsv
        integer :: kernelused
        real(kind(0.d0)), optional :: strength2

        smoothfunc:select case(kernelused)
        case (1)
        call s_gaussian(updatedvar, center, cell, strength, stddsv, strength2)
        case (2)
        call s_deltafunc(updatedvar, cell, strength)
        end select smoothfunc

    end subroutine s_smoothfunction

    !> The purpose of this procedure is to smear the strength of the bubbles
        !!      in the Eulerian framework using a gaussian kernel function.
    subroutine s_gaussian(updatedvar, center, cell, strength, stddsv, strength2)

        real(kind(0.d0)), dimension(3) :: center
        real(kind(0.d0)) :: strength, stddsv
        integer, dimension(3) :: cell, epsilonbaux
        type(scalar_field) :: updatedvar
        real(kind(0.d0)), optional :: strength2

        ! For now conservative for restart when epsilonbaux < buff_size
        epsilonbaux(:) = 3

        call s_applygaussian(updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2)

        !For symmetric BC
        if ((bc_x%beg == -2 .or. bc_x%end == -2 .or. bc_y%beg == -2 .or. bc_y%end == -2 .or. bc_z%beg == -2 .or. &
             bc_z%end == -2) .or. (bc_x%beg == proc_rank .or. bc_x%end == proc_rank .or. bc_y%beg == proc_rank .or. &
                                   bc_y%end == proc_rank .or. bc_z%beg == proc_rank .or. bc_z%end == proc_rank)) then

            call s_gaussian_symmetric_bc(updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2)

        end if

    end subroutine s_gaussian

    !> This subroutine smeares the strength of the fictitious bubbles, across the
        !!      boundary when symmetric boundary condition is used.
    subroutine s_gaussian_symmetric_bc(updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2)

        real(kind(0.d0)), dimension(3) :: center, centertmp
        real(kind(0.d0)) :: strength, stddsv
        integer, dimension(3) :: cell, celltmp, epsilonbaux
        type(scalar_field) :: updatedvar
        real(kind(0.d0)), optional :: strength2

        celltmp = cell
        centertmp = center

        if (((bc_x%beg == -2) .or. (bc_x%beg == proc_rank)) .and. cell(1) + 1 <= epsilonbaux(1)) then

            call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 1)
            call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
            !Apply(x-d,y,z)

            if (((bc_y%beg == -2) .or. (bc_y%beg == proc_rank)) .and. cell(2) + 1 <= epsilonbaux(2) .and. (.not. cyl_coord)) then

                call s_kernel_shift_bc(center, centertmp, cell, celltmp, 2, 1)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x-d,y-d,z)

                centertmp(1) = center(1)
                celltmp(1) = cell(1)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x,y-d,z)

                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y-d,z-d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y-d,z-d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y,z-d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y-d,z+d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y-d,z+d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y,z+d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if

                end if

            else if (((bc_y%end == -2) .or. (bc_y%end == proc_rank)) .and. cell(2) >= n + 1 - epsilonbaux(2) .and. (.not. cyl_coord)) then

                call s_kernel_shift_bc(center, centertmp, cell, celltmp, 2, 2)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x-d,y+d,z)

                centertmp(1) = center(1)
                celltmp(1) = cell(1)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x,y+d,z)

                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y+d,z-d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y+d,z-d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y,z-d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y+d,z+d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y+d,z+d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y,z+d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            else
                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y,z-d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) + 1 >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x-d,y,z+d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            end if
        else if (((bc_x%end == -2) .or. (bc_x%end == proc_rank)) .and. cell(1) >= m + 1 - epsilonbaux(1)) then

            call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 2)
            call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
            !Apply(x+d,y,z)

            if (((bc_y%beg == -2) .or. (bc_y%beg == proc_rank)) .and. cell(2) + 1 <= epsilonbaux(2) .and. (.not. cyl_coord)) then

                call s_kernel_shift_bc(center, centertmp, cell, celltmp, 2, 1)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x+d,y-d,z)

                centertmp(1) = center(1)
                celltmp(1) = cell(1)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x,y-d,z)

                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y-d,z-d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y-d,z-d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y,z-d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y-d,z+d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y-d,z+d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y,z+d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            else if (((bc_y%end == -2) .or. (bc_y%end == proc_rank)) .and. cell(2) >= n + 1 - epsilonbaux(2) .and. (.not. cyl_coord)) then

                call s_kernel_shift_bc(center, centertmp, cell, celltmp, 2, 2)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x+d,y+d,z)

                centertmp(1) = center(1)
                celltmp(1) = cell(1)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x,y+d,z)

                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y+d,z-d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y+d,z-d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y,z-d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) + 1 >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y+d,z+d)

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 1, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y+d,z+d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y,z+d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            else
                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y,z-d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x+d,y,z+d)

                        centertmp(1) = center(1)
                        celltmp(1) = cell(1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            end if
        else
            if (((bc_y%beg == -2) .or. (bc_y%beg == proc_rank)) .and. cell(2) + 1 <= epsilonbaux(2) .and. (.not. cyl_coord)) then

                call s_kernel_shift_bc(center, centertmp, cell, celltmp, 2, 1)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x,y-d,z)

                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y-d,z-d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y-d,z+d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            else if (((bc_y%end == -2) .or. (bc_y%end == proc_rank)) .and. cell(2) >= n + 1 - epsilonbaux(2) .and. (.not. cyl_coord)) then
                call s_kernel_shift_bc(center, centertmp, cell, celltmp, 2, 2)
                call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                !Apply(x,y+d,z)

                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y+d,z-d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y+d,z+d)

                        centertmp(2) = center(2)
                        celltmp(2) = cell(2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            else
                if (p > 0) then
                    if (((bc_z%beg == -2) .or. (bc_z%beg == proc_rank)) .and. cell(3) + 1 <= epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 1)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z-d)

                    else if (((bc_z%end == -2) .or. (bc_z%end == proc_rank)) .and. cell(3) >= p + 1 - epsilonbaux(3)) then

                        call s_kernel_shift_bc(center, centertmp, cell, celltmp, 3, 2)
                        call s_applygaussian(updatedvar, centertmp, celltmp, strength, epsilonbaux, stddsv, strength2)
                        !Apply(x,y,z+d)

                    end if
                end if
            end if
        end if

    end subroutine s_gaussian_symmetric_bc

    !> The purpose of this procedure is to apply the gaussian kernel function.
    subroutine s_applygaussian(updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2)

        real(kind(0.d0)), dimension(3) :: nodecoord, center, auxvect
        real(kind(0.d0)) :: strength, func, distance, stddsv
        integer, dimension(3) :: cell, cellaux, cellaux2
        integer :: idir, idir2, i, j, k, ini, iend, jini, jend
        integer, dimension(3) :: epsilonbaux
        type(scalar_field) :: updatedvar
        logical :: celloutside
        real(kind(0.d0)) :: theta, dtheta, L2, dz, Lz2
        integer :: Nr, Nr_count
        real(kind(0.d0)), optional :: strength2 ! for product of two smeared functions

        nodecoord = 0

        if (num_dims == 3) then
            k = -epsilonbaux(3)
        else
            k = 0
        end if

        i = -epsilonbaux(1); j = -epsilonbaux(2)

3001    if ((i <= epsilonbaux(1)) .and. (j <= epsilonbaux(2))) then

            celloutside = .false.

            cellaux(1) = cell(1) + i
            cellaux(2) = cell(2) + j
            cellaux(3) = cell(3) + k

            if (cellaux(1) < -buff_size) then
                celloutside = .true.
                i = i + 1
            end if

            !check ghost part in y direction
            if (cellaux(2) < -buff_size) then
                celloutside = .true.
                j = j + 1
            end if

            if (cyl_coord .and. num_dims /= 3) then
                if ((cellaux(2) < n+buff_size) .and. (.not. celloutside)) then 
                    if (y_cc(cellaux(2)) < 0d0) then
                        celloutside = .true.
                        j = j + 1
                    end if
                end if
            end if

            ! Temp
            if (num_dims == 3) then
                if (cellaux(3) < -buff_size) then
                    celloutside = .true.
                    k = k + 1
                end if
            end if

            if (celloutside) goto 3001
            ! Temp
            if (cellaux(3) > p + buff_size) return
            if (cellaux(2) > n + buff_size) goto 3002
            if (cellaux(1) > m + buff_size) goto 3003 !return

            !z direction, periodic
            nodecoord(1) = x_cc(cellaux(1))
            nodecoord(2) = y_cc(cellaux(2))
            if (p > 0) nodecoord(3) = z_cc(cellaux(3))
            auxvect(:) = center(:) - nodecoord(:)
            distance = sqrt(auxvect(1)**2 + auxvect(2)**2 + auxvect(3)**2)
            if (num_dims == 3) then
                func = exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
            else

                ! for 2D cylindrical coordinate we smear particles in the azimuthal
                ! direction for given r
                if (cyl_coord) then
                    theta = 0d0
                    Nr = ceiling(2d0*PI*nodecoord(2)/(y_cb(cellaux(2)) - y_cb(cellaux(2) - 1))) ! number of cells in the 3D ring
                    dtheta = 2d0*PI/Nr ! dtheta of the ring divisions
                    L2 = center(2)**2d0 + nodecoord(2)**2d0 - 2d0*center(2)*nodecoord(2)*cos(theta)
                    distance = DSQRT(auxvect(1)**2d0 + L2) ! Distance from bubble to axisymmetric plane
                    ! Factor 2d0 is for symmetry (upper half of the 2D field (+r) is considered)
                    func = dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
                    Nr_count = 0
                    do while (Nr_count < Nr - 1)
                        Nr_count = Nr_count + 1
                        theta = Nr_count*dtheta
                        ! trigonometric relation
                        L2 = center(2)**2d0 + nodecoord(2)**2d0 - 2d0*center(2)*nodecoord(2)*cos(theta)
                        distance = DSQRT(auxvect(1)**2d0 + L2)
                        ! nodecoord(2)*dtheta is the azimuthal width of the cell
                        if (present(strength2)) then
                            ! Product of two Gaussians
                            func = func + &
                                   dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**6
                        else
                            func = func + &
                                   dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 ! eqn. 60 to map bubble volumet onto the axisymmetric grid
                        end if
                    end do

                    ! 2D with virtual depth
                else

                    theta = 0d0
                    Nr = ceiling(charwidth/(y_cb(cellaux(2)) - y_cb(cellaux(2) - 1)))
                    Nr_count = 1 - epsilonbaux(3)
                    dz = y_cb(cellaux(2) + 1) - y_cb(cellaux(2))
                    Lz2 = (center(3) - (dz*(5d-1 + Nr_count) - charwidth/2d0))**2d0
                    distance = DSQRT(auxvect(1)**2d0 + auxvect(2)**2d0 + Lz2)
                    func = dz/charwidth*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
                    do while (Nr_count < Nr - 1 + (epsilonbaux(3) - 1))
                        Nr_count = Nr_count + 1
                        Lz2 = (center(3) - (dz*(5d-1 + Nr_count) - charwidth/2d0))**2d0
                        distance = DSQRT(auxvect(1)**2d0 + auxvect(2)**2d0 + Lz2)
                        if (present(strength2)) then
                            ! Product of two Gaussians
                            func = func + &
                                   dz/charwidth*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**6
                        else
                            func = func + &
                                   dz/charwidth*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
                        end if
                    end do
                end if
            end if

            if (present(strength2)) then
                updatedvar%sf(cellaux(1), cellaux(2), cellaux(3)) = updatedvar%sf(cellaux(1), cellaux(2), cellaux(3)) + func*strength*strength2
            else
                updatedvar%sf(cellaux(1), cellaux(2), cellaux(3)) = updatedvar%sf(cellaux(1), cellaux(2), cellaux(3)) + func*strength
            end if

            if (j < epsilonbaux(2)) then
                j = j + 1
                goto 3001
            end if

3002        j = -epsilonbaux(2)
            i = i + 1
            goto 3001

        end if

3003    if ((num_dims == 3) .and. (k < epsilonbaux(3))) then
            k = k + 1
            i = -epsilonbaux(1); j = -epsilonbaux(2)
            goto 3001
        end if

    end subroutine s_applygaussian

    !> The purpose of this procedure is to smear the strength of the bubbles
        !!      in the Eulerian framework using a delta kernel function.
    subroutine s_deltafunc(updatedvar, cell, strength)

        real(kind(0.d0)) :: Dchar, strength, func, distance, Vol
        integer, dimension(3) :: cell
        integer :: idir, i, j
        type(scalar_field) :: updatedvar

        if (num_dims == 2) then
            Vol = dx(cell(1))*dy(cell(2))*charwidth
            if (cyl_coord) Vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2d0*PI
        else
            Vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
        end if

        updatedvar%sf(cell(1), cell(2), cell(3)) = updatedvar%sf(cell(1), cell(2), cell(3)) + strength/Vol

    end subroutine s_deltafunc

    !> The purpose of this procedure is to apply the delta kernel function.
        !! @param updatedvar Variable to be updated
        !! @param strength Variable to be smeared
        !! @param cell  Computational coordinates of the cell that contains the bubble.
        !! @param psi Local psi coordinates
    subroutine s_applydelta(updatedvar, strength, cell, psi)

        type(scalar_field) :: updatedvar
        real(kind(0.d0)), dimension(2, 2, 2) :: stencil !local coordinates
        real(kind(0.d0)), dimension(3) :: psi !local coordinates
        real(kind(0.d0)) :: strength
        integer, dimension(3) :: cell, cellaux
        integer :: i, j, k

        stencil(1, 1, 1) = (1.0d0 - psi(1))*(1.0d0 - psi(2))*(1.0d0 - psi(3))
        stencil(2, 1, 1) = psi(1)*(1.0d0 - psi(2))*(1.0d0 - psi(3))
        stencil(1, 2, 1) = (1.0d0 - psi(1))*psi(2)*(1.0d0 - psi(3))
        stencil(2, 2, 1) = psi(1)*psi(2)*(1.0d0 - psi(3))

        if (p > 0) then
            stencil(1, 1, 2) = (1.0d0 - psi(1))*(1.0d0 - psi(2))*psi(3)
            stencil(2, 1, 2) = psi(1)*(1.0d0 - psi(2))*psi(3)
            stencil(1, 2, 2) = (1.0d0 - psi(1))*psi(2)*psi(3)
            stencil(2, 2, 2) = psi(1)*psi(2)*psi(3)
            do i = 0, 1; do j = 0, 1; do k = 0, 1
                        cellaux(1) = min(max(cell(1) + i, -buff_size), m + buff_size)
                        cellaux(2) = min(max(cell(2) + j, -buff_size), n + buff_size)
                        cellaux(3) = min(max(cell(3) + k, -buff_size), p + buff_size)
                        updatedvar%sf(cellaux(1), cellaux(2), cellaux(3)) = updatedvar%sf(cellaux(1), cellaux(2), cellaux(3)) + stencil(i + 1, j + 1, k + 1)*strength
                    end do; end do; end do
        else
            do i = 0, 1; do j = 0, 1
                    cellaux(1) = min(max(cell(1) + i, -buff_size), m + buff_size)
                    cellaux(2) = min(max(cell(2) + j, -buff_size), n + buff_size)
                    updatedvar%sf(cellaux(1), cellaux(2), 0) = updatedvar%sf(cellaux(1), cellaux(2), 0) + stencil(i + 1, j + 1, 1)*strength
                end do; end do
        end if

    end subroutine s_applydelta

    subroutine s_remeshdelta(updatedvar, node, strength, rangecells)

        real(kind(0.d0)) :: strength, func, distance
        real(kind(0.d0)), dimension(3) :: psi
        integer, dimension(3) :: cell
        type(particledata) :: node
        type(scalar_field) :: updatedvar
        integer :: i
        integer, dimension(3, 2), optional :: rangecells

        call s_get_psi(node%tmp%s, psi, cell)
        call s_applydelta(updatedvar, strength, cell, psi)

        if (present(rangecells)) call s_update_rangecells(cell, rangecells)

    end subroutine s_remeshdelta

    subroutine s_update_rangecells(cell, rangecells)

        integer, dimension(3, 2), optional :: rangecells
        integer, dimension(3) :: cell
        integer :: i

        do i = 1, 3
            rangecells(i, 1) = max(min(rangecells(i, 1), cell(i) - 1), -buff_size)
            rangecells(i, 2) = max(rangecells(i, 2), cell(i) + 1)
        end do
        rangecells(1, 2) = min(rangecells(1, 2), m + buff_size)
        rangecells(2, 2) = min(rangecells(2, 2), n + buff_size)
        if (p == 0) then
            rangecells(3, :) = 0
        else
            rangecells(3, 2) = min(rangecells(3, 2), p + buff_size)
        end if

    end subroutine s_update_rangecells

    !> The purpose of this subroutine is to locate temporal spatial and computational coordinates
        !!      of the bubles when applying symmetric boundary condition.
        !! @param center Real coordinates of the lagragian bubble
        !! @param centertmp Temporal coordinates of the lagragian bubble for symmetric bc
        !! @param cell  Real computational coordinates of the cell that contains the bubble
        !! @param celltmp Temporal computational coordinates for symmetric bc
        !! @param i Coordinate index (i=1,2,3 -> x,y,z)
        !! @param j Boundary index (j=1,2 -> beg,end)
    subroutine s_kernel_shift_bc(center, centertmp, cell, celltmp, i, j)

        real(kind(0.d0)), dimension(3) :: center, centertmp
        real(kind(0.d0)) :: strength, stddsv
        integer, dimension(3) :: cell, celltmp
        integer :: i, j

        if (i == 1) then

            if (j == 1) then
                if (bc_x%beg == proc_rank) then
                    centertmp(1) = center(1) + x_cb(m) - x_cb(-1)
                    celltmp(1) = cell(1) + m + 1
                else
                    centertmp(1) = 2*x_cb(-1) - center(1)
                    celltmp(1) = -1 - cell(1)
                end if
            else
                if (bc_x%end == proc_rank) then
                    centertmp(1) = center(1) - x_cb(m) + x_cb(-1)
                    celltmp(1) = cell(1) - m - 1
                else
                    centertmp(1) = 2*x_cb(m) - center(1)
                    celltmp(1) = 2*m + 1 - cell(1)
                end if
            end if
        else if (i == 2) then
            if (j == 1) then
                if (bc_y%beg == proc_rank) then
                    centertmp(2) = center(2) + y_cb(n) - y_cb(-1)
                    celltmp(2) = cell(2) + n + 1
                else
                    centertmp(2) = 2*y_cb(-1) - center(2)
                    celltmp(2) = -1 - cell(2)
                end if
            else
                if (bc_y%end == proc_rank) then
                    centertmp(2) = center(2) - y_cb(n) + y_cb(-1)
                    celltmp(2) = cell(2) - n - 1
                else
                    centertmp(2) = 2*y_cb(n) - center(2)
                    celltmp(2) = 2*n + 1 - cell(2)
                end if
            end if
        else if (i == 3) then
            if (j == 1) then
                if (bc_z%beg == proc_rank) then
                    centertmp(3) = center(3) + z_cb(p) - z_cb(-1)
                    celltmp(3) = cell(3) + p + 1
                else
                    centertmp(3) = 2*z_cb(-1) - center(3)
                    celltmp(3) = -1 - cell(3)
                end if
            else
                if (bc_z%end == proc_rank) then
                    centertmp(3) = center(3) - z_cb(p) + z_cb(-1)
                    celltmp(3) = cell(3) - p - 1
                else
                    centertmp(3) = 2*z_cb(p) - center(3)
                    celltmp(3) = 2*p + 1 - cell(3)
                end if
            end if
        end if

    end subroutine s_kernel_shift_bc

    !> This subroutine returns the local psi coordinates, it also gives
            !!        the cell with respect these coordinates are defined.
            !! @param coord Input coordintes
            !! @param psi   Local psi coordinates
            !! @param cell  Computational coordinate of the cell
    subroutine s_get_psi(coord, psi, cell)

        real(kind(0.d0)), dimension(3), intent(in) :: coord
        real(kind(0.d0)), dimension(3) :: psi !local coordinates
        integer, dimension(3) :: cell

        call s_get_cell(coord, cell)

        !obtain psi(1)
        psi(1) = (coord(1) - real(cell(1)))*dx(cell(1)) + x_cb(cell(1) - 1)
        if (cell(1) == (m + buff_size)) then
            cell(1) = cell(1) - 1
            psi(1) = 1.0d0
        else if (cell(1) == (-buff_size)) then
            psi(1) = 0.0d0
        else
            if (psi(1) < x_cc(cell(1))) cell(1) = cell(1) - 1
            psi(1) = abs((psi(1) - x_cc(cell(1)))/(x_cc(cell(1) + 1) - x_cc(cell(1))))
        end if

        !obtain psi(2)
        psi(2) = (coord(2) - real(cell(2)))*dy(cell(2)) + y_cb(cell(2) - 1)
        if (cell(2) == (n + buff_size)) then
            cell(2) = cell(2) - 1
            psi(2) = 1.0d0
        else if (cell(2) == (-buff_size)) then
            psi(2) = 0.0d0
        else
            if (psi(2) < y_cc(cell(2))) cell(2) = cell(2) - 1
            psi(2) = abs((psi(2) - y_cc(cell(2)))/(y_cc(cell(2) + 1) - y_cc(cell(2))))
        end if

        !obtain psi(3)
        if (p > 0) then
            psi(3) = (coord(3) - real(cell(3)))*dz(cell(3)) + z_cb(cell(3) - 1)
            if (cell(3) == (p + buff_size)) then
                cell(3) = cell(3) - 1
                psi(3) = 1.0d0
            else if (cell(3) == (-buff_size)) then
                psi(3) = 0.0d0
            else
                if (psi(3) < z_cc(cell(3))) cell(3) = cell(3) - 1
                psi(3) = abs((psi(3) - z_cc(cell(3)))/(z_cc(cell(3) + 1) - z_cc(cell(3))))
            end if
        else
            psi(3) = 0.0d0
        end if

    end subroutine s_get_psi

    !> This subroutine finds the computational coordintates of the cell
            !!        where the given coordinates lies.
            !! @param s         Coordintes of the lagrangian particle
            !! @param get_cell  Cell in the computational coordinates
    subroutine s_get_cell(s, get_cell)

        integer, dimension(3) :: get_cell
        integer :: i
        real(kind(0.d0)), dimension(3) :: s

        get_cell(:) = int(s(:))

        do i = 1, num_dims
            if (s(i) < 0.0d0) get_cell(i) = get_cell(i) - 1
        end do

    end subroutine s_get_cell

end module m_kernel_functions
