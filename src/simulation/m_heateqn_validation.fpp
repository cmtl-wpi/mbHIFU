!>
!! @file m_heateqn.f90
!! @brief Contains module m_viscous

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to study HIFU
module m_heateqn

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_viscous

    use m_hifu

    ! ==========================================================================
    implicit none
    private; public :: s_cbc_heatEqn, s_bessel_I0, s_heatEqn_analyticalSol

contains

    ! ==========================================================================
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! READ ME: !!!!!!!!
    ! Important variables and where they are stored:
    ! q_cons_hifu(1)%sf(j,k,l): Heat intensity (generated). Summation over time (finite number of samples) OR
    !                           Error when running validation
    ! q_cons_hifu(2)%sf(j,k,l): Number of finite samples to compute q_cons_hifu(1).
    ! q_cons_hifu(3)%sf(j,k,l): Temperature.
    ! q_cons_hifu(4)%sf(j,k,l): RHS value from heat transfer eqn (finite volume discretization).
    ! q_cons_hifu(5)%sf(j,k,l): Analytical Solution for the specific problem.
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
    !Analytical solution to test the following problem ONLY
    ! # 2D heat difffusion problem. Clindircal coord.
    ! # The aim is to validate the heat transfer solver, no heat source, with 4 given temperature boundaries (see plot below)
    ! # Description:
    ! #
    ! #         ========= T1 = 2*T0 = 600 K =========
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #      T0 = 300 K                           T0 = 300 K
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #        ||                                    ||
    ! #         ==== T0 = finite | rotational axis ====
    ! #
    subroutine s_heatEqn_analyticalSol(q_cons_hifu)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_hifu
        integer :: j, k, h, i, imax !< Generic loop iterators
        real(kind(0d0)) :: T0, T1 !bc left, right, upper, bottom
        real(kind(0d0)) :: Cn, r0, L, bess_r0, bess_r, val1, val2
        integer :: val_idx

        val_idx = u_hifu_idx+1

        if (proc_rank==0) print*, 'Populating the analytical solution'

        !INPUT VALUES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        T0 = hifu_Tref
        T1 = 2*hifu_Tref
        r0 = 25.0
        L = 50.0
        imax = 400

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        h = 0
        do  j = 0, m  ! In x-dir
            do k = 0, n ! In y-dir	    
                val1 = 0.0
                do i = 1, imax
                    call s_bessel_I0(bess_r0,(i*pi*r0/L))
                    call s_bessel_I0(bess_r,(i*pi*y_cc(k)/L)) ! r = y_cc
                    Cn = (-2/(i*pi*bess_r0))*(cos(i*pi)-1)
                    val1 = val1 + Cn*bess_r*sin(i*pi*x_cc(j)/L)
                end do
                q_cons_hifu(val_idx)%sf(j,k,h) = T0 + (T1-T0) * val1

                !if (q_cons_hifu(5)%sf(j,k,h) /= 0) print*, 'Analytic solution not equal to zero'
                !if (proc_rank==0 .and. j==10 .and. k==10) print*, 'Analytic temperature probe is ', q_cons_hifu(5)%sf(j,k,h)
            end do
        end do
        
        call s_populate_HIFU_variables_buffers(q_cons_hifu)

    end subroutine s_heatEqn_analyticalSol ! =====================================


    ! ! ==========================================================================
    !Update temperature boundary condition (Intended only for validation of heatEqn solver with 2D diffusion rod problem)
    subroutine s_cbc_heatEqn(q_cons_hifu)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_hifu
        integer :: j, k, l !< Generic loop iterators
        real(kind(0d0)) :: T_L, T_R, T_U, T_B !bc left, right, upper, bottom
        integer :: val_idx

        val_idx = u_hifu_idx+1
        T_L = hifu_Tref
        T_R = hifu_Tref
        T_U = 2*hifu_Tref
        !T_B is symmetric
        
        l = 0
        do  j = 0, m            ! In x-dir
            do k = 0, n         ! In y-dir		

                ! Left boundary, 2*T0 at x=0 mm         
                if (abs(x_cc(j)) < 2*dx(j)/3) then                                     
                    !At the left-upper corner
                    if (abs(y_cc(k)-25.0) < 2*dy(k)/3) then
                        q_cons_hifu(val_idx)%sf(j,k,l) = (T_L+T_U)/2
                    else
                        q_cons_hifu(val_idx)%sf(j,k,l) = T_L
                    end if

                ! Right boundary, T0 at x=50 mm
                else if (abs(x_cc(j)-50.0) < 2*dx(j)/3) then                                      
                    !At the right-upper corner
                    if (abs(y_cc(k)-25.0) < 2*dy(k)/3) then
                        q_cons_hifu(val_idx)%sf(j,k,l) = (T_L+T_U)/2
                    else
                        q_cons_hifu(val_idx)%sf(j,k,l) = T_R
                    end if

                ! Upper boundary, T0 at y=25 mm
                else if (abs(y_cc(k)-25.0) < 2*dy(k)/3) then
                    q_cons_hifu(val_idx)%sf(j,k,l) = T_U

                ! Bottom boundary, T0 at y=0 mm (symmetric boundary)
                !else if (abs(y_cc(k)) < 2*dy(k)/3) then
                    !q_cons_hifu(3)%sf(j,k,l) = hifu_Tref

                end if 
            end do
        end do

        call s_populate_HIFU_variables_buffers(q_cons_hifu)

    end subroutine s_cbc_heatEqn ! =============================================

    ! ==========================================================================
    !Modified Bessel function for analytical solution I0(x)
    subroutine s_bessel_I0(yy, xx)

        real(kind(0d0)), intent(IN) :: xx
        real(kind(0d0)), intent(OUT) :: yy
        integer :: i, j, k !< Generic loop iterators
        real(kind(0d0)) :: var1, var2, var3

        yy   = 1
        var1 = 1
        var2 = 1

        do i = 1,10
            var1 = i*2.0;
            var2 = var2*(var1**2.0);
            yy = yy + (xx**var1)/var2;
        end do
    end subroutine s_bessel_I0 ! ===============================================

end module
