!>
!! @file m_hifu.f90
!! @brief Contains module m_hifu

#:include 'macros.fpp'

!> @brief HIFU post_processing
module m_hifu

    ! Dependencies ============================================================
    use m_derived_types      !< Definitions of the derived types
    use m_global_parameters  !< Definitions of the global parameters
    use m_mpi_proxy          !< Message passing interface (MPI) module proxy

    ! ==========================================================================

    implicit none

contains

    !> Initializes the hifu stages
    subroutine s_HIFU_start_stages()

        ! Post process only the most advanced hifu stage
        if (hifu_params%stg1 .and. hifu_params%stg2 .and. .not. hifu_params%stg3) then
            hifu_params%stg1 = .false.
        else if (hifu_params%stg3) then
            hifu_params%stg1 = .false.
            hifu_params%stg2 = .false.
        end if

        ! Set time params
        if (hifu_params%stg1) then
            if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 1 -> Post_process'
            hifu = .false.
            return
        else if (hifu_params%stg2) then
            if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 2 -> Post_process'
            if (cfl_dt) then
                if (n_start == 0) n_start = int(hifu_params%t_stop_stg1/t_save)
                t_stop = hifu_params%t_stop_stg2
            else
                if (mod(hifu_params%t_step_stop_stg1, t_step_save) == 0) then
                    t_step_start = hifu_params%t_step_stop_stg1
                else
                    t_step_start = (hifu_params%t_step_stop_stg1/t_step_save + 1)*t_step_save
                end if
                t_step_stop = hifu_params%t_step_stop_stg2
            end if
            if (bubbles_lagrange) lag_hifu_wrt = .true.
        else if (hifu_params%stg3) then
            if (proc_rank == 0) print *, 'WARNING :: HIFU -> Stage 3 -> Post_process'
            cfl_dt = .false.
            t_step_start = hifu_params%t_step_save_stg3
            t_step_save = hifu_params%t_step_save_stg3
            ! t_step_stop = hifu_params%t_step_stop_stg3 - 1
            t_step_stop = hifu_params%t_step_save_stg3

            ! if (hifu_params%stg3_3d) then
            !     if (hifu_params%cartesian) then
            !         m = hifu_params%m
            !         n = hifu_params%n
            !         p = hifu_params%p
            !         cyl_coord = .false.
            !         m_glb = m
            !         n_glb = n
            !         p_glb = p
            !         num_dims = 3
            !         grid_geometry = 1
            !         nGlobal = (m_glb + 1)*(n_glb + 1)*(p_glb + 1)

            !         bc_x%beg = -6; bc_x%end = -6
            !         bc_y%beg = -6; bc_y%end = -6
            !         bc_z%beg = -6; bc_z%end = -6
            !     else
            !         p = hifu_params%p
            !         p_glb = p
            !         nGlobal = (m_glb + 1)*(n_glb + 1)*(p_glb + 1)
            !         num_dims = 3
            !         if (bc_x%beg == -20) bc_x%beg = -6  ! from -20: acoustic bc
            !         bc_z%beg = -1; bc_z%end = -1  ! Assume entire cylindrical ring is taking care by one processor
            !         if (bc_y%beg == -2) bc_y%beg = -21  !   from -2: reflective boundary
            !     end if
            ! else
            !     if (p > 0) then  ! Full 3D
            !         ! if (bc_x%beg == BC_ROT_PERIODIC) bc_x%beg = -6 if (bc_y%beg == BC_ROT_PERIODIC) bc_y%beg = -6
            !         if (bc_z%beg == -20) bc_z%beg = -6
            !     end if
            ! end if
        end if

    end subroutine s_HIFU_start_stages

    !> Initializes the hifu stages
    subroutine s_HIFU_indexes()

        if (hifu_params%stg2) then
            hifu_params%qac_idx = 1
            hifu_params%qac_prms_idx = 2
            hifu_params%tsamp_idx = 3
            hifu_params%P_idx = 4
        end if

        if (hifu_params%stg3) then
            hifu_params%T_idx = 1
            hifu_params%qac_idx = 3
            hifu_params%qvis_idx = 4
            hifu_params%qth_idx = 5
        end if

    end subroutine s_HIFU_indexes

end module m_hifu
