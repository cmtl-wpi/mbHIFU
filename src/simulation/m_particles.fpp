!>
!! @file m_particles.f90
!! @brief Contains module m_particles

#:include 'macros.fpp'

!> @brief This module is used to add the lagrangian subgrid bubble model
module m_particles

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif
    use m_global_parameters     !< Definitions of the global parameters

    use m_derived_types         !< Definitions of the derived types

    use m_rhs                   !< Right-hand-side (RHS) evaluation procedures

    use m_data_output           !< Run-time info & solution data output procedures

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_kernel_functions      !< Definitions of the kernel functions

    ! ==========================================================================

    implicit none

    type(list3D), target            :: qbl
    type(cellListinfo), pointer     :: cellwbList
    type(bounds_info), dimension(3) :: subrange
    type(particleListinfo), pointer :: particlesubList

    real(kind(0.d0)) :: Rmax = 0.
    real(kind(0.d0)) :: Rmin = 1000000

    integer :: id

contains

    !> Initializes the lagrangian subgrid bubble solver
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_initialize_lagrangian_solver(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0.d0)) :: dtoutput, tend
        integer :: i, imax, save_count

        dt0=dt

        ! Initializing particle modules
        if (.not. avgdensFlag) solverapproach = 0

        if ((solverapproach == 2) .and. avgdensflag) then
            ! comp 1: (1 - beta)
            ! comp 2: dbetadt
            ! comp 3 - imax : auxiliary variables
            imax = 4
            if (clusterflag >= 4) imax = 10 !Subgrid noise model
            bubblesources = .true.
        else if ((solverapproach == 0) .and. avgdensflag) then
            imax = 3
            bubblesources = .false.
        else
            imax = 2
            bubblesources = .false.
        end if
        allocate (q_particle(1:imax))
        ! 1: one minus the voidfraction (1-beta)
        ! 2: Temporal derivative of the void fraction
        ! 3-5: Extra-allocated variables in those cases where source terms are required
        do i = 1, imax
            if (p > 0) then
                allocate (q_particle(i)%sf(-buff_size:m + buff_size, &
                                           -buff_size:n + buff_size, -buff_size:p + buff_size))
            else
                allocate (q_particle(i)%sf(-buff_size:m + buff_size, &
                                           -buff_size:n + buff_size, 0:0))
            end if
        end do
        q_particle(1)%sf = 1.d0 !represents 1-beta
        do i = 2, imax
            q_particle(i)%sf = 0.d0
        end do

        call s_read_input_particles(q_cons_vf, q_prim_vf)

        if (cfl_dt) then
            save_count = int(mytime/t_save)
        else
            save_count = t_step_start
        end if

        ! if (save_count == 0) then
            call s_write_data_files(q_cons_vf, q_prim_vf, save_count, q_particle(1)) !Parameter 'parallel_io' must be True
            call s_write_restart_particles(save_count) !parallel
            if (avgdensflag) call s_write_void_evol(mytime)
        ! end if
        call s_populate_variables_buffers(q_cons_vf, q_particle=q_particle)

    end subroutine s_initialize_lagrangian_solver

    !> The purpose of this procedure is to read the input file with the particles' information
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_read_input_particles(q_cons_vp, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vp
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf

        character(len=len_trim(case_dir) + 3*name_len) :: file_loc
        type(particlenode), pointer :: particle
        real(kind(0.d0)), dimension(8) :: inputparticle
        real(kind(0.d0)) :: qtime
        integer :: i, j, k, l, nparticles, save_count
        logical :: file_exist, indomain
        type(particlenode), pointer :: node
        type(cellwb), pointer :: cellwbaux

        ! To get the pressure to be used in s_add_particle()
        call s_populate_variables_buffers(q_cons_vp, q_particle=q_particle)
        call s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp%vf, ix, iy, iz, q_particle(1))

        ! Initialize particle lists and allocate list space
        nparticles = 0
        allocate (qbl%fp(-buff_size:m + buff_size, -buff_size:n + buff_size, -buff_size:p + buff_size))
        do k = -buff_size, p + buff_size
            do j = -buff_size, n + buff_size
                do i = -buff_size, m + buff_size
                    qbl%fp(i, j, k)%nb = 0
                    allocate (qbl%fp(i, j, k)%List)
                    nullify (qbl%fp(i, j, k)%List%next)
                    nullify (qbl%fp(i, j, k)%List%prev)
                    nullify (qbl%fp(i, j, k)%List%data)
                end do
            end do
        end do
        allocate (cellwbList)
        allocate (cellwbList%List)
        nullify (cellwbList%List%next)
        nullify (cellwbList%List%prev)
        nullify (cellwbList%List%data)
        cellwbList%nb = 0

        ! Read the input particle file or restart point
        if (cfl_dt) then
            save_count = n_start
        else
            save_count = t_step_start
        end if

        write (file_loc, '(a,i0,a)') 'particle_mpi_io', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=file_exist)

        if (save_count == 0 .or. (.not. file_exist)) then
            if (proc_rank==0) print*, 'Reading particles input file at save_count: ', save_count
            inquire (file='input/particles.dat', exist=file_exist)
            if (file_exist) then
                open (unit=85, file='input/particles.dat', form='formatted')
101             read (85, *, end=102) (inputparticle(i), i=1, 8)
                indomain = particle_in_domain(inputparticle(1:3), inputCoord = .true.)
                id = id + 1
                if (indomain) then
                    nparticles = nparticles + 1
                    call s_add_particle(inputparticle, q_cons_vp, q_prim_vf)
                end if
                goto 101
102             continue
            else
                stop "if you include particles, you have to initializate them in input/particles.dat"
            end if
        else
            if (proc_rank==0) print*, 'Restarting particles at save_count: ', save_count
            call s_add_particle_restart(nparticles, save_count) !Parameter 'parallel_io' must be True
        end if

        print *, " PARTICLES RUNNING, in proc", proc_rank, "number:", nparticles, "/", id

        ! Allocate needed sublists space
        allocate (particlesubList)
        allocate (particlesubList%List)
        nullify (particlesubList%next)
        nullify (particlesubList%prev)
        nullify (particlesubList%List%next)
        nullify (particlesubList%List%prev)
        nullify (particlesubList%List%data)
        particlesubList%nb = 0
        cellwbaux => cellwbList%List%next
        do while (associated(cellwbaux))
            node => qbl%fp(cellwbaux%data%coord(1), cellwbaux%data%coord(2), cellwbaux%data%coord(3))%List%next
            do while (associated(node))
                call s_add_particle_to_list(node%data, particlesubList)
                node => node%next
            end do
            cellwbaux => cellwbaux%next
        end do

        ! Apply density correction
        if (avgdensFlag) then
            particle => particlesubList%List%next
            do while (associated(particle))
                call s_transfer_data_to_tmp(particle%data)
                particle%data%tmp%shell = particle%data%shell
                particle => particle%next
            end do
            call s_smear_voidfraction(q_cons_vp)
            if (solverapproach == 1) then
                !Definition averaged quantities
                q_cons_vp(E_idx)%sf(0:m, 0:n, 0:p) = q_cons_vp(E_idx)%sf(0:m, 0:n, 0:p)*q_particle(1)%sf(0:m, 0:n, 0:p)
                do i = 1, cont_idx%end
                    q_cons_vp(i)%sf(0:m, 0:n, 0:p) = q_cons_vp(i)%sf(0:m, 0:n, 0:p)*q_particle(1)%sf(0:m, 0:n, 0:p)
                end do
                do i = mom_idx%beg, mom_idx%end
                    q_cons_vp(i)%sf(0:m, 0:n, 0:p) = q_cons_vp(i)%sf(0:m, 0:n, 0:p)*q_particle(1)%sf(0:m, 0:n, 0:p)
                end do
            end if
        end if

        qtime = 0.0d0
        if (cfl_dt) then
            qtime = n_start*t_save
        else
            qtime = t_step_start*dt
        end if
        mytime = qtime
        
        if (particleoutFlag) call s_write_particles(qtime)
        call s_populate_variables_buffers(q_cons_vp, q_particle=q_particle)

        call s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp%vf, &
                                                           ix, iy, iz, q_particle(1))

    end subroutine s_read_input_particles

    !> The purpose of this procedure is to add information of the particles when starting fresh
        !! @param inputparticle Particle number
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_add_particle(inputparticle, q_cons_vp, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vp
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0.d0)), dimension(8), intent(IN) :: inputparticle
        type(cellwb), pointer :: cellnode
        type(cellwbcoord), pointer :: cellinfo
        type(particledata), pointer :: particleinfo
        type(particleListinfo), pointer :: particleListaux
        real(kind(0.d0)) :: pliq, volparticle, concvap, totalmass, kparticle, cpparticle
        real(kind(0.d0)) :: omegaN, PeG, PeT, rhol, cson, pcrit, qv, gamma, pi_inf
        integer, dimension(3) :: cell
        real(kind(0.d0)), dimension(2) :: Re
        real(kind(0.d0)) :: sigmaTotal, bubbleTemp, act_massflag

        allocate (particleinfo)
        particleinfo%id = id
        particleinfo%x(:) = inputparticle(1:3)
        particleinfo%xprev(:) = inputparticle(1:3)
        particleinfo%u(:) = inputparticle(4:6)
        particleinfo%y(1) = inputparticle(7)  ! particle radius
        particleinfo%R0 = inputparticle(7)
        particleinfo%y(2) = inputparticle(8)  ! interface velocity
        particleinfo%Rmax = 1.0d0
        particleinfo%Rmin = 1.0d0
        particleinfo%dphidt = 0.0d0

        if (cyl_coord .and. p == 0) then
            particleinfo%x(2) = dsqrt(particleinfo%x(2)**2d0 + particleinfo%x(3)**2d0)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            particleinfo%x(3) = atan2(inputparticle(3), inputparticle(2))
            particleinfo%xprev = particleinfo%x
        end if
        cell = -buff_size
        call s_locate_cell(particleinfo%x, cell, particleinfo%tmp%s, particleinfo%x(2))
        if (Pbase_bc == dflt_real) then
            pliq = f_interpolate(particleinfo%tmp%s, q_prim_vf(E_idx))
        else
            pliq = Pbase_bc
        end if
        
        if (pliq < 0) print *, "Negative pressure", proc_rank, &
            q_cons_vp(E_idx)%sf(cell(1), cell(2), cell(3)), q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)), cell
        call s_convert_to_mixture_variables(q_cons_vp, cell(1), cell(2), cell(3), rhol, gamma, pi_inf, qv, Re)

        ! Marmotant model parameters
        particleinfo%shell = lipidCoatingModel
        particleinfo%Rbuck = dble(particleinfo%shell) * (particleinfo%R0/sqrt(1.0d0+(sigma0_lipidCoat/surfaceElast_lipidCoat)))
        particleinfo%Rrupt = particleinfo%Rbuck * sqrt(1.0d0+(sigmabubble/surfaceElast_lipidCoat))
        act_massflag = max(dble(massflag) - dble(particleinfo%shell), 0d0)

        if (particleinfo%shell .eq. 1) then !Marmotant model
            sigmaTotal = sigma0_lipidCoat
        else
            sigmaTotal = sigmabubble
        end if

        ! Intial particle pressure
        particleinfo%p = pliq + 2.0d0*sigmaTotal/particleinfo%R0
        if (sigmaTotal /= 0.0d0) then
            pcrit = pvap - 4.0d0*sigmaTotal/(3.d0*sqrt(3.0d0*particleinfo%p*particleinfo%R0**3/(2.0d0*sigmaTotal)))
            pref = particleinfo%p
        else
            pcrit = 0.0d0
        end if

        particleinfo%equilibrium = .false.

        ! Initial particle mass
        volparticle = 4.0d0/3.0d0*pi*particleinfo%R0**3d0 ! volume
        particleinfo%mv = pvap*volparticle*(1.0d0/(Rvap*Thost))*act_massflag ! vapermass if no shell
        particleinfo%mg = (particleinfo%p - pvap*act_massflag)*volparticle*(1.0d0/(Rgas*Thost)) ! gasmass
        if (particleinfo%mg <= 0.0d0) stop 'the initial mass of gas inside the bubble is negative. Check your initial conditions'
        totalmass = particleinfo%mg + particleinfo%mv ! totalmass

        ! Bubble natural frequency
        concvap = particleinfo%mv/(particleinfo%mv + particleinfo%mg)
        omegaN = (3.0d0*(particleinfo%p - pvap*act_massflag) + 4.0d0*sigmaTotal/particleinfo%R0)/rhol
        if (pvap*act_massflag > particleinfo%p) then
            print *, 'Not allowed: bubble initially located in a region with pressure below the vapor pressure'
            print *, 'location:', particleinfo%x(1:3)
            stop
        end if
        omegaN = dsqrt(omegaN/particleinfo%R0**2)

        cpparticle = concvap*cpvapor + (1.0d0 - concvap)*cpgas
        kparticle = concvap*kvapor + (1.0d0 - concvap)*kgas

        ! Mass and heat transfer coefficients (based on Preston 2007)
        PeT = totalmass/volparticle*cpparticle*particleinfo%R0**2d0*omegaN/kparticle
        particleinfo%betaT = f_transfercoeff(PeT, 1.0d0)*dble(heatflag)
        PeG = particleinfo%R0**2d0*omegaN/diffcoefvap
        particleinfo%betaC = f_transfercoeff(PeG, 1.0d0)*act_massflag

        ! Terms to work out directly the heat flux in getfluxes
        particleinfo%betaT = particleinfo%betaT*kparticle

        if (particleinfo%mg <= 0.0d0) stop 'PROBLEM WITH THE MASS OF THE particle, CHECK PARTICLES ARE INSIDE THE doMAIN'

        ! volparticle = 4.0d0/3.0d0*pi*particleinfo%y(1)**3
        ! bubbleTemp = (particleinfo%p-pvap)*volparticle/(particleinfo%mg*Rgas)
        ! print '(A,7E24.17)', 'ADDparticle temperature', particleinfo%p, sigmaTotal

        particleListaux => qbl%fp(cell(1), cell(2), cell(3))
        if (particleListaux%nb == 0) then
            !Adding bubble to the head of the list
            allocate (cellnode)
            allocate (cellinfo)
            cellwbList%nb = cellwbList%nb + 1
            cellinfo%coord = cell
            cellnode%data => cellinfo
            cellnode%next => cellwbList%List%next
            nullify (cellnode%prev)
            qbl%fp(cell(1), cell(2), cell(3))%cellpointer => cellnode
            if (associated(cellwbList%List%next)) cellwbList%List%next%prev => cellnode
            cellwbList%List%next => cellnode
        end if
        call s_transfer_data_to_tmp(particleinfo)
        particleinfo%tmp%shell = particleinfo%shell
        call s_add_particle_to_list(particleinfo, particleListaux)

    end subroutine s_add_particle

    !> The purpose of this procedure is to add information of the particles
        !!      from a restart point in parallel
        !! @param nparticles Particle number in the domain
    subroutine s_add_particle_restart(nparticles, save_count)

        integer :: nparticles

        character(len=len_trim(case_dir) + 2*name_len) :: t_step_dir
        character(len=len_trim(case_dir) + 3*name_len) :: file_loc
        logical :: dir_check
#ifdef MFC_MPI
        real(kind(0.d0)), dimension(23) :: inputvals
        real(kind(0.d0)) :: id_real
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer :: view

        type(particledata), pointer :: particleinfo
        type(particlelistinfo), pointer :: particlelistaux

        type(cellwb), pointer :: cellnode
        type(cellwbcoord), pointer :: cellinfo

        integer, dimension(3) :: cell
        logical :: indomain, particle_file, file_exist

        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer :: ifile, ireq, ierr, data_size, tot_data
        integer :: i

        integer :: save_count
        integer :: varsMarmotant = 3

        write (file_loc, '(a,i0,a)') 'particle_mpi_io', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=file_exist)

        if (file_exist) then
            if (proc_rank == 0) then
                open (9, file=trim(file_loc), form='unformatted', status='unknown')
                !read (9) tot_data, time_real, dt_next_inp, tot_step
                read (9) tot_data, mytime, dt
                close (9)
                print*, 'Reading particle_mpi_io: ', tot_data, mytime, dt
            end if
        else
            print '(a)', trim(file_loc)//' is missing. exiting ...'
            call s_mpi_abort
        end if

        call MPI_BCAST(tot_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mytime, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        gsizes(1) = tot_data
        gsizes(2) = 21+varsMarmotant
        lsizes(1) = tot_data
        lsizes(2) = 21+varsMarmotant
        start_idx_part(1) = 0
        start_idx_part(2) = 0

        call MPI_type_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, MPI_doUBLE_PRECISION, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        ! Open the file to write all flow variables
        write (file_loc, '(a,i0,a)') 'particle', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=particle_file)

        if (particle_file) then
            call MPI_FILE_open(MPI_COMM_WORLD, file_loc, MPI_MODE_RdoNLY, &
                               mpi_info_int, ifile, ierr)
            disp = 0d0
            call MPI_FILE_SET_VIEW(ifile, disp, MPI_doUBLE_PRECISION, view, &
                                   'native', mpi_info_null, ierr)
            allocate (MPI_IO_DATA_particle(tot_data, 1:(21+varsMarmotant)))
            call MPI_FILE_read_ALL(ifile, MPI_IO_DATA_particle, (21+varsMarmotant)*tot_data, &
                                   MPI_doUBLE_PRECISION, status, ierr)
            do i = 1, tot_data
                id = int(MPI_IO_DATA_particle(i, 1))
                inputvals(1:(20+varsMarmotant)) = MPI_IO_DATA_particle(i, 2:(21+varsMarmotant))
                indomain = particle_in_domain(inputvals(1:3))
                if (indomain .and. (id > 0)) then
                    allocate (particleinfo)
                    nparticles = nparticles + 1
                    particleinfo%id = id
                    particleinfo%x(1:3) = inputvals(1:3)
                    particleinfo%xprev(1:3) = inputvals(4:6)
                    particleinfo%u(1:3) = inputvals(7:9)
                    particleinfo%y(1:2) = inputvals(10:11)
                    particleinfo%R0 = inputvals(12)
                    particleinfo%Rmax = inputvals(13)
                    particleinfo%Rmin = inputvals(14)
                    particleinfo%dphidt = inputvals(15)
                    particleinfo%p = inputvals(16)
                    particleinfo%mv = inputvals(17)
                    particleinfo%mg = inputvals(18)
                    particleinfo%betaT = inputvals(19)
                    particleinfo%betaC = inputvals(20)
                    ! varsMarmotant
                    particleinfo%shell = inputvals(21)
                    particleinfo%Rbuck = inputvals(22)
                    particleinfo%Rrupt = inputvals(23)

                    particleinfo%equilibrium = .false.
                    cell = -buff_size
                    call s_locate_cell(particleinfo%x, cell, particleinfo%tmp%s)
                    particleListaux => qbl%fp(cell(1), cell(2), cell(3))
                    if (particleListaux%nb == 0) then
                        !Adding bubble to the head of the list
                        allocate (cellnode); allocate (cellinfo)
                        cellwbList%nb = cellwbList%nb + 1
                        cellinfo%coord = cell
                        cellnode%data => cellinfo
                        cellnode%next => cellwbList%List%next
                        nullify (cellnode%prev)
                        qbl%fp(cell(1), cell(2), cell(3))%cellpointer => cellnode
                        if (associated(cellwbList%List%next)) cellwbList%List%next%prev => cellnode
                        cellwbList%List%next => cellnode
                    end if
                    call s_transfer_data_to_tmp(particleinfo)
                    particleinfo%tmp%shell = particleinfo%shell
                    call s_add_particle_to_list(particleinfo, particleListaux)
                end if
            end do
            deallocate (MPI_IO_DATA_particle)
        end if
        call MPI_FILE_CLOSE(ifile, ierr)
#endif

    end subroutine s_add_particle_restart

    !> The purpose of this subroutine is to add the particles information
          !!        in the allocated list space
          !! @param particle     Particle data
          !! @param particleList Particle list
    subroutine s_add_particle_to_list(particle, particleList)

        !Adding particle to the head of the list
        type(particledata), pointer :: particle
        type(particlenode), pointer :: node
        type(particleListinfo), pointer :: particleList

        allocate (node)
        particleList%nb = particleList%nb + 1
        node%data => particle
        node%next => particleList%List%next
        nullify (node%prev)
        if (associated(particleList%List%next)) particleList%List%next%prev => node
        particleList%List%next => node

    end subroutine s_add_particle_to_list

    !>  The purpose of this subroutine is to smear the particles in the Eulerian frame
        !! @param q Smeared particles in the conservative variables
    subroutine s_smear_voidfraction(q)

        type(scalar_field), dimension(sys_size), intent(IN) :: q !remove this q

        type(particlenode), pointer :: particle
        real(kind(0.d0)) :: volpart, totmass, stddsv, volpart2
        real(kind(0.d0)), dimension(3) :: nodecoord
        integer, dimension(3) :: cell
        integer :: i, j, k, l, kernel
        integer, dimension(3, 2) :: rangecells

        q_particle(1)%sf = 0.0d0
        q_particle(2)%sf = 0.0d0
        nodecoord(3) = 0

        rangecells(1:3, 2) = -buff_size
        rangecells(1, 1) = m
        rangecells(2, 1) = n
        rangecells(3, 1) = p

        if (projectiontype == 0) then
            particle => particlesubList%List%next
            do while (associated(particle))
                volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
                call s_get_cell(particle%data%tmp%s, cell)
                call s_update_rangecells(cell, rangecells)
                nodecoord(1) = particle%data%tmp%x(1)
                nodecoord(2) = particle%data%tmp%x(2)
                if (p > 0) nodecoord(3) = particle%data%tmp%x(3)
                call s_compute_stddsv(cell, kernel, volpart, stddsv)
                call s_smoothfunction(q_particle(1), nodecoord, cell, volpart, kernel, stddsv)
                particle => particle%next
            end do
            subrange(1)%beg = max(rangecells(1, 1), -buff_size)
            subrange(1)%end = min(rangecells(1, 2), m + buff_size)
            subrange(2)%beg = max(rangecells(2, 1), -buff_size)
            subrange(2)%end = min(rangecells(2, 2), n + buff_size)
            subrange(3)%beg = max(rangecells(3, 1), -buff_size)
            subrange(3)%end = min(rangecells(3, 2), p + buff_size)
        else
            particle => particlesubList%List%next
            do while (associated(particle))
                volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
                call s_remeshdelta(q_particle(2), particle%data, volpart, rangecells)
                particle => particle%next
            end do
            if (bubblesources) q_particle(4)%sf = q_particle(2)%sf
            subrange(1)%beg = max(rangecells(1, 1), -buff_size)
            subrange(1)%end = min(rangecells(1, 2), m + buff_size)
            subrange(2)%beg = max(rangecells(2, 1), -buff_size)
            subrange(2)%end = min(rangecells(2, 2), n + buff_size)
            subrange(3)%beg = max(rangecells(3, 1), -buff_size)
            subrange(3)%end = min(rangecells(3, 2), p + buff_size)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        cell(1) = i
                        cell(2) = j
                        cell(3) = k
                        nodecoord(1) = x_cc(cell(1))
                        nodecoord(2) = y_cc(cell(2))
                        if (p > 0) nodecoord(3) = z_cc(cell(3))

                        if (q_particle(2)%sf(i, j, k) /= 0.0d0) then
                            call s_compute_stddsv(cell, kernel, q_particle(2)%sf(i, j, k), stddsv)
                            call s_smoothfunction(q_particle(1), nodecoord, cell, q_particle(2)%sf(i, j, k), kernel, stddsv)
                        end if
                    end do
                end do
            end do

        end if

        !I store 1-beta (I should probably do it initializing
        !                q_particle(1)%fp=1 and using a negative volume)
        q_particle(1)%sf = 1.-q_particle(1)%sf
        subrange(1)%beg = max(subrange(1)%beg - ceiling(5*epsilonb), -buff_size*1.0d0)
        subrange(1)%end = min(subrange(1)%end + ceiling(5*epsilonb), (m + buff_size)*1.0d0)
        subrange(2)%beg = max(subrange(2)%beg - ceiling(5*epsilonb), -buff_size*1.0d0)
        subrange(2)%end = min(subrange(2)%end + ceiling(5*epsilonb), (n + buff_size)*1.0d0)
        if (p > 3) then
            subrange(3)%beg = max(subrange(3)%beg - ceiling(5*epsilonb), -buff_size*1.0d0)
            subrange(3)%end = min(subrange(3)%end + ceiling(5*epsilonb), (p + buff_size)*1.0d0)
        end if

        !=== dbetadt==========
        if (bubblesources) then
            q_particle(2)%sf = 0.0d0
            q_particle(3)%sf = 0.0d0
            if ((clusterflag >= 4)) then
                q_particle(5)%sf = 0.0d0
            end if
            if (projectiontype == 0) then
                particle => particlesubList%List%next
                do while (associated(particle))
                    volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
                    volpart2 = volpart
                    call s_get_cell(particle%data%tmp%s, cell)
                    call s_compute_stddsv(cell, kernel, volpart, stddsv)
                    volpart = 4.0d0*pi*particle%data%tmp%y(1)**2*particle%data%tmp%y(2)
                    nodecoord(1) = particle%data%tmp%x(1)
                    nodecoord(2) = particle%data%tmp%x(2)
                    if (p > 0) nodecoord(3) = particle%data%tmp%x(3)
                    if (clusterflag >= 4) call s_smoothfunction(q_particle(5), nodecoord, cell, volpart, kernel, stddsv, volpart2)
                    call s_smoothfunction(q_particle(2), nodecoord, cell, volpart, kernel, stddsv)
                    particle => particle%next
                end do
            else
                nodecoord(3) = 0
                particle => particlesubList%List%next
                do while (associated(particle))
                    volpart = 4.0d0*pi*particle%data%tmp%y(1)**2*particle%data%tmp%y(2)
                    call s_remeshdelta(q_particle(3), particle%data, volpart)
                    particle => particle%next
                end do
                do i = 0, m
                    do j = 0, n
                        do k = 0, p
                            cell(1) = i
                            cell(2) = j
                            cell(3) = k
                            nodecoord(1) = x_cc(cell(1))
                            nodecoord(2) = y_cc(cell(2))
                            if (p > 0) nodecoord(3) = z_cc(cell(3))
                            if (q_particle(3)%sf(i, j, k) /= 0.0d0) then
                                call s_compute_stddsv(cell, kernel, q_particle(4)%sf(i, j, k), stddsv)
                                call s_smoothfunction(q_particle(2), nodecoord, cell, q_particle(3)%sf(i, j, k), kernel, stddsv)
                            end if
                        end do
                    end do
                end do
            end if
        end if
        subrange(1)%beg = max(subrange(1)%beg, 0*1.0d0)
        subrange(1)%end = min(subrange(1)%end, m*1.0d0)
        subrange(2)%beg = max(subrange(2)%beg, 0*1.0d0)
        subrange(2)%end = min(subrange(2)%end, n*1.0d0)
        if (p > 0) then
            subrange(3)%beg = max(subrange(3)%beg, 0*1.0d0)
            subrange(3)%end = min(subrange(3)%end, p*1.0d0)
        end if

        ! Limiting void fraction given max value
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    q_particle(1)%sf(i, j, k) = max(q_particle(1)%sf(i, j, k), 1.d0 - valmaxvoid)
                end do
            end do
        end do

    end subroutine s_smear_voidfraction

    !> Calculates the standard deviation of the particle being smeared in the Eulerian frame
        !! @param cell Cell where the particle is located
        !! @param kernel Kernel type or smoothening function
        !! @param volpart Volume of the particle
        !! @param stddsv Standard deviaton
    subroutine s_compute_stddsv(cell, kernel, volpart, stddsv)

        integer :: kernel
        integer, dimension(3) :: cell
        real(kind(0.d0)) :: volpart, stddsv

        real(kind(0.d0)) :: chardist, rad, charvol

        call s_get_char_dist(cell, chardist)
        call s_get_char_vol(cell, charvol)
        kernel = smoothtype

        if (((volpart/charvol) > 0.5d0*valmaxvoid) .or. (smoothtype == 1)) then
            kernel = 1
            rad = (3.0d0*volpart/(4.0d0*pi))**(1.0d0/3.0d0)
            stddsv = 1.0d0*epsilonb*max(chardist, rad)
        else
            stddsv = 0.0d0
        end if

    end subroutine s_compute_stddsv

    !>  The purpose of this subroutine is to add the particle source terms following the formulation of Kazuki and Colonius (2018)
        !! @param q Conservative variables
        !! @param dq Calculated change of conservative variables
        !! @param q_prim Conservative variables
    subroutine s_add_sources(q, dq, q_prim)

        type(scalar_field), dimension(sys_size), intent(IN) :: q
        type(scalar_field), dimension(sys_size), intent(IN) :: dq
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim

        integer :: i, j, k, l

        real(kind(0d0)) :: qmin, qmax

        qmin = 100.0d0
        qmax = -100.0d0

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    if (q_particle(1)%sf(i, j, k) > (1.0d0 - valmaxvoid)) then
                        do l = 1, E_idx
                            if (clusterflag >= 4) then
                                dq(l)%sf(i, j, k) = dq(l)%sf(i, j, k) + q(l)%sf(i, j, k)*(q_particle(2)%sf(i, j, k) + q_particle(5)%sf(i, j, k)) ! Add the residual (g_res) to account the non-uniform bubble cloud (2D approximation)
                            else
                                dq(l)%sf(i, j, k) = dq(l)%sf(i, j, k) + q(l)%sf(i, j, k)/q_particle(1)%sf(i, j, k)*q_particle(2)%sf(i, j, k)
                            end if
                        end do
                    else
                        print*, 'correction applied', q_particle(1)%sf(i, j, k)
                    end if
                    qmin = min(qmin, q_particle(1)%sf(i, j, k))
                    qmax = max(qmax, q_particle(1)%sf(i, j, k))
                end do
            end do
        end do

        !print*, '1-beta, maximum and minimum', qmin, qmax, proc_rank

        do l = 1, num_dims
            call s_gradient_dir(q_prim(E_idx), q_particle(3), l)

            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (q_particle(1)%sf(i, j, k) > (1.0d0 - valmaxvoid)) then
                            dq(mom_idx%beg + l - 1)%sf(i, j, k) = dq(mom_idx%beg + l - 1)%sf(i, j, k) - (1.0d0 - q_particle(1)%sf(i, j, k))/q_particle(1)%sf(i, j, k)*q_particle(3)%sf(i, j, k)
                        end if
                    end do
                end do
            end do

            !source in energy
            q_particle(3)%sf = q_prim(E_idx)%sf*q_prim(mom_idx%beg + l - 1)%sf
            call s_gradient_dir(q_particle(3), q_particle(4), l)

            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (q_particle(1)%sf(i, j, k) > (1.0d0 - valmaxvoid)) then
                            dq(E_idx)%sf(i, j, k) = dq(E_idx)%sf(i, j, k) - q_particle(4)%sf(i, j, k)*(1.0d0 - q_particle(1)%sf(i, j, k))/q_particle(1)%sf(i, j, k)
                        end if
                    end do
                end do
            end do
        end do

    end subroutine s_add_sources

    !>  Contains the bubble dynamics subroutines
        !! @param qtime Current time from the adaptative 4th/5th order Runge-Kutta-Cash-Karp time stepper
        !! @param step Current time step in the adaptative stepper
        !! @param q Conservative variables
        !! @param t_step Current global time step
        !! @param q_prim Primitive variables
        !! @param dq Calculated change of conservative variables
    subroutine s_RK_particle_dynamics(qtime, step, q, t_step, q_prim, dq)

        type(scalar_field), dimension(sys_size) :: q, q_prim
        type(scalar_field), dimension(sys_size), optional :: dq
        integer, intent(IN) :: step, t_step
        real(kind(0.d0)) :: qtime

        type(particlenode), pointer :: particle
        real(kind(0.d0)), dimension(5) :: intvalues
        real(kind(0.d0)), dimension(3) :: totalforce, DupDt
        real(kind(0.d0)) :: gammaparticle, vaporflux, heatflux, time_avg
        integer, dimension(3) :: cell
        integer :: i, j, k, l
        type(cellwb), pointer :: cellwbaux
        real(kind(0.d0)) :: preterm1, term2, paux, pint, Romega, term1_fac, Rb
        real(kind(0.d0)), dimension(3) :: scoord

        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size
        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg

        if (avgdensFlag) then
            do i = 1, sys_size
                dq(i)%sf = 0.0d0
            end do
        end if

        !Update eulerian framework
        if (avgdensFlag) call s_smear_voidfraction(q)
        ! do i = 1, cont_idx%end
        !     do l = iz%beg, iz%end
        !         do k = iy%beg, iy%end
        !             do j = ix%beg, ix%end
        !                 q_prim(i)%sf(j, k, l) = q(i)%sf(j, k, l)
        !             end do
        !         end do
        !     end do
        ! end do
        ! do i = adv_idx%beg, sys_size
        !     do l = iz%beg, iz%end
        !         do k = iy%beg, iy%end
        !             do j = ix%beg, ix%end
        !                 q_prim(i)%sf(j, k, l) = q(i)%sf(j, k, l)
        !             end do
        !         end do
        !     end do
        ! end do

        if (coupledFlag) then
            call s_compute_rhs(q, q_prim, dq, t_step=t_step, time_avg=time_avg, qtime=qtime)
        end if

        if (avgdensFlag) then
            call s_populate_variables_buffers(q, q_particle=q_particle)
            ! ix%beg = -buff_size; ix%end = m + buff_size
            ! iy%beg = -buff_size; iy%end = n + buff_size
            ! if (p > 0) iz%beg = -buff_size; iz%end = p + buff_size
            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
            if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size
            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            call s_convert_conservative_to_primitive_variables(q, q_prim, gm_alpha_qp%vf, ix, iy, iz, q_particle(1))

            ! For conservative to primitive below
            do i = 1, sys_size
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            q_cons_qp%vf(i)%sf(j, k, l) = q(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            if ((clusterflag > 0) .or. (correctpresFlag)) then ! subgrid p_inf model from kazuki and Colonius

                ! Calculate potentials
                cellwbaux => cellwbList%List%next
                do while (associated(cellwbaux))
                    cell = cellwbaux%data%coord

                    ! Solve cell
                    scoord(:) = cell(:) + 0.5d0
                    paux = f_pressure_inf(scoord, q_prim(E_idx), 2, preterm1, term2, Romega)
                    particle => qbl%fp(cell(1), cell(2), cell(3))%List%next
                    do while (associated(particle))
                        call s_compute_pressureliq_int(pint, particle%data%tmp%p, particle%data%tmp%y(1), particle%data%tmp%y(2), &
                                                       particle%data%tmp%shell, particle%data%Rbuck, particle%data%Rrupt)
                        !pint = pint + 0.5d0*particle%data%tmp%y(2)**2

                        if (clusterflag == 2) then
                            !particle%data%dphidt = (paux - pint) + term2
                            particle%data%dphidt = paux - pint - (term2 - 0.5d0*particle%data%tmp%y(2)**2) ! need ( ) * density
                            ! Accouting for the potential induced by the bubble averaged over the control volume
                            ! Note that this is based on the incompressible flow assumption near the bubble.
                            Rb = particle%data%tmp%y(1)
                            term1_fac = 3.0d0/2.0d0*(Rb*(Romega**2d0 - Rb**2d0))/(Romega**3d0 - Rb**3d0)
                            particle%data%dphidt = particle%data%dphidt/(1.0d0 - term1_fac)
                        end if
                        particle => particle%next
                    end do
                    cellwbaux => cellwbaux%next
                end do
            end if
        end if

        if (bubblesources) call s_add_sources(q, dq, q_prim)

        particle => particlesubList%List%next
        do while (associated(particle))
            if (.not. particle%data%equilibrium) then
                call s_compute_interface_fluxes(particle%data, vaporflux, heatflux, gammaparticle)
                particle%data%dbdt(step)%dpbdt = -3.0d0*gammaparticle/particle%data%tmp%y(1)* &
                                                  (particle%data%tmp%p*particle%data%tmp%y(2) - &
                                                   heatflux - (Rvap*Thost)*vaporflux)
                particle%data%dbdt(step)%dmvdt = 4.0d0*pi*particle%data%tmp%y(1)**2*vaporflux
                !print*, 'heat and vapor fluxes: ', heatflux, vaporflux
            else
                particle%data%dbdt(step)%dpbdt = 0.0d0
                particle%data%dbdt(step)%dmvdt = 0.0d0
            end if
            particle%data%dbdt(step)%dxdt(:) = 0.0d0
            particle%data%dbdt(step)%dudt(:) = 0.0d0
            particle => particle%next
        end do

        !Radial motion model
        if (RPflag) then
            particle => particlesubList%List%next
            do while (associated(particle))
                if (.not. particle%data%equilibrium) then
                    call s_compute_RP(particle%data, step, q, q_prim, qtime)
                else
                    particle%data%dbdt(step)%dydt(2) = 0.
                end if
                particle%data%dbdt(step)%dydt(1) = particle%data%tmp%y(2)
                particle => particle%next
            end do
        else
            particle => particlesubList%List%next
            do while (associated(particle))
                particle%data%dbdt(step)%dydt(2) = 0.0d0
                particle%data%dbdt(step)%dydt(1) = 0.0d0
                particle => particle%next
            end do
        end if

    end subroutine s_RK_particle_dynamics

    !>  This subroutine computes the Keller-Miksis equation
        !! @param bubbletmp Current particle information
        !! @param step Current time step in the adaptative stepper
        !! @param q Conservative variables
        !! @param q_prim Primitive variables
        !! @param qtime Current time from the adaptative 4th/5th order Runge-Kutta-Cash-Karp time stepper
    subroutine s_compute_RP(bubbletmp, step, q, q_prim, qtime)

        type(scalar_field), dimension(sys_size) :: q
        type(scalar_field), dimension(sys_size) :: q_prim
        type(particledata) :: bubbletmp
        type(particlederivative) :: dbdt
        real(kind(0.d0)) :: pliqint, pbubble, deltaP, pinf, termI, aux1, &
                            aux2, velint, rhol, cson, E, H, qv, gamma, pi_inf
        integer, dimension(3) :: cell
        integer :: step, i
        real(kind(0.d0)), dimension(E_idx - mom_idx%beg) :: vel
        real(kind(0.d0)), dimension(2) :: Re
        real(kind(0.d0)), optional :: qtime
        real(kind(0.d0)) :: temp

        pbubble = bubbletmp%tmp%p ! pres in the bubble
        call s_compute_pressureliq_int(pliqint, pbubble, bubbletmp%tmp%y(1), bubbletmp%tmp%y(2), &
                                    bubbletmp%tmp%shell, bubbletmp%Rbuck, bubbletmp%Rrupt)
        call s_get_cell(bubbletmp%tmp%s, cell)
        pinf = f_pressure_inf(bubbletmp%tmp%s, q_prim(E_idx), 1, aux1, aux2) ! getting p_inf
        call s_convert_to_mixture_variables(q, cell(1), cell(2), cell(3), rhol, gamma, pi_inf, qv, Re)
        if (solverapproach == 1) rhol = rhol/q_particle(1)%sf(cell(1), cell(2), cell(3))

        ! Computing speed of sound
        do i = 1, E_idx - mom_idx%beg
            vel(i) = q_prim(i + cont_idx%end)%sf(cell(1), cell(2), cell(3))
        end do
        E = gamma*pinf + pi_inf + 0.5d0*rhol*dot_product(vel, vel)
        H = (E + pinf)/rhol
        cson = sqrt((H - 0.5d0*dot_product(vel, vel))/gamma)

        deltaP = pliqint - pinf
        !print*, 'delta pressure:', deltaP, q_prim(E_idx)%sf(cell(1), cell(2), cell(3))
        termI = 0.0d0
        velint = bubbletmp%tmp%y(2) - bubbletmp%dbdt(step)%dmvdt/(4.0d0*pi*bubbletmp%tmp%y(1)**2*rhol)
        bubbletmp%dbdt(step)%dydt(2) = ((1.0d0 + velint/cson)*deltaP/rhol + termI &
                                        + bubbletmp%dbdt(step)%dpbdt*bubbletmp%tmp%y(1)/rhol/cson &
                                             - velint**2*3.0d0/2.0d0*(1.0d0 - velint/3.0d0/cson)) &
                                                         /(bubbletmp%tmp%y(1)*(1.0d0 - velint/cson))

    end subroutine s_compute_RP

    !>  This subroutine computes the fluxes at the bubbles' interface
        !! @param bubbletmp Current bubble information
        !! @param vaporflux Mass flux
        !! @param heatflux Heat flux
        !! @param gammabubble Specific heat of the vapor-gas mixture in the bubble
    subroutine s_compute_interface_fluxes(bubbletmp, vaporflux, heatflux, gammabubble)

        real(kind(0.d0)) :: vaporflux, heatflux, concvapint, bubbleTemp, volbubble, kbubble, &
                            avgconc, Rmixt, gammabubble, rhogas
        real(kind(0.d0)) :: deltaTemp, deltaConc, act_massflag
        type(particledata) :: bubbletmp

        act_massflag = max(dble(massflag) - dble(bubbletmp%tmp%shell), 0d0)

        if (act_massflag == 0) then
            concvapint = 0.d0
        else
            concvapint = (Rvap/Rgas)*(bubbletmp%tmp%p/pvap - 1.0d0)
            concvapint = 1.0d0/(1.0d0 + concvapint)
        end if

        volbubble = 4.0d0/3.0d0*pi*bubbletmp%tmp%y(1)**3d0
        bubbleTemp = bubbletmp%tmp%p*volbubble/(bubbletmp%mg*Rgas + bubbletmp%tmp%mv*Rvap)
        deltaTemp = bubbleTemp - Thost
        gammabubble = concvapint*gammavapor + (1.0d0 - concvapint)*gammagas
        heatflux = -(gammabubble - 1.0d0)/gammabubble*bubbletmp%betaT*deltaTemp/bubbletmp%tmp%y(1)

        avgconc = bubbletmp%tmp%mv/(bubbletmp%mg + bubbletmp%tmp%mv)
        Rmixt = concvapint*Rvap + (1.0d0 - concvapint)*Rgas
        concvapint = min(concvapint, 0.99d0)
        vaporflux = (1.0d0 - concvapint)*bubbletmp%tmp%y(1)
        rhogas = (bubbletmp%mg + bubbletmp%tmp%mv)/(4.0d0/3.0d0*pi*bubbletmp%tmp%y(1)**3d0)
        deltaConc = avgconc - concvapint
        vaporflux = -diffcoefvap*bubbletmp%betaC*deltaConc*rhogas/vaporflux

        ! print '(A,7E24.17)', 'temperature', deltaTemp, deltaConc, bubbletmp%tmp%y(1)

    end subroutine s_compute_interface_fluxes

    !>  This subroutine updates the particle variables
        !! @param dt Time step
        !! @param RKstep Current time step in the adaptative stepper
        !! @param RKcoef 4th/5th order Runge-Kutta-Cash-Karp coefficient
        !! @param largestep Logical variable to determine if the adaptative time step is too large
        !! @param q Conservative variables
        !! @param dq Calculated change of conservative variables
        !! @param q_prim Primitive variables
    subroutine s_update_particle(dt, RKstep, RKcoef, largestep, q, dq, q_prim, shellFlagtmp)

        type(particlenode), pointer :: particle
        real(kind(0.d0)) :: dt
        real(kind(0.d0)), dimension(6), intent(IN) :: RKcoef
        integer, intent(IN) :: RKstep
        logical, optional :: shellFlagtmp
        type(vector_field), dimension(:), optional :: q
        type(vector_field), dimension(:), optional :: dq
        type(scalar_field), dimension(:), optional :: q_prim

        integer, dimension(3) :: oldcell, newcell
        integer :: i, j, ierr
        logical :: largestep, change, indomain
        type(cellwb), pointer :: cellwbaux
        type(particlenode), pointer :: nodeaux
        integer, dimension(3) :: cell
        real(kind(0.d0)) radiusOld, velOld
        logical :: aux

        particle => particlesubList%List%next

        do while (associated(particle))

            call s_get_cell(particle%data%tmp%s, oldcell)
            call s_transfer_data_to_tmp(particle%data)
            if ( present(shellFlagtmp)) particle%data%tmp%shell=particle%data%shell

            radiusOld = particle%data%tmp%y(1)
            velOld = particle%data%tmp%y(2)

            do i = 1, RKstep
                particle%data%tmp%y(1:2) = particle%data%tmp%y(1:2) + dt*RKcoef(i)*particle%data%dbdt(i)%dydt(1:2)
                particle%data%tmp%x(1:3) = particle%data%tmp%x(1:3) + dt*RKcoef(i)*particle%data%dbdt(i)%dxdt(1:3)
                particle%data%tmp%u(1:3) = particle%data%tmp%u(1:3) + dt*RKcoef(i)*particle%data%dbdt(i)%dudt(1:3)
                particle%data%tmp%p = particle%data%tmp%p + dt*RKcoef(i)*particle%data%dbdt(i)%dpbdt
                particle%data%tmp%mv = particle%data%tmp%mv + dt*RKcoef(i)*particle%data%dbdt(i)%dmvdt
            end do

            if ( &!(abs(particle%data%tmp%y(2))/csonhost >= 1.0d0) .or. &  ! no supersonic interface speed
                (particle%data%tmp%y(1) <= 0.0d0) .or. &                 ! no negative radius
                (particle%data%tmp%x(1) /= particle%data%tmp%x(1))) then ! finite bubble location

                print*, 'Large time step. Radius from:', radiusOld, ' to :', particle%data%tmp%y(1), &
                                        '; and velocity from', velOld, ' to :', particle%data%tmp%y(2)
                largestep = .true.
                if (dt < 5.d-14) then
                    print *, 'WARNING large step: removing particle', dt, particle%data%id
                    call s_remove_particle(particle, 1)
                    call s_remove_particle(particle, 2)
                    goto 710
                end if
                goto 711
            end if

            indomain = particle_in_domain(particle%data%tmp%x)

            if (.not. indomain) then
                print *, 'not in domain', particle%data%id, particle%data%tmp%x(1), particle%data%xprev(1), x_cb(-buff_size - 1), x_cb(n + buff_size)
                call s_remove_particle(particle, 2)
                goto 710
            end if

            if (particle%data%equilibrium) then
                call s_equilibrium_state(particle%data, largestep, q_prim(E_idx))
                if (largestep) goto 711
            end if

            particle => particle%next
710         continue
        end do

711     if (num_procs > 1) then
            call MPI_ALLREDUCE (largestep, aux, 1 , MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
            largestep = aux
        end if

        if (largestep) return

        ! Update fluid variables
        if (present(q)) then
            do i = 1, sys_size
                q(2)%vf(i)%sf(0:m, 0:n, 0:p) = q(1)%vf(i)%sf(0:m, 0:n, 0:p)
                do j = 1, RKstep
                    q(2)%vf(i)%sf(0:m, 0:n, 0:p) = q(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                                                   + dt*RKcoef(j)*dq(j)%vf(i)%sf(0:m, 0:n, 0:p)
                end do
            end do
        end if

    end subroutine s_update_particle

    !>  This subroutine calculates the equilibrium state of the lagrangian bubbles
        !! @param particle Variables of the particle
        !! @param largestep Logical variable to determine if the adaptative time step is too large
        !! @param pres Pressure surrounding the bubble
    subroutine s_equilibrium_state(particle, largestep, pres)

        type(particledata) :: particle
        real(kind(0.d0)) :: req(1)
        real(kind(0.d0)) :: pinf, aux1, aux2, sigmaTotal, prinVar1, prinVar2
        integer, dimension(3) :: cell
        logical :: largestep, cond
        type(scalar_field) :: pres

        call s_get_cell(particle%tmp%s, cell)
        pinf = f_interpolate(particle%tmp%s, pres)

        req(1) = particle%tmp%y(1)
        prinVar1 = particle%tmp%y(1)
        prinVar2 = particle%tmp%p
        cond = .false.
        call s_get_equilibrium_radius(10000, req, 1.0d-10, 1.0d-10, particle, pinf, cond)

        particle%tmp%y(1) = req(1)
        if ((req(1) <= 0.0d0) .or. cond) then
            print *, 'released', particle%id
            print*, 'In equilibrium state : RELEASING PARTICLE ', prinVar1, particle%tmp%y(1), prinVar2, particle%tmp%p, pinf
            particle%equilibrium = .false.
            largestep = .true.
            print*, 'Large step in equilibrium state subroutine'
            return
        end if
        if (particle%tmp%shell .eq. 1) then    !Marmotant model
            sigmaTotal = surfaceElast_lipidCoat*((particle%tmp%y(1)/particle%Rbuck)**2-1.0d0)
            if (particle%tmp%y(1)<particle%Rbuck) sigmaTotal = 0.d0
        else
            sigmaTotal = sigmabubble
        end if
        particle%tmp%p = pinf - pvap + 2.0d0*sigmaTotal/req(1)
        !particle%tmp%mv = pvap*4.0d0/3.0d0*pi*req(1)**3*(1/Rvap)
        particle%tmp%mv = pvap*4.0d0/3.0d0*pi*req(1)**3*(1.0d0/(Rvap*Thost))
        
        if (prinVar2 /= particle%tmp%p) print*, 'In equilibrium state : PRESSURE MISMATCH ', prinVar1, particle%tmp%y(1), prinVar2, particle%tmp%p, pinf
        if (prinVar1 /= particle%tmp%y(1)) print*, 'In equilibrium state : RADIUS MISMATCH ', prinVar1, particle%tmp%y(1), prinVar2, particle%tmp%p, pinf

    end subroutine s_equilibrium_state

    !> The purpose of this procedure is to calculate the equilibrium radius of the bubble.
        !! @param ntrial Number of Newton-Rhapson steps to improve the root
        !! @param x Initial guess
        !! @param tolx Tolerance error of calculating x
        !! @param tolf Tolerance error of calculating the bubble equation
        !! @param particle Bubble data
        !! @param pliq Pressure at infinity
        !! @param largestep Adaptative time step switch
    subroutine s_get_equilibrium_radius(ntrial, x, tolx, tolf, particle, pliq, largestep)

        ! uses lubksb, ludcmp usrfun
        !       Given an initial guess x for a root in n dimensions, take ntrial
        !       Newton-Rhapson steps to improve the root. Step if the root
        !       converges in either summed absolute variable increments tolx or
        !       summed absolute function values tolf

        integer, parameter :: np = 15, n = 1
        integer :: ntrial, i, k, indx(np)
        type(particledata) :: particle
        real(kind(0.d0)) :: tolf, tolx, x(n), pliq, sigmaTotal
        real(kind(0.d0)) :: d, errf, errx, fjac(np, np), fvec(np), p(np)
        logical :: largestep

        if (particle%tmp%shell .eq. 1) then    !Marmotant model
            sigmaTotal = surfaceElast_lipidCoat*((particle%tmp%y(1)/particle%Rbuck)**2-1.0d0)
            if (particle%tmp%y(1)<particle%Rbuck) sigmaTotal = 0.d0
        else
            sigmaTotal = sigmabubble
        end if

        do k = 1, ntrial
            ! Supply the values of the function at x in fvec and the Jacobian Matrix at fjac
            !fvec(1) = 4.0d0*pi*x(1)*(pliq - pvap) - 3.0d0*particle%mg*(Rgas)/x(1)**2 + 8.0d0*pi*sigmaTotal
            fvec(1) = 4.0d0*pi*x(1)*(pliq - pvap) - 3.0d0*particle%mg*(Rgas*Thost)/x(1)**2 + 8.0d0*pi*sigmaTotal
            !fjac(1, 1) = 4.0d0*pi*(pliq - pvap) + 6.0d0*particle%mg*(Rgas)/x(1)**3
            fjac(1, 1) = 4.0d0*pi*(pliq - pvap) + 6.0d0*particle%mg*(Rgas*Thost)/x(1)**3

            errf = 0.
            do i = 1, n
                errf = errf + abs(fvec(i))
            end do
            if (errf <= tolf) return
            do i = 1, n
                p(i) = -fvec(i)
            end do
            call ludcmp(fjac, indx, d, largestep)
            if (largestep) return
            call lubksb(fjac, indx, p)
            errx = 0.
            do i = 1, n
                errx = errx + abs(p(i))
                x(i) = x(i) + p(i)
            end do
            if (errx <= tolx) return
        end do

    end subroutine s_get_equilibrium_radius

    !> The purpose of the next two procedures is to solve a nonlinear system of equations which is 
        !!      utilized to get the equilibrium radius of the lagrangian bubbles. It solves 
        !!      for the radius in the Laplace equation.
        !! @param a Input matrix
        !! @param indx Dimensions of matrix a
        !! @param d Output
    subroutine ludcmp(a, indx, d, cond)

        integer :: indx(2)
        integer, parameter :: nmax = 500, n = 1, np = 15
        real(kind(0.d0)), parameter :: tiny = 1.0d-20
        real(kind(0.d0)) :: d, a(15, 15)
        integer :: i, imax, j, k
        real(kind(0.d0)) :: aamax, dum, sum, vv(nmax)
        logical :: cond

        d = 1.
        do i = 1, n
            aamax = 0.
            do j = 1, n
                if (abs(a(i, j)) > aamax) aamax = abs(a(i, j))
            end do
            if (aamax == 0.) then
                cond = .true.
                return
            end if
            vv(i) = 1./aamax
        end do

        do j = 1, n
            do i = 1, j - 1
                sum = a(i, j)
                do k = 1, i - 1
                    sum = sum - a(i, k)*a(k, j)
                end do
                a(i, j) = sum
            end do
            aamax = 0.
            do i = j, n
                sum = a(i, j)
                do k = 1, j - 1
                    sum = sum - a(i, k)*a(k, j)
                end do
                a(i, j) = sum
                dum = vv(i)*abs(sum)
                if (dum >= aamax) then
                    imax = i
                    aamax = dum
                end if
            end do
            if (j /= imax) then
                do k = 1, n
                    dum = a(imax, k)
                    a(imax, k) = a(j, k)
                    a(j, k) = dum
                end do
                d = -d
                vv(imax) = vv(j)
            end if
            indx(j) = imax
            if (a(j, j) == 0.) a(j, j) = tiny
            if (j /= n) then
                dum = 1./a(j, j)
                do i = j + 1, n
                    a(i, j) = a(i, j)*dum
                end do
            end if
        end do

    end subroutine ludcmp

    subroutine lubksb(a, indx, b)

        integer, parameter :: np = 15, n = 1
        integer :: indx(n)
        real(kind(0.d0)) :: a(np, np), b(n)
        integer :: i, ii, j, ll
        real(kind(0.d0)) :: sum

        ii = 0

        do i = 1, n
            ll = indx(i)
            sum = b(ll)
            b(ll) = b(i)
            if (ii /= 0) then
                do j = ii, i - 1
                    sum = sum - a(i, j)*b(j)
                end do
            elseif (sum /= 0.) then
                ii = i
            end if
            b(i) = sum
        end do

        do i = n, 1, -1
            sum = b(i)
            do j = i + 1, n
                sum = sum - a(i, j)*b(j)
            end do
            b(i) = sum/a(i, i)
        end do

    end subroutine lubksb

    !>  This subroutine calculates the maximum error between the 4th and 5th order Runge-Kutta-Cash-Karp algorithm.
        !!      If the errors are smaller than a tolerance with the same time step size, then the algorithm employs
        !!      the 5th order solution, while if not, both eULERIAN/lagrangian variables are re-calculated with a 
        !!      smaller time step size.
        !! @param timetmp Current time
        !! @param dt Time discretization
        !! @param RKcoef 4th/5th order Runge-Kutta-Cash-Karp coefficient
        !! @param errmax Maximum error
        !! @param t_step Time step
        !! @param q Conservative variables
        !! @param q_prim Primitive variables
        !! @param dq Calculated change of conservative variables
    subroutine s_calculate_RKerror(timetmp, dt, RKcoef, errmax, t_step, q, q_prim, dq)

        type(particlenode), pointer :: particle
        real(kind(0.d0)), intent(IN) :: dt, timetmp
        real(kind(0.d0)), dimension(6), intent(IN) :: RKcoef
        real(kind(0.d0)) :: errmax, erraux, errb
        type(vector_field), dimension(:), optional :: q
        type(scalar_field), dimension(:), optional :: q_prim
        type(vector_field), dimension(:), optional :: dq
        integer :: i, j, k, l, l1, nb
        integer, intent(IN) :: t_step
        real(kind(0.d0)) :: time_avg

        errmax = 0.0d0
        erraux = 0.0d0

        particle => particlesubList%List%next
        do while (associated(particle))
            errb = 0.0d0
            if (.not. particle%data%equilibrium) then
                !Bubble radius error
                do i = 1, 6
                    erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dydt(1)
                end do
                errb = max(errb, abs(erraux)*dt/particle%data%R0)
                !if (errb/RKeps > 1.0d0) print*, 'Truncation error Bubble radius:', particle%data%R0, particle%data%x(1), particle%data%x(2), errb, particle%data%id

                !Interface velocity error
                erraux = 0.0d0
                do i = 1, 6
                    erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dydt(2)
                end do
                errb = max(errb, abs(erraux)*dt)
                !if (errb/RKeps > 1.0d0) print*, 'Truncation error Interface vel:', particle%data%R0, particle%data%x(1), particle%data%x(2), errb, particle%data%id

                !Bubble velocity error
                do j = 1, 3
                    erraux = 0.0d0
                    do i = 1, 6
                        erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dxdt(j)
                    end do
                    errb = max(errb, abs(erraux)*dt/(abs(particle%data%tmp%u(j)) + 1.0d-4))
                end do
            end if
            !if (errb/RKeps > 1.0d0) print*, 'Truncation error bubble velocity:', particle%data%R0, particle%data%x(1), particle%data%x(2), errb, particle%data%id
            errmax = max(errmax, errb)
            
            particle => particle%next
        end do

        if (present(q)) then
            do l1 = 1, cont_idx%end
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            erraux = q(1)%vf(l1)%sf(i, j, k)
                            do l = 1, 6
                                erraux = erraux + dt*RKcoef(l)*dq(l)%vf(l1)%sf(i, j, k)
                            end do
                            erraux = max(errmax, erraux)
                        end do
                    end do
                end do
            end do
            do l1 = mom_idx%beg, mom_idx%beg + num_dims - 1
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            erraux = q(1)%vf(l1)%sf(i, j, k)
                            do l = 1, 6
                                erraux = erraux + dt*RKcoef(l)*dq(l)%vf(l1)%sf(i, j, k)
                            end do
                            erraux = max(errmax, erraux)
                        end do
                    end do
                end do
            end do
            call s_compute_rhs(q(2)%vf, q_prim, dq(6)%vf, t_step=t_step, time_avg=time_avg, qtime=timetmp)
        end if

    end subroutine s_calculate_RKerror

    !>  This subroutine updates the conservative fields after performing the adaptative Runge-Kutta-Cash-Karp time stepper.
        !! @param q Conservative variables
        !! @param update_fields Flag to update the fields
        !! @param q_prim Primitive variables
    subroutine s_update_RK(q, update_fields, q_prim, q_cons_hifu, hdid)

        type(particlenode), pointer :: particle
        type(vector_field), dimension(:), optional :: q
        type(scalar_field), dimension(:), optional :: q_prim
        integer :: i, j, k
        real(kind(0.d0)) :: pinf, pcrit, sigmaTotal
        logical :: update_fields, release_part
        !hifu vars
        type(scalar_field), dimension(sys_size_hifu), intent(inout), optional :: q_cons_hifu
        real(kind(0.d0)), intent(in), optional :: hdid
        real(kind(0.d0)) :: gammabubble, heatflux, bubbletemp, concvapint, qther, qvis, volbubble, stddsv
        real(kind(0.d0)), dimension(3) :: nodecoord
        real(kind(0.d0)) :: tmp, qvis_beforeKernel, qvis_afterKernel, qth_beforeKernel, qth_afterKernel, radBubble, velBubble
        integer, dimension(3) :: cell
        integer :: kernel

        radBubble = 0.0d0
        velBubble = 0.0d0
        qvis_beforeKernel = 0.0d0
        qth_beforeKernel = 0.0d0
        qvis_afterKernel = 0.0d0
        qth_afterKernel = 0.0d0


        if (present(q_cons_hifu)) then !zeroing before smoothening viscous and thermal damping from the bubbles 
            do i = 0,m
                do j = 0,n
                    do k = 0,p
                        q_cons_hifu(qvis_hifu_idx + 1)%sf(i,j,k) = 0.0d0
                        q_cons_hifu(qth_hifu_idx + 1)%sf(i,j,k) = 0.0d0
                    end do
                end do
            end do
        end if

        particle => particlesubList%List%next
        release_part = .false.
        do while (associated(particle))
            if (.not. particle%data%equilibrium) then
                pinf = f_pressure_inf(particle%data%tmp%s, q_prim(E_idx), 1)
                if (particle%data%tmp%shell .eq. 1) then    !Marmotant model
                    sigmaTotal = surfaceElast_lipidCoat*((particle%data%tmp%y(1)/particle%data%Rbuck)**2-1.0d0)
                    if (particle%data%tmp%y(1)<particle%data%Rbuck) sigmaTotal = 0.d0
                else
                    sigmaTotal = sigmabubble
                end if
                pcrit = pvap - &
                        4.0d0*sigmaTotal/(3.d0*sqrt(3.0d0*(pref + .0*sigmaTotal/particle%data%R0)*particle%data%R0**3/(2.0d0*sigmaTotal)))
                pcrit = min(pcrit, -pref)
                if (abs((pcrit - pinf)/pcrit) < 0.5d0) release_part = .true.
            end if
            particle%data%x = particle%data%tmp%x
            particle%data%u = particle%data%tmp%u
            particle%data%y = particle%data%tmp%y
            particle%data%p = particle%data%tmp%p
            particle%data%mv = particle%data%tmp%mv
            if (particle%data%tmp%shell /= particle%data%shell) print*, 'Shell broke from particle, rank', particle%data%id, proc_rank
            particle%data%shell = particle%data%tmp%shell

            !HIFU calculating heat sources
            if (present(q_cons_hifu)) then

                radBubble = particle%data%y(1)
                velBubble = particle%data%y(2)

                volbubble = 4.0d0/3.0d0*pi*particle%data%y(1)**3
                call s_get_cell(particle%data%tmp%s, cell)
                nodecoord(1) = particle%data%x(1)
                nodecoord(2) = particle%data%x(2)
                if (p > 0) nodecoord(3) = particle%data%x(3)
                call s_compute_stddsv(cell, kernel, volbubble, stddsv) ! Kernel function based on the bubble volume

                ! Viscous damping of the bubbles
                qvis = (4.0d0*pi*particle%data%y(1)**2)*(4.0d0*vischost*(particle%data%y(2)**2)/(particle%data%y(1)))
                qvis_beforeKernel = qvis
                call s_smoothfunction ( q_cons_hifu(qvis_hifu_idx+1), nodecoord, cell , qvis, kernel, stddsv) ! Sampling viscous intensity
                

                ! Thermal damping of the bubbles
                if (particle%data%shell .eq. 1) then !no mass transfer
                    concvapint = 0.d0
                else
                    concvapint = (Rvap/Rgas)*(particle%data%p/pvap - 1.0d0)
                    concvapint = 1.0d0/(1.0d0 + concvapint)
                end if
                bubbleTemp = particle%data%mg*Rgas + particle%data%mv*Rvap
                bubbleTemp = particle%data%p*volbubble/bubbleTemp
                gammabubble = concvapint*gammavapor + (1.0d0 - concvapint)*gammagas
                heatflux = -(gammabubble - 1.0d0)/gammabubble*particle%data%betaT*(bubbleTemp - Thost)/particle%data%y(1)
                qther = (4.0d0*pi*particle%data%y(1)**2)*(heatFlux)
                qth_beforeKernel = qther
                call s_smoothfunction ( q_cons_hifu(qth_hifu_idx+1), nodecoord, cell , qther, kernel, stddsv) ! Sampling thermal intensity
            end if

            particle => particle%next
        end do

        if (present(q_cons_hifu)) then !updating viscous and thermal intensity fields
            do i = 0,m
                do j = 0,n
                    do k = 0,p
                        q_cons_hifu(qvis_hifu_idx)%sf(i,j,k) = q_cons_hifu(qvis_hifu_idx)%sf(i,j,k) + hdid * q_cons_hifu(qvis_hifu_idx + 1)%sf(i,j,k)
                        q_cons_hifu(qth_hifu_idx)%sf(i,j,k)  = q_cons_hifu(qth_hifu_idx)%sf(i,j,k)  + hdid * q_cons_hifu(qth_hifu_idx + 1)%sf(i,j,k)
                    end do
                end do
            end do
        end if

        !Writing viscous and thermal damping after and before kernel (intendend for only 1 bubble in the domain in one processor only)
        ! if (present(q_cons_hifu)) then
        !     tmp = radBubble
        !     call s_mpi_allreduce_sum(tmp, radBubble)
        !     tmp = velBubble
        !     call s_mpi_allreduce_sum(tmp, velBubble)

        !     !Before kernel
        !     tmp = qvis_beforeKernel
        !     call s_mpi_allreduce_sum(tmp, qvis_beforeKernel)
        !     tmp = qth_beforeKernel
        !     call s_mpi_allreduce_sum(tmp, qth_beforeKernel)

        !     !After kernel
        !     do i = 0,m
        !         do j = 0,n
        !             do k = 0,p
        !                 qvis_afterKernel = qvis_afterKernel + q_cons_hifu(qvis_hifu_idx+1)%sf(i,j,k) * dx(i)*dy(j)*y_cc(j)*2d0*pi !W
        !                 qth_afterKernel = qth_afterKernel + q_cons_hifu(qth_hifu_idx+1)%sf(i,j,k) * dx(i)*dy(j)*y_cc(j)*2d0*pi    !W 
        !             end do
        !         end do
        !     end do

        !     tmp = qvis_afterKernel
        !     call s_mpi_allreduce_sum(tmp, qvis_afterKernel)
        !     tmp = qth_afterKernel
        !     call s_mpi_allreduce_sum(tmp, qth_afterKernel)

        !     if (proc_rank==0) write (98, '(6x,8e24.8)') &
        !         hdid, &
        !         q_cons_hifu(tt_hifu_idx)%sf(0,0,0) + hdid, &
        !         qvis_beforeKernel, &
        !         qvis_afterKernel, &
        !         qth_beforeKernel, &
        !         qth_afterKernel, &
        !         radBubble, &
        !         velBubble

        ! end if

        ! Releasing particles
        particle => particlesubList%List%next
        if (release_part) then
            do while (associated(particle))
                particle%data%equilibrium = .false.
                particle => particle%next
            end do
        end if

        if ((coupledflag .or. bubblesources) .and. (update_fields)) then
            do i = 1, sys_size
                q(1)%vf(i)%sf = q(2)%vf(i)%sf
            end do
        end if

        if (avgdensFlag) then
            call s_smear_voidfraction(q(1)%vf)
        end if


    end subroutine s_update_RK

    !> This function calculates the heat and mass transfer coefficients following the formulation proposed by Preston et al. (2007)
        !! @param Pe Peclet number of heat or mass transfer
        !! @param omegaN Characteristic frequency
    function f_transfercoeff(Pe, omegaN)

        real(kind(0.d0)) :: f_transfercoeff, Pe, omegaN
        complex :: transferfunc, auxc

        transferfunc = csqrt(cmplx(0.0d0, Pe*omegaN))
        auxc = (cexp(-cmplx(2.0d0, 0.0d0)*transferfunc) + 1.0d0)/(-cexp(-cmplx(2.0d0, 0.0d0)*transferfunc) + 1.0d0)
        transferfunc = transferfunc*auxc - 1.0d0
        transferfunc = 1.0d0/transferfunc - 3.0d0/cmplx(0.0d0, Pe*omegaN)
        f_transfercoeff = dble(1.0d0/transferfunc)

    end function f_transfercoeff

    !> This function calculates the pressure at the bubble wall (bubble-liquid interface).
        !! @param pbubble Pressure inside the bubble
        !! @param radius Bubble radius
        !! @param bubblevel Interface velocity
    subroutine s_compute_pressureliq_int(pressureliq_int, pbubble, radius, bubblevel, shell, Rbuck, Rrupt)

        real(kind(0.d0)), intent(out)   :: pressureliq_int
        real(kind(0.d0)), intent(in)    :: pbubble, radius, bubblevel, Rbuck, Rrupt
        integer, intent(inout) :: shell

        real(kind(0.d0)) :: sigma

        pressureliq_int = 0.d0

        !Adding Marmottant model for lipid coated bubbles
        if (shell==0) then !No shell
            sigma = sigmabubble
        else if (radius <= Rbuck ) then
            sigma = 0.0d0
        else if (radius > Rrupt) then
            sigma = sigmabubble
            shell = 0
        else
            sigma = surfaceElast_lipidCoat*((radius/Rbuck)**2-1.0d0)
        end if


        pressureliq_int = pbubble - 4.0d0*vischost*bubblevel/radius &
                                - 2.0d0*sigma/radius - 4.0d0*shell*surfaceDilatVisc_lipidCoat*bubblevel/(radius**2)
        
        !print*, sigma, pbubble, bubblevel, shell


    end subroutine s_compute_pressureliq_int

    !> The purpose of this procedure is obtain the pressure that drives the bubble oscillations p_inf
        !!      utilizing the model proposed by Maeda and Colonius (2018) to separate p_inf from the
        !!      pressure radiated by the bubble due to its own oscillations.
        !! @param scoord Spatial coordinates
        !! @param pres  Eulerian pressure
        !! @param ptype 1: p at infinity, 2: averaged P at the bubble location
    function f_pressure_inf(scoord, pres, ptype, preterm1, term2, Romega)

        type(particlenode), pointer :: bubble
        type(scalar_field) :: pres
        real(kind(0.d0)), dimension(3) :: distance, center
        real(kind(0.d0)) :: f_pressure_inf, jac, dij, dpotjdt, dc, vol, aux, &
                            volgas, term1, Rbeq, denom, Rmax, stddsv, &
                            charvol, charpres, charvol2, charpres2
        real(kind(0.d0)), dimension(3) :: scoord
        real(kind(0.d0)), optional :: preterm1, term2, Romega
        integer, dimension(3) :: cell, cellaux
        integer, dimension(3) :: epsilonbaux
        integer :: ptype !1=p at infinity, 2= averaged P at the bubble location
        integer :: i, j, k, dir
        logical :: celloutside

        f_pressure_inf = 0.0d0
        call s_get_cell(scoord, cell)

        if ((clusterflag == 0)) then
            !getting p_cell in terms of only the current cell by interpolation

            call s_get_char_vol(cell, vol)
            bubble => qbl%fp(cell(1), cell(2), cell(3))%List%next
            Rmax = 0.0d0
            do while (associated(bubble))
                Rmax = Rmax + bubble%data%tmp%y(1)**3
                bubble => bubble%next
            end do

            ! Surrogate bubble radius
            Rmax = Rmax**(1.0d0/3.0d0)

            ! Getting the cell volulme as Omega
            call s_get_char_dist(cell, stddsv)

            !p_cell (interpolated)
            f_pressure_inf = f_interpolate(scoord, pres)

            !R_Omega
            dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)

        else if (clusterflag == 1) then
            ! Just making the characteristic volume 3^3 times bigger than the cell volume
            ! Not recomended to use for production

            ! Getting Omega (not Omega_L)
            call s_get_char_vol(cell, vol)
            vol = vol*2.7d1 !Just multiplying by 3^3=27
            bubble => qbl%fp(cell(1), cell(2), cell(3))%List%next
            Rmax = 0.0d0
            do while (associated(bubble))
                Rmax = Rmax + bubble%data%tmp%y(1)**3
                bubble => bubble%next
            end do

            ! Surrogate bubble radius
            Rmax = Rmax**(1.0d0/3.0d0)

            ! Getting the cell volulme as Omega
            call s_get_char_dist(cell, stddsv)

            !p_cell to get interpolate
            f_pressure_inf = f_interpolate(scoord, pres)

            !R_Omega
            dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)

        else if (clusterflag >= 2) then
            ! Bubble dynamic closure from Kazuki and Colonius (2018)

            ! Range of cells included in Omega
            if (smoothtype == 1) then
                epsilonbaux(:) = 3
            else if (projectiontype == 1) then
                epsilonbaux(:) = 3
            end if

            charvol = 0.d0
            charpres = 0.d0
            charvol2 = 0.d0
            charpres2 = 0.d0
            vol = 0.d0
            if (num_dims == 3) then
                k = -epsilonbaux(3)
            else
                k = 0
            end if
            i = -epsilonbaux(1); j = -epsilonbaux(2)

3001        if ((i <= epsilonbaux(1)) .and. (j <= epsilonbaux(2))) then
                celloutside = .false.
                cellaux(1) = cell(1) + i
                cellaux(2) = cell(2) + j
                cellaux(3) = cell(3) + k

                !Check ghost part in x-direction
                if (cellaux(1) < -buff_size) then
                    celloutside = .true.
                    i = i + 1
                end if

                !Check ghost part in y-direction
                if (cellaux(2) < -buff_size) then
                    celloutside = .true.
                    j = j + 1
                end if
                if (cyl_coord .and. (num_dims /= 3)) then
                    if ((cellaux(2) < n+buff_size) .and. (.not. celloutside)) then 
                        if (y_cc(cellaux(2)) < 0d0) then
                            celloutside = .true.
                            j = j + 1
                        end if
                    end if
                end if

                !Check ghost part in z-direction
                if (num_dims == 3) then
                    if (cellaux(3) < -buff_size) then
                        celloutside = .true.
                        k = k + 1
                    end if
                end if

                if (cellaux(1) > m + buff_size) celloutside = .true.
                if (cellaux(2) > n + buff_size) celloutside = .true.
                if (cellaux(3) > p + buff_size) celloutside = .true.

                if (.not. celloutside) then
                    call s_get_char_vol(cellaux, vol)
                    charvol = charvol + vol
                    charpres = charpres + pres%sf(cellaux(1), cellaux(2), cellaux(3))*vol
                    charvol2 = charvol2 + vol*q_particle(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                    charpres2 = charpres2 + pres%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                *vol*q_particle(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                end if

                if (j < epsilonbaux(2)) then
                    j = j + 1
                    goto 3001
                end if

3002            j = -epsilonbaux(2)
                i = i + 1
                goto 3001

            end if

3003        if ((num_dims == 3) .and. (k < epsilonbaux(3))) then
                k = k + 1
                i = -epsilonbaux(1); j = -epsilonbaux(2)
                goto 3001
            end if
            f_pressure_inf = charpres2/charvol2 ! p_cell eqn 30 (kazuki's paper)
            vol = charvol
            dc = (3.0d0*abs(vol)/(4.0d0*pi))**(1.0d0/3.0d0)  !positive volume
        else

            print '(A)', 'Check clusterflag. Exiting ...'
            call s_mpi_abort()

        end if

        if (correctpresFlag .and. present(preterm1)) then
            dpotjdt = 0.0d0 !potential derivative contribution from other bubbles
            volgas = 0.0d0
            term1 = 0.0d0
            term2 = 0.0d0
            denom = 0.0d0
            bubble => qbl%fp(cell(1), cell(2), cell(3))%List%next

            do while (associated(bubble))
                volgas = volgas + bubble%data%tmp%y(1)**3 !surrogate bubble volume
                denom = denom + bubble%data%tmp%y(1)**2
                term1 = term1 + bubble%data%dphidt*bubble%data%tmp%y(1)**2
                term2 = term2 + bubble%data%tmp%y(2)*bubble%data%tmp%y(1)**2
                bubble => bubble%next
            end do

            Rbeq = volgas**(1.0d0/3.0d0) !surrogate bubble radius
            aux = dc**3 - Rbeq**3
            term2 = term2/denom
            !term2 = 3.0d0/2.0d0*term2**2*Rbeq**3*(1.0d0 - Rbeq/dc)/aux
            term2 = 3.0d0/4.0d0*term2**2*Rbeq**3*(1.0d0 - Rbeq/dc)/aux
            preterm1 = 3.0d0/2.0d0*Rbeq*(dc**2 - Rbeq**2)/(aux*denom)

            !Control volume radius
            if (present(Romega)) Romega = dc

            ! Getting p_inf
            if (ptype == 1) then
                !print*, f_pressure_inf, preterm1, term1, term2, f_pressure_inf + preterm1*term1 + term2
                f_pressure_inf = f_pressure_inf + preterm1*term1 + term2 ! eqn 80 (kazuki's pape) need: term2 * rho
            end if

        end if

    end function f_pressure_inf

    !> This function performs a bilinear interpolation.
          !! @param coord Interpolation coordintes
          !! @param q     Input scalar field
    function f_interpolate(coord, q)

        type(scalar_field), intent(in) :: q
        real(kind(0.d0)), dimension(3), intent(in) :: coord
        real(kind(0.d0)) :: f_interpolate, tmp
        real(kind(0.d0)), dimension(3) :: psi !local coordinates
        integer, dimension(3) :: cell

        call s_get_psi(coord, psi, cell)

        if (p == 0) then  !2D
            tmp = q%sf(cell(1), cell(2), cell(3))*(1.0d0 - psi(1))*(1.0d0 - psi(2))
            tmp = tmp + q%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1.0d0 - psi(2))
            tmp = tmp + q%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)
            tmp = tmp + q%sf(cell(1), cell(2) + 1, cell(3))*(1.0d0 - psi(1))*psi(2)
        else              !3D
            tmp = q%sf(cell(1), cell(2), cell(3))*(1.0d0 - psi(1))*(1.0d0 - psi(2))*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1.0d0 - psi(2))*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1), cell(2) + 1, cell(3))*(1.0d0 - psi(1))*psi(2)*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1), cell(2), cell(3) + 1)*(1.0d0 - psi(1))*(1.0d0 - psi(2))*psi(3)
            tmp = tmp + q%sf(cell(1) + 1, cell(2), cell(3) + 1)*psi(1)*(1.0d0 - psi(2))*psi(3)
            tmp = tmp + q%sf(cell(1) + 1, cell(2) + 1, cell(3) + 1)*psi(1)*psi(2)*psi(3)
            tmp = tmp + q%sf(cell(1), cell(2) + 1, cell(3) + 1)*(1.0d0 - psi(1))*psi(2)*psi(3)
        end if

        f_interpolate = tmp

    end function f_interpolate

    !> This subroutine returns the computational coordinate of the cell for the given position.
          !! @param pos     Input coordintes
          !! @param cell    Computational coordinate of the cell
          !! @param scoord  Calculated particle coordintes
    subroutine s_locate_cell(pos, cell, scoord, radPos)

        real(kind(0.d0)), dimension(3) :: pos
        real(kind(0.d0)), dimension(3), optional :: scoord
        real(kind(0d0)), optional :: radPos
        integer, dimension(3) :: cell
        integer :: i, j, k

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
            if (cell(2) > n + buff_size) print*, 'in locate cell:', radPos, y_cb(n + buff_size), y_cb(-buff_size - 1)
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
        if (present(scoord)) then
            scoord(1) = cell(1) + (pos(1) - x_cb(cell(1) - 1))/dx(cell(1))
            scoord(2) = cell(2) + (pos(2) - y_cb(cell(2) - 1))/dy(cell(2))
            scoord(3) = 0.0d0
            if (p > 0) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3) - 1))/dz(cell(3))
            call s_get_cell(scoord, cell)
        end if

    end subroutine s_locate_cell

    !> This subroutine transfer data into the temporal list in the lagrangian solver
          !! @param particle_data  Particle data
    subroutine s_transfer_data_to_tmp(particle_data)

        type(particledata) :: particle_data

        particle_data%tmp%x = particle_data%x
        particle_data%tmp%u = particle_data%u
        particle_data%tmp%y = particle_data%y
        particle_data%tmp%p = particle_data%p
        particle_data%tmp%mv = particle_data%mv
        !particle_data%tmp%shell = particle_data%shell

    end subroutine s_transfer_data_to_tmp

    !> The purpose of this procedure is to calculate the characteristic cell distance
        !! @param cell Computational coordinates
        !! @param Chardist Characteristic distance
    subroutine s_get_char_dist(cell, Chardist)

        real(kind(0.d0)) :: Chardist
        integer, dimension(3) :: cell

        if (p > 0) then
            Chardist = (dx(cell(1))*dy(cell(2))*dz(cell(3)))**(1./3.)
        else
            Chardist = sqrt(dx(cell(1))*dy(cell(2)))
        end if

    end subroutine s_get_char_dist

    !> The purpose of this procedure is to calculate the characteristic cell volume
        !! @param cell Computational coordinates
        !! @param Chardist Characteristic volume
    subroutine s_get_char_vol(cell, Charvol)

        real(kind(0.d0)) :: Charvol
        integer, dimension(3) :: cell

        if (p > 0) then
            Charvol = dx(cell(1))*dy(cell(2))*dz(cell(3))
        else
            if (cyl_coord) then
                Charvol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2d0*pi
            else
                Charvol = dx(cell(1))*dy(cell(2))*charwidth
            end if

        end if

    end subroutine s_get_char_vol

    !> The purpose of this procedure is to determine if the global coordinates of the bubbles are present in the current MPI processor
        !! @param pos_part Spatial coordinates of the bubble
    function particle_in_domain(pos_part, inputCoord)

        logical :: particle_in_domain
        real(kind(0.d0)), dimension(3) :: pos_part
        logical, optional :: inputCoord
        real(kind(0d0)) :: radialPosition

        ! 2D
        if (p == 0 .and. cyl_coord .neqv. .true.) then
            ! Defining a virtual z-axis that has the same dimensions as y-axis
            ! defined in the input file
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < charwidth/2d0) .and. (pos_part(3) >= -charwidth/2d0))
        end if
        
        ! 2D cyl_coord
        if (p == 0 .and. cyl_coord .eqv. .true.) then
            
            !particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
            !                      (abs(pos_part(2)) < y_cb(n + buff_size)) .and. (abs(pos_part(2)) >= max(y_cb(-buff_size - 1), 0d0)))
            
            if (present(inputCoord)) then
                radialPosition = sqrt(pos_part(2)**2 + pos_part(3)**2)
            else
                radialPosition = pos_part(2)
            end if

            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                        (radialPosition < y_cb(n + buff_size)) .and. (radialPosition >= max(y_cb(-buff_size - 1), 0d0)))

            if (particle_in_domain .and. (radialPosition > y_cb(n + buff_size))) print*, 'Error in particle in domain:', radialPosition, y_cb(n + buff_size), y_cb(-buff_size - 1)
        end if

        ! 3D
        if (p > 0) then
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < z_cb(p + buff_size)) .and. (pos_part(3) >= z_cb(-buff_size - 1)))
        end if

        ! For symmetric boundary condition
        if (bc_x%beg == -2) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) >= x_cb(-1)))
        end if
        if (bc_x%end == -2) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) < x_cb(m)))
        end if
        if (bc_y%beg == -2 .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) >= y_cb(-1)))
        end if
        if (bc_y%end == -2 .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) < y_cb(n)))
        end if

        if (p > 0) then
            if (bc_z%beg == -2) then
                particle_in_domain = (particle_in_domain .and. (pos_part(3) >= z_cb(-1)))
            end if
            if (bc_z%end == -2) then
                particle_in_domain = (particle_in_domain .and. (pos_part(3) < z_cb(p)))
            end if
        end if

    end function particle_in_domain

    function particle_in_domain_physical(pos_part)
        ! This subroutine is used for mpi_parallel_io
        logical :: particle_in_domain_physical
        real(kind(0.d0)), dimension(3) :: pos_part

        particle_in_domain_physical = ((pos_part(1) < x_cb(m)) .and. (pos_part(1) >= x_cb(-1)) .and. &
                                       (pos_part(2) < y_cb(n)) .and. (pos_part(2) >= y_cb(-1)))

        if (p > 0) then
            particle_in_domain_physical = (particle_in_domain_physical .and. (pos_part(3) < z_cb(p)) .and. (pos_part(3) >= z_cb(-1)))
        end if

    end function particle_in_domain_physical

    !> The purpose of this procedure is the gradient of a scalar field along the x, y and z directions following a 
        !!      second-order central difference considering uneven widths
        !! @param q Input scalar field
        !! @param dq Output gradient of q
        !! @param dir Gradient spatial direction
    subroutine s_gradient_dir(q, dq, dir)

        type(scalar_field) :: q, dq
        integer :: dir, i, j, k, l, lmax
        real(kind(0.d0)) :: aux1, aux2

        if (dir == 1) then
            ! Gradient in x dir.
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        aux1 = dx(i) + dx(i - 1)
                        aux2 = dx(i) + dx(i + 1)
                        dq%sf(i, j, k) = q%sf(i, j, k)*(dx(i + 1) - dx(i - 1)) &
                                         + q%sf(i + 1, j, k)*aux1 &
                                         - q%sf(i - 1, j, k)*aux2
                        dq%sf(i, j, k) = dq%sf(i, j, k)/(aux1*aux2)
                    end do
                end do
            end do
        else
            if (dir == 2) then
                ! Gradient in y dir.
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            aux1 = dy(j) + dy(j - 1)
                            aux2 = dy(j) + dy(j + 1)
                            dq%sf(i, j, k) = q%sf(i, j, k)*(dy(j + 1) - dy(j - 1)) &
                                             + q%sf(i, j + 1, k)*aux1 &
                                             - q%sf(i, j - 1, k)*aux2
                            dq%sf(i, j, k) = dq%sf(i, j, k)/(aux1*aux2)
                        end do
                    end do
                end do
            else
                ! Gradient in z dir.
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            aux1 = dz(k) + dz(k - 1)
                            aux2 = dz(k) + dz(k + 1)
                            dq%sf(i, j, k) = q%sf(i, j, k)*(dz(k + 1) - dz(k - 1)) &
                                             + q%sf(i, j, k + 1)*aux1 &
                                             - q%sf(i, j, k - 1)*aux2
                            dq%sf(i, j, k) = dq%sf(i, j, k)/(aux1*aux2)
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_gradient_dir

    !> Subroutine that writes on each time step the changes of the lagrangian bubbles, and is activated 
        !!  with particleoutFlag.
        !!  @param q_time Current time
    subroutine s_write_particles(qtime)

        type(particlenode), pointer :: particle
        real(kind(0.d0)) :: qtime, pinf
        integer :: i
        integer, dimension(3) :: cell

        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'particles_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        if (qtime == 0.d0) then
            open (11, FILE=trim(file_loc), FORM='formatted', position='rewind')
            write (11, *) 'currentTime, particleID, x, y, z, ', &
                'coreVaporMass, coreVaporConcentration, radius, interfaceVelocity, ', &
                'corePressure'
        else
            open (11, FILE=trim(file_loc), FORM='formatted', position='append')
        end if

        ! Cycle through list
        particle => particlesubList%List%next
        do while (associated(particle))

            write (11, '(6X,E12.6,I24.8,8E24.8,I24.8)') &
                    qtime, &
                    particle%data%id, &
                    particle%data%tmp%x(1), &
                    particle%data%tmp%x(2), &
                    particle%data%tmp%x(3), &
                    particle%data%tmp%mv, &
                    particle%data%tmp%mv/(particle%data%tmp%mv + particle%data%mg), &
                    particle%data%tmp%y(1), &
                    particle%data%tmp%y(2), &
                    particle%data%tmp%p,    &
                    particle%data%tmp%shell

            particle => particle%next
        end do

        close (11)

    end subroutine s_write_particles

    !>  Subroutine that writes some useful statistics related to the volume fraction
            !!       of the particles (void fraction) in the computatioational domain
            !!       on each time step.
            !!  @param q_time Current time
    subroutine s_write_void_evol(qtime)

        real(kind(0.d0)) :: qtime, voidmax, voidavg, vol, volcell, voltot
        real(kind(0.d0)) :: voidmax_glb, voidavg_glb, vol_glb
        integer :: i, j, k
        integer, dimension(3) :: cell
        logical :: prevfile

        character(LEN=path_len + 2*name_len) :: file_loc

        if (proc_rank == 0) then
            write (file_loc, '(A)') 'voidfraction.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            if (qtime == 0.d0) then
                open (12, FILE=trim(file_loc), FORM='formatted', position='rewind')
                write (12, *) 'currentTime, averageVoidFraction, ', &
                    'maximumVoidFraction, totalParticlesVolume'
                write (12, *) 'The averageVoidFraction value does ', &
                    'not reflect the real void fraction in the cloud since the ', &
                    'cells which do not have bubbles are not accounted'
            else
                open (12, FILE=trim(file_loc), FORM='formatted', position='append')
            end if
        end if

        voidmax = 0.0d0
        voidavg = 0.0d0
        vol = 0.0d0

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    cell(1) = i
                    cell(2) = j
                    cell(3) = k
                    voidmax = max(voidmax, 1 - q_particle(1)%sf(i, j, k))
                    call s_get_char_vol(cell, volcell)
                    if ((1.-q_particle(1)%sf(i, j, k)) > 5e-11) then
                        voidavg = voidavg + (1 - q_particle(1)%sf(i, j, k))*volcell
                        vol = vol + volcell
                    end if
                end do
            end do
        end do

        if (num_procs > 1) then
            call s_mpi_allreduce_max(voidmax, voidmax_glb)
            voidmax = voidmax_glb
            call s_mpi_allreduce_sum(vol, vol_glb)
            vol = vol_glb
            call s_mpi_allreduce_sum(voidavg, voidavg_glb)
            voidavg = voidavg_glb
        end if

        voltot = voidavg
        ! This voidavg value does not reflect the real void fraction in the cloud
        ! since the cell which does not have bubbles are not accounted
        if (vol > 0.) voidavg = voidavg/vol

        if (proc_rank == 0) then
            write (12, '(6X,f12.6,3f24.8)') &
                    qtime, &
                    voidavg, &
                    voidmax, &
                    voltot
            close (12)
        end if

    end subroutine s_write_void_evol

    !>  Subroutine that writes the restarting files for the particles in the lagrangian solver.
        !!  @param t_step Current time step
    subroutine s_write_restart_particles(t_step) !_parallel

        ! Location of time-step folder corresponding to time-step, t_step
        character(LEN=len_trim(case_dir) + 2*name_len) :: t_step_dir

        ! Generic string used to store the address of a particular file
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc

        ! Generic logical used for purpose of asserting whether a particular
        ! directory is or is not located in the designated location
        logical :: dir_check
        logical :: file_exist

        type(particlenode), pointer :: particle
        integer :: i, t_step
        integer :: nparticles, tot_part, tot_part_wrtn, npart_wrtn

#ifdef MFC_MPI
        ! For Parallel I/O
        integer :: ifile, ireq, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer :: view
        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer, dimension(num_procs) :: part_order, part_ord_mpi
        integer :: varsMarmotant = 3

        nparticles = 0d0
        if (particlesubList%nb /= 0) then
            particle => particlesubList%List%next
            do while (associated(particle))
                if (particle_in_domain_physical(particle%data%x(1:3))) then
                    nparticles = nparticles + 1
                end if
                particle => particle%next
            end do
        end if

        ! Total number of particles
        call MPI_ALLREDUCE(nparticles, tot_part, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        ! Total number of particles written so far
        call MPI_ALLREDUCE(npart_wrtn, tot_part_wrtn, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        lsizes(1) = max(1, nparticles)
        lsizes(2) = 21 + varsMarmotant

        ! if the partcle number is zero, put 1 since MPI cannot deal with writing
        ! zero particle
        part_order(:) = 1
        part_order(proc_rank + 1) = max(1, nparticles)

        call MPI_ALLREDUCE(part_order, part_ord_mpi, num_procs, MPI_integer, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

        gsizes(1) = sum(part_ord_mpi(1:num_procs))
        gsizes(2) = 21 + varsMarmotant

        start_idx_part(1) = sum(part_ord_mpi(1:proc_rank + 1)) - part_ord_mpi(proc_rank + 1)!-MAX(1,nparticles)
        start_idx_part(2) = 0

        write (file_loc, '(A,I0,A)') 'particle_mpi_io', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if

        ! Writing down the total number of particles
        if (proc_rank == 0) then
            open (9, FILE=trim(file_loc), FORM='unformatted', STATUS='unknown')
            !write (9) gsizes(1), time_real, dt_next_inp, tot_step
            write (9) gsizes(1), mytime, dt
            close (9)
        end if

        call MPI_type_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, MPI_doUBLE_PRECISION, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        allocate (MPI_IO_DATA_particle(1:max(1, nparticles), 1:(21+varsMarmotant)))

        ! Open the file to write all flow variables
        write (file_loc, '(A,I0,A)') 'particle', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                           mpi_info_int, ifile, ierr)

        disp = 0d0

        call MPI_FILE_SET_VIEW(ifile, disp, MPI_doUBLE_PRECISION, view, &
                               'native', mpi_info_null, ierr)

        ! Cycle through list
        i = 1

        if (nparticles == 0) then
            MPI_IO_DATA_particle(1, 1:(21+varsMarmotant)) = 0d0
        else

            particle => particlesubList%List%next
            do while (associated(particle))

                if (particle_in_domain_physical(particle%data%x(1:3))) then

                    MPI_IO_DATA_particle(i, 1) = real(particle%data%id)
                    MPI_IO_DATA_particle(i, 2:4) = particle%data%x(1:3)
                    MPI_IO_DATA_particle(i, 5:7) = particle%data%xprev(1:3)
                    MPI_IO_DATA_particle(i, 8:10) = particle%data%u(1:3)
                    MPI_IO_DATA_particle(i, 11:12) = particle%data%y(1:2)
                    MPI_IO_DATA_particle(i, 13) = particle%data%R0
                    MPI_IO_DATA_particle(i, 14) = particle%data%Rmax
                    MPI_IO_DATA_particle(i, 15) = particle%data%Rmin
                    MPI_IO_DATA_particle(i, 16) = particle%data%dphidt
                    MPI_IO_DATA_particle(i, 17) = particle%data%p
                    MPI_IO_DATA_particle(i, 18) = particle%data%mv
                    MPI_IO_DATA_particle(i, 19) = particle%data%mg
                    MPI_IO_DATA_particle(i, 20) = particle%data%betaT
                    MPI_IO_DATA_particle(i, 21) = particle%data%betaC
                    ! varsMarmotant
                    MPI_IO_DATA_particle(i, 22) = particle%data%shell
                    MPI_IO_DATA_particle(i, 23) = particle%data%Rbuck
                    MPI_IO_DATA_particle(i, 24) = particle%data%Rrupt

                    i = i + 1

                end if

                particle => particle%next

            end do

        end if

        call MPI_FILE_write_ALL(ifile, MPI_IO_DATA_particle, (21+varsMarmotant)*max(1, nparticles), &
                                MPI_doUBLE_PRECISION, status, ierr)

        call MPI_FILE_CLOSE(ifile, ierr)

        deallocate (MPI_IO_DATA_particle)

#endif

    end subroutine s_write_restart_particles!_parallel

    subroutine s_calculate_particle_stats()

        type(particlenode), pointer :: particle

        ! Cycle through list
        particle => particlesubList%List%next
        do while (associated(particle))
            Rmax = max(Rmax, particle%data%y(1)/particle%data%R0)
            Rmin = min(Rmin, particle%data%y(1)/particle%data%R0)
            particle%data%Rmax = max(particle%data%Rmax, particle%data%y(1)/particle%data%R0)
            particle%data%Rmin = min(particle%data%Rmin, particle%data%y(1)/particle%data%R0)
            particle => particle%next
        end do

    end subroutine s_calculate_particle_stats

    subroutine s_write_particle_stats

        type(particlenode), pointer :: particle
        integer :: i
        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'stats_particles_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        open (13, FILE=trim(file_loc), FORM='formatted', position='rewind')
        write (13, *) 'proc_rank, particleID, x, y, z, Rmax, Rmin'

        particle => particlesubList%List%next
        do while (associated(particle))
            write (13, *) proc_rank, particle%data%id, (particle%data%x(i), i=1, 3), particle%data%Rmax, particle%data%Rmin
            particle => particle%next
        end do

        close (13)

    end subroutine s_write_particle_stats

    !> The purpose of this subroutine is to remove the particles information
          !!        from the sublist space
          !! @param particle     Particle data
          !! @param originlist   Particle list id: 1=from qbl, 2=SubList
    subroutine s_remove_particle(particle, originlist)

        type(particlenode), pointer :: sublistnode, particle, nodeaux, oldnode
        type(cellwb), pointer :: cellwbaux, oldcellnode
        integer, dimension(3) :: cell
        integer :: originlist

        if (originlist == 1) then
            sublistnode => particlesubList%List%next
            do while (.not. associated(sublistnode%data, particle%data))
                sublistnode => sublistnode%next
            end do
        else
            sublistnode => particle
        end if

        !! Removing element from the 3D list structure
        call s_get_cell(particle%data%tmp%s, cell)
        nodeaux => qbl%fp(cell(1), cell(2), cell(3))%List%next
        !! Removing node from qbl
801     if (associated(nodeaux%data, particle%data)) then
            !! Removing from 3D list
            oldnode => nodeaux
            if (associated(nodeaux%next)) then
                if (associated(nodeaux%prev)) nodeaux%prev%next => nodeaux%next
                nodeaux%next%prev => nodeaux%prev
            else
                if (associated(nodeaux%prev)) nullify (nodeaux%prev%next) !last element of the list
            end if
            if (associated(nodeaux, qbl%fp(cell(1), cell(2), cell(3))%List%next)) then !first element
                if (associated(nodeaux%next)) then
                    qbl%fp(cell(1), cell(2), cell(3))%List%next => nodeaux%next
                else
                    nullify (qbl%fp(cell(1), cell(2), cell(3))%List%next)
                end if
            end if
            deallocate (oldnode)
            qbl%fp(cell(1), cell(2), cell(3))%nb = qbl%fp(cell(1), cell(2), cell(3))%nb - 1

            if (qbl%fp(cell(1), cell(2), cell(3))%nb == 0) then
                cellwbaux => qbl%fp(cell(1), cell(2), cell(3))%cellpointer
            !!  Removing from cell list
                oldcellnode => cellwbaux
                if (associated(cellwbaux%next)) then
                    if (associated(cellwbaux%prev)) cellwbaux%prev%next => cellwbaux%next
                    cellwbaux%next%prev => cellwbaux%prev
                else
                    if (associated(cellwbaux%prev)) nullify (cellwbaux%prev%next) !last element of the list
                end if
                if (associated(cellwbaux, cellwbList%List%next)) then !first element
                    if (associated(cellwbaux%next)) then
                        cellwbList%List%next => cellwbaux%next
                    else
                        nullify (cellwbList%List%next)
                    end if
                end if
                cellwbList%nb = cellwbList%nb - 1
                deallocate (oldcellnode%data); deallocate (oldcellnode)
            end if
        else
            nodeaux => nodeaux%next
            goto 801
        end if

        !deallocating particle sublist node (if it is the first)
        !deallocating data
        if (associated(sublistnode%next)) then
            if (associated(sublistnode%prev)) then
                sublistnode%prev%next => sublistnode%next
                sublistnode%next%prev => sublistnode%prev
            else
                nullify (sublistnode%next%prev)
            end if
        else
            if (associated(sublistnode%prev)) nullify (sublistnode%prev%next) !required?
        end if

        if (associated(sublistnode%data, particlesubList%List%next%data)) then !first element
            if (associated(sublistnode%next)) then
                particlesubList%List%next => sublistnode%next
            else
                nullify (particlesubList%List%next)
            end if
        end if
        particlesubList%nb = particlesubList%nb - 1

        !next particle
        if (originlist == 1) then
            particle => qbl%fp(cell(1), cell(2), cell(3))%List%next
        else
            if (associated(particle%next)) then
                particle => particle%next
            else
                nullify (particle)
            end if
        end if

        deallocate (sublistnode%data)
        deallocate (sublistnode)

    end subroutine s_remove_particle

    subroutine s_deallocate_particles()

        integer :: i, j, k, imax
        type(particlenode), pointer :: particle

        if ((solverapproach == 2) .and. avgdensflag) then
            imax = 4
            if (clusterflag >= 4) imax = 10 !subgrid noise model
        else if ((solverapproach == 0) .and. avgdensflag) then
            imax = 3
        else
            imax = 2
        end if
        do i = 1, imax
            deallocate (q_particle(i)%sf)
        end do
        deallocate (q_particle)
        particle => particlesublist%list%next
        do while (associated(particle))
            call s_remove_particle(particle, 2)
        end do
        do k = -buff_size, p + buff_size
            do j = -buff_size, n + buff_size
                do i = -buff_size, m + buff_size
                    deallocate (qbl%fp(i, j, k)%list)
                end do
            end do
        end do
        deallocate (qbl%fp)
        deallocate (cellwblist%list)
        deallocate (cellwblist)

    end subroutine s_deallocate_particles

end module m_particles