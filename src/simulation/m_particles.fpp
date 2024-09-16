!>
!! @file m_particles.f90
!! @brief Contains module m_particles

#:include 'macros.fpp'

!> @brief This module is used to compute the Euler-Lagrangian sub-grid bubble dynamic variables
module m_particles

    ! Dependencies =============================================================

    use m_global_parameters     !< Definitions of the global parameters

    use m_derived_types         !< Definitions of the derived types

    use m_particles_types       !< Definitions of the derived particle types

    use m_rhs                   !< Right-hand-side (RHS) evaluation procedures

    use m_mpi_particles         !< Message passing interface (MPI) module for the particles

    use m_data_output           !< Run-time info & solution data output procedures

    use m_particles_output      !< Run-time info & solution data output procedures for particles

    use m_mpi_common
    
    ! ==========================================================================

    implicit none

    contains

    !> Initializes the lagrangian solver
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_initialize_lagrangian_solver(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0.d0)) :: dtoutput, tend
        integer :: i,imax

        ! Time-step iterator to the first time-step
        dt0 = dt
        !dtoutput = dt * t_step_save
        time_real = t_step_start * dt
        !tend = time_real + ( t_step_stop - t_step_start) * dt
        dt_next_inp = dt0

        ! Initializing particle modules 
        if((solverapproach.eq.2).and.avgdensflag) then
        ! comp 1: (1 - beta)
        ! comp 2: dbetadt
        ! comp 3 - imax : auxiliary variables
            imax=4
            if(clusterflag.ge.4) imax=10 !Subgrid noise model
            bubblesources=.true.
        else if ((solverapproach.eq.0).and.avgdensflag) then
            imax = 3
            bubblesources=.false.
        else
            imax=2
            bubblesources=.false.
        end if
        allocate(q_particle(1:imax))
        do i=1,imax
            if(p > 0) then
                dim=3
                allocate(q_particle(i)%sf(-buff_size : m+buff_size, &
                    -buff_size : n+buff_size, -buff_size : p+buff_size ))
            else
                dim=2
                allocate(q_particle(i)%sf(-buff_size : m+buff_size, &
                                       -buff_size : n+buff_size, 0 : 0 ))
            end if
        end do
        q_particle(1)%sf = 1.d0 !represents 1-beta
        do i=2,imax
            q_particle(i)%sf = 0.d0
        end do

        call s_read_input_particles(q_cons_vf, q_prim_vf)

        if(t_step_start.eq.0) then
            call s_write_data_files(q_cons_vf, q_prim_vf, t_step_start, q_particle(1)) !Parameter 'parallel_io' must be True 
            call write_restart_particles_parallel (t_step_start)
            if(avgdensflag) call write_void_evol (time_real)
        end if
        call s_populate_primitive_variables_buffers(q_cons_vf,q_particle=q_particle)

        ! Recovering dtnext
        if(t_step_start .ne. 0.0d0) then
            dtnext = dt_next_inp
        else if(dtmaxpart.gt.0.0d0) then
            dtnext = min(dt0, dtmaxpart)
        else
            dtnext = dt0
        end if


    end subroutine s_initialize_lagrangian_solver

    !> Non-dimensionalize inputs (NEED TO BE FIXED)
    subroutine s_particles_nondimensionalize_inputs

        ! Non-dimensionalizing inputs
        if (particleflag) then
            epsilonb = 1
            Runiv = Runiv*Tini/(csonref**2)
            pvap = pvap/(rholiqref*csonref**2)
            cpgas = cpgas*Tini/(csonref**2)
            cpvapor = cpvapor*Tini/(csonref**2)
            kgas = kgas*Tini/(Lref*csonref**3*rholiqref)
            kvapor = kvapor*Tini/(Lref*csonref**3*rholiqref)
            diffcoefvap = diffcoefvap/(Lref*csonref)
            sigmabubble = sigmabubble/(rholiqref*csonref**2*Lref)
            viscref = viscref/(rholiqref*csonref*Lref)
            if (.not.avgdensFlag) solverapproach = 0

        ! Default values
        else
            avgdensFlag = .false.
            particleoutFlag = .false.
            particlestatFlag = .false.
            RPflag = .false.
            clusterflag = dflt_int
            stillparticlesflag= .false.
            heatflag= dflt_int
            massflag= dflt_int
            csonref = dflt_real
            rholiqref = dflt_real
            Lref = dflt_real
            Tini = dflt_real
            Runiv = dflt_real
            gammagas = dflt_real
            gammavapor = dflt_real
            pvap = dflt_real
            cpgas = dflt_real
            cpvapor = dflt_real
            kgas = dflt_real
            kvapor = dflt_real
            MWgas = dflt_real
            MWvap = dflt_real
            diffcoefvap = dflt_real
            sigmabubble = dflt_real
            viscref = dflt_real
            RKeps = dflt_real
            ratiodt = dflt_int
            projectiontype = dflt_int
            smoothtype = dflt_int
            epsilonb= dflt_real
            coupledFlag = .false.
            solverapproach = 2
            correctpresFlag = .false.
            charwidth = dflt_real
            valmaxvoid = dflt_real
            dtmaxpart = dflt_real
            do_particles = .false.
            bubblesources = .false.
        end if

    end subroutine s_particles_nondimensionalize_inputs

    !> The purpose of this procedure is to read the input file with the particles' information
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_read_input_particles(q_cons_vp, q_prim_vf)

        use m_particles_output

        ! Conservative variables
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vp
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(particlenode), pointer           :: particle
        real(kind(0.d0)), dimension(8)   :: inputparticle
        real(kind(0.d0))                 :: qtime
        integer  :: i,j,k,l,nparticles
        logical  :: file_exist,indomain
        
        call s_populate_primitive_variables_buffers(q_cons_vp,q_particle=q_particle)

        if(model_eqns == 2 .and. (adv_alphan .neqv. .true.)) then        
            q_cons_vp(sys_size)%sf = 1d0
            do i = adv_idx%beg, adv_idx%end
                q_cons_vp(sys_size)%sf = &
                q_cons_vp(sys_size)%sf - &
                q_cons_vp(i)%sf
            end do
        end if

        ix%beg = -buff_size; ix%end = m + buff_size
        iy%beg = -buff_size; iy%end = n + buff_size
        if(p > 0) then
            iz%beg = -buff_size; iz%end = p + buff_size
        end if

        ! To get the pressure to be used in s_add_particle()
        call s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp%vf, ix, iy, iz, q_particle(1))

        call compute_cell_centers()
        nparticles = 0
        if(num_procs > 1) call initialize_particles_mpi()
        call initialize_particlelists
        if (t_step_start.eq.0) then
            inquire (file='input/particles.dat'   ,exist=file_exist    )
            if (file_exist) then
                open(unit=85,file='input/particles.dat',form='formatted')
                101 read(85 ,*,end=102)  (inputparticle(i), i=1,8)
                indomain = particle_in_domain(inputparticle(1:3))
                id = id + 1 
                if (indomain) then
                  nparticles = nparticles + 1
                  call s_add_particle(inputparticle, q_cons_vp, q_prim_vf)
                end if
                goto 101
                102 continue
            else
                stop "if you include particles, you have to initializate them in input/particles.dat"
            end if
        else
            call s_add_particle_restart (nparticles) !Parameter 'parallel_io' must be True 
        end if

        print *, " PARTICLES RUNNING, in proc", proc_rank, "number:", nparticles, "/",id
        call Create_sublists ()

        ! Apply density correction
        if (avgdensFlag) then
            particle  => particlesubList%List%next
            do while(Associated(particle))
                call transfertotmp (particle%data)
                particle => particle%next
            end do
            call s_smear_voidfraction(q_cons_vp)
            if ( solverapproach.eq.1 ) then
            !Definition averaged quantities
                q_cons_vp(E_idx)%sf(0:m,0:n,0:p) = q_cons_vp(E_idx)%sf(0:m,0:n,0:p) * q_particle(1)%sf(0:m,0:n,0:p)
                do i = 1, cont_idx%end 
                    q_cons_vp(i)%sf(0:m,0:n,0:p)   = q_cons_vp(i)%sf(0:m,0:n,0:p)     * q_particle(1)%sf(0:m,0:n,0:p)
                end do
                do i = mom_idx%beg, mom_idx%end
                    q_cons_vp(i)%sf(0:m,0:n,0:p)   = q_cons_vp(i)%sf(0:m,0:n,0:p)     * q_particle(1)%sf(0:m,0:n,0:p)
                end do
            end if
        end if

        qtime = 0.0d0
        if (particleoutFlag)  call write_particles (qtime)
        call s_populate_primitive_variables_buffers(q_cons_vp,q_particle=q_particle)
        if(model_eqns == 2 .and. (adv_alphan .neqv. .true.)) then
            q_cons_vp(sys_size)%sf = 1d0
            do i = adv_idx%beg, adv_idx%end
                q_cons_vp(sys_size)%sf = &
                q_cons_vp(sys_size)%sf - &
                q_cons_vp(i)%sf
            end do
        end if

        ix%beg = -buff_size; ix%end = m + buff_size
        iy%beg = -buff_size; iy%end = n + buff_size
        if(p > 0) then
            iz%beg = -buff_size; iz%end = p + buff_size
        end if

        call s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp%vf, &
                                                                    ix, iy, iz, q_particle(1))
  
    end subroutine s_read_input_particles

    !> The purpose of this procedure is to add information of the particles when starting fresh
        !! @param inputparticle Particle number
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_add_particle(inputparticle, q_cons_vp, q_prim_vf)

        type(particledata)    , pointer :: particleinfo
        type(particleListinfo), pointer :: particleListaux
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vp
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0.d0)), dimension(8),intent(IN)  :: inputparticle
        real(kind(0.d0))  :: pliq, volparticle, concvap, totalmass, kparticle, &
                            cpparticle,omegaN,PeG,PeT,rhol, cson, pcrit
        integer, dimension(3) :: cell
        real(kind(0.d0)), dimension(2) :: Re
        real(kind(0.d0)), dimension( num_fluids, num_fluids ) :: We

        allocate(particleinfo)
        particleinfo%id        = id
        particleinfo%x(:)      = inputparticle(1:3)
        particleinfo%xprev(:)  = inputparticle(1:3)
        particleinfo%u(:)      = inputparticle(4:6)
        particleinfo%y(1)      = inputparticle(7)  ! particle radius
        particleinfo%R0        = inputparticle(7)
        particleinfo%y(2)      = inputparticle(8)  ! interface velocity
        particleinfo%Rmax      = 1.
        particleinfo%Rmin      = 1.
        particleinfo%dphidt    = 0.0d0

        if(cyl_coord .and. p.eq.0) then
            particleinfo%x(2) = dsqrt(particleinfo%x(2)**2d0+particleinfo%x(3)**2d0)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            particleinfo%x(3) = ATAN2(inputparticle(3),inputparticle(2))
            particleinfo%xprev = particleinfo%x
        end if
        cell    = -buff_size
        call locate_cell ( particleinfo%x,  cell, particleinfo%tmp%s )
        pliq = Interpolate( particleinfo%tmp%s, q_prim_vf(E_idx))
        if(pliq<0) print *, "Negative pressure", proc_rank, &
            q_cons_vp(E_idx)%sf(cell(1),cell(2),cell(3)),q_prim_vf(E_idx)%sf(cell(1),cell(2),cell(3)), cell
        call get_mixture_variables(q_prim_vf, q_cons_vp, pliq, cell(1), cell(2), cell(3), rhol, cson, Re, We)

        ! Intial particle pressure
        particleinfo%p =  pliq + 2.0*sigmabubble/particleinfo%R0
        if(sigmabubble.NE.0.0d0) then
            pcrit = pvap - 4.0d0*sigmabubble/(3.d0*sqrt(3.0d0*particleinfo%p*particleinfo%R0**3/(2.0d0*sigmabubble)))
            pref  = particleinfo%p 
        else
            pcrit = 0.0d0
        end if

        particleinfo%equilibrium = .false.

        ! Initial particle mass
        volparticle = 4.0d0/3.0d0*pi*particleinfo%R0**3 ! volume
        particleinfo%mv = pvap*volparticle*MWvap/Runiv*DBLE(massflag) ! vapermass
        particleinfo%mg = (particleinfo%p-pvap*DBLE(massflag))*volparticle*MWgas/Runiv ! gasmass
        if (particleinfo%mg.LE.0.0d0) stop 'the initial mass of gas inside the bubble is negative. Check your initial conditions'
        totalmass = particleinfo%mg + particleinfo%mv ! totalmass

        ! Bubble natural frequency
        concvap = particleinfo%mv/(particleinfo%mv+particleinfo%mg)
        omegaN = (3.0d0*(particleinfo%p-pvap*real(massflag))+4.0d0*sigmabubble/particleinfo%R0)/rhol
        if (pvap*real(massflag).gt.particleinfo%p) then
            print *, 'Not allowed: bubble initially located in a region with pressure below the vapor pressure'
            print *, 'location:', particleinfo%x(1:3)
            stop
        end if
        omegaN = dsqrt(omegaN/particleinfo%R0**2)

        cpparticle = concvap*cpvapor    + (1.0d0 - concvap)*cpgas
        kparticle  = concvap*kvapor     + (1.0d0 - concvap)*kgas
        ! Mass and heat transfer coefficients (based on Preston 2007)

        PeT   = totalmass/volparticle * cpparticle * particleinfo%R0**2*omegaN/kparticle
        particleinfo%betaT = transfercoeff(PeT,1.0d0)*real(heatflag)
        PeG   = particleinfo%R0**2*omegaN/diffcoefvap
        particleinfo%betaC = transfercoeff(PeG,1.0d0)*real(massflag)

        ! terms to work out directly the heat flux in getfluxes
        particleinfo%betaT = particleinfo%betaT*kparticle

        if (particleinfo%mg.LE.0.0d0) stop 'PROBLEM WITH THE MASS OF THE particle, CHECK PARTICLES ARE INSIDE THE doMAIN'

        particleListaux => qbl%fp(cell(1),cell(2),cell(3))
        if (particleListaux%nb.eq.0) call addtocell_list ( cell )
        call transfertotmp (particleinfo) 
        call addparticletolist (particleinfo,particleListaux)

    end subroutine s_add_particle

    !> The purpose of this procedure is to add information of the particles from a restart point in parallel
        !! @param nparticles Particle number in the domain
    subroutine s_add_particle_restart (nparticles)

        character(len = len_trim(case_dir) + 2*name_len) :: t_step_dir
        character(len = len_trim(case_dir) + 3*name_len) :: file_loc
        logical :: dir_check
        integer :: id, nparticles
#ifdef MFC_MPI
        real(kind(0.d0)), dimension(20)   :: inputvals
        real(kind(0.d0)) :: id_real
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer :: view

        type(particledata)    , pointer :: particleinfo
        type(particlelistinfo), pointer :: particlelistaux

        integer, dimension(3)   :: cell
        logical                 :: indomain, particle_file, file_exist

        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer :: ifile, ireq, ierr, data_size, tot_data
        integer :: i
    
        write(file_loc, '(a,i0,a)') 'particle_mpi_io' , t_step_start, '.dat'
          file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
        inquire(file = trim(file_loc),exist = file_exist)
        
        if (file_exist) then
            if(proc_rank .eq. 0) then
                open(9, file = trim(file_loc), form = 'unformatted', status = 'unknown')
                read(9) tot_data, time_real, dt_next_inp, tot_step
                close(9)
            end if
        else
            print '(a)', trim(file_loc) // ' is missing. exiting ...'
            call s_mpi_abort
        end if
    
        call MPI_BCAST(tot_data, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(time_real, 1, MPI_doUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt_next_inp, 1, MPI_doUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        gsizes(1)=tot_data
        gsizes(2)=21
        lsizes(1)=tot_data
        lsizes(2)=21
        start_idx_part(1)=0
        start_idx_part(2)=0

        call MPI_type_CREATE_SUBARRAY(2,gsizes,lsizes,start_idx_part,&
                      MPI_ORDER_FORTRAN,MPI_doUBLE_PRECISION,view,ierr)
        call MPI_type_COMMIT(view,ierr)

        ! open the file to write all flow variables
        write(file_loc, '(a,i0,a)') 'particle' , t_step_start, '.dat'
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
        inquire(file = trim(file_loc),exist = particle_file)
  
        if (particle_file) then
            call MPI_FILE_open(MPI_COMM_WORLD,file_loc,MPI_MODE_RdoNLY, &
                                 mpi_info_int,ifile,ierr)
            disp = 0d0
            call MPI_FILE_SET_VIEW(ifile,disp,MPI_doUBLE_PRECISION,view, &
                                              'native',mpi_info_null,ierr)
            allocate(MPI_IO_DATA_particle(tot_data,1:21))
            call MPI_FILE_read_ALL(ifile,MPI_IO_DATA_particle,21*tot_data, &
                                          MPI_doUBLE_PRECISION,status,ierr)
            do i = 1, tot_data
                id = int(MPI_IO_DATA_particle(i,1))
                inputvals(1:20) = MPI_IO_DATA_particle(i,2:21)
                indomain = particle_in_domain(inputvals(1:3))
                if (indomain.and.(id.gt.0)) then
                    allocate(particleinfo)
                    nparticles = nparticles + 1
                    particleinfo%id        = id
                    particleinfo%x(1:3)    = inputvals(1:3) 
                    particleinfo%xprev(1:3)= inputvals(4:6) 
                    particleinfo%u(1:3)    = inputvals(7:9) 
                    particleinfo%y(1:2)    = inputvals(10:11) 
                    particleinfo%R0        = inputvals(12) 
                    particleinfo%Rmax      = inputvals(13) 
                    particleinfo%Rmin      = inputvals(14) 
                    particleinfo%dphidt    = inputvals(15) 
                    particleinfo%p         = inputvals(16) 
                    particleinfo%mv        = inputvals(17) 
                    particleinfo%mg        = inputvals(18) 
                    particleinfo%betaT     = inputvals(19) 
                    particleinfo%betaC     = inputvals(20) 
                    particleinfo%equilibrium = .false.
                    cell    = -buff_size
                    call locate_cell ( particleinfo%x,  cell, particleinfo%tmp%s )
                    particleListaux => qbl%fp(cell(1), cell(2), cell(3))
                    if (particleListaux%nb.eq.0) call addtocell_list ( cell )
                    call transfertotmp (particleinfo) 
                    call addparticletolist (particleinfo,particleListaux)
                end if
            end do
            deallocate(MPI_IO_DATA_particle)
        end if
        call MPI_FILE_CLOSE(ifile,ierr)
#endif

    end subroutine s_add_particle_restart

    !> Calculates the standard deviation of the particle being smeared in the Eulerian frame
        !! @param cell Cell where the particle is located
        !! @param kernel Kernel type or smootheing function
        !! @param volpart Volume of the particle
        !! @param stddsv Standard deviaton
    subroutine s_compute_stddsv(cell, kernel, volpart, stddsv)

      real(kind(0.d0))         :: chardist, volpart, rad, stddsv, charvol
      integer, dimension(3)    :: cell
      integer                  :: kernel

      call get_char_dist(cell,chardist)
      call get_char_vol(cell,charvol)
      kernel = smoothtype
    
      if (((volpart/charvol).gt.0.5d0*valmaxvoid).or.(smoothtype.eq.1)) then
        kernel = 1
        rad    = (3.0d0*volpart/(4.0d0*pi))**(1.0d0/3.0d0)
        stddsv = 1.0d0*epsilonb*MAX(chardist,rad)
      else
        stddsv = 0.0d0
      end if

    end subroutine s_compute_stddsv

    !>  The purpose of this subroutine is to smear the particles in the Eulerian frame
        !! @param q Smeared particles in the conservative variables
    subroutine s_smear_voidfraction (q)

        use m_kernel_functions

        type(particlenode),pointer       :: particle
        type(scalar_field), dimension(sys_size), intent(IN) :: q
        real(kind(0.d0))               :: volpart, totmass, stddsv, volpart2
        real(kind(0.d0)), dimension(3) :: nodecoord
        integer, dimension(3)          :: cell
        integer                        :: i, j, k, l, kernel
        integer, dimension(3,2)        :: rangecells

        q_particle(1)%sf = 0.0d0
        q_particle(2)%sf = 0.0d0
        nodecoord(3)     = 0

        rangecells(1:3,2) = -buff_size
        rangecells(1,1) = m
        rangecells(2,1) = n
        rangecells(3,1) = p

        if (projectiontype.eq.0) then   
            particle  => particlesubList%List%next
            do while(Associated(particle))
                volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
                cell = get_cell_from_s(particle%data%tmp%s)
                call update_rangecells (cell, rangecells)
                nodecoord(1) = particle%data%tmp%x(1)
                nodecoord(2) = particle%data%tmp%x(2)
                if (p > 0) nodecoord(3) = particle%data%tmp%x(3)
                call s_compute_stddsv(cell, kernel, volpart, stddsv)
                call smoothfunction ( q_particle(1), nodecoord, cell , volpart, kernel, stddsv)
                particle => particle%next
            end do
            subrange(1)%beg=max(rangecells(1,1), -buff_size)
            subrange(1)%end=min(rangecells(1,2),m+buff_size)
            subrange(2)%beg=max(rangecells(2,1), -buff_size)
            subrange(2)%end=min(rangecells(2,2),n+buff_size)
            subrange(3)%beg=max(rangecells(3,1), -buff_size)
            subrange(3)%end=min(rangecells(3,2),p+buff_size)
        else
            particle  => particlesubList%List%next
            do while(Associated(particle))
                volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
                call remeshdelta ( q_particle(2), particle%data, volpart, rangecells )     
                particle => particle%next
            end do
            if (bubblesources) q_particle(4)%sf = q_particle(2)%sf
            subrange(1)%beg=max(rangecells(1,1), -buff_size)
            subrange(1)%end=min(rangecells(1,2),m+buff_size)
            subrange(2)%beg=max(rangecells(2,1), -buff_size)
            subrange(2)%end=min(rangecells(2,2),n+buff_size)
            subrange(3)%beg=max(rangecells(3,1), -buff_size)
            subrange(3)%end=min(rangecells(3,2),p+buff_size)
            do k=0,p
                do j=0,n
                    do i=0,m
                        cell(1) = i 
                        cell(2) = j
                        cell(3) = k
                        nodecoord(1)=x_cc_lp(cell(1))
                        nodecoord(2)=y_cc_lp(cell(2))
                        if (p.gt.0) nodecoord(3)=z_cc_lp(cell(3))

                        if (q_particle(2)%sf(i,j,k).NE.0.0d0) then
                            call s_compute_stddsv(cell, kernel, q_particle(2)%sf(i,j,k), stddsv)
                            call smoothfunction ( q_particle(1), nodecoord, cell , q_particle(2)%sf(i,j,k), kernel, stddsv)
                        end if
                    end do
                end do
            end do

        end if

        !I store 1-beta (I should probably do it initializing
        !                q_particle(1)%fp=1 and using a negative volume)
        q_particle(1)%sf = 1. - q_particle(1)%sf
        subrange(1)%beg=max(subrange(1)%beg-CEILING(5*epsilonb), -buff_size*1.0d0)
        subrange(1)%end=min(subrange(1)%end+CEILING(5*epsilonb),(m+buff_size)*1.0d0)
        subrange(2)%beg=max(subrange(2)%beg-CEILING(5*epsilonb), -buff_size*1.0d0)
        subrange(2)%end=min(subrange(2)%end+CEILING(5*epsilonb),(n+buff_size)*1.0d0)
        if (p > 3) then
          subrange(3)%beg=max(subrange(3)%beg-CEILING(5*epsilonb), -buff_size*1.0d0)
          subrange(3)%end=min(subrange(3)%end+CEILING(5*epsilonb),(p+buff_size)*1.0d0)
        end if

        !=== dbetadt==========
        if (bubblesources) then
            q_particle(2)%sf = 0.0d0
            q_particle(3)%sf = 0.0d0
            if((clusterflag.GE.4)) then
                q_particle(5)%sf = 0.0d0
            end if
            if (projectiontype.eq.0) then
                particle  => particlesubList%List%next
                do while(Associated(particle))
                    volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
                    volpart2 = volpart
                    cell = get_cell_from_s(particle%data%tmp%s)
                    call s_compute_stddsv(cell, kernel, volpart, stddsv)
                    volpart = 4.0d0*pi*particle%data%tmp%y(1)**2*particle%data%tmp%y(2)
                    nodecoord(1) = particle%data%tmp%x(1)
                    nodecoord(2) = particle%data%tmp%x(2)
                    if (p > 0) nodecoord(3) = particle%data%tmp%x(3)
                    if(clusterflag.GE.4) call smoothfunction ( q_particle(5), nodecoord, cell , volpart, kernel, stddsv, volpart2)
                    call smoothfunction ( q_particle(2), nodecoord, cell , volpart, kernel, stddsv)
                    particle => particle%next
                end do
            else
                nodecoord(3)=0
                particle  => particlesubList%List%next
                do while(Associated(particle))
                  volpart = 4.0d0*pi*particle%data%tmp%y(1)**2*particle%data%tmp%y(2)
                  call remeshdelta ( q_particle(3), particle%data, volpart )     
                  particle => particle%next
                end do
                do i=0,m
                    do j=0,n
                        do k=0,p
                            cell(1) = i 
                            cell(2) = j
                            cell(3) = k
                            nodecoord(1)=x_cc_lp(cell(1))
                            nodecoord(2)=y_cc_lp(cell(2))
                            if (p > 0) nodecoord(3)=z_cc_lp(cell(3))
                            if (q_particle(3)%sf(i,j,k).NE.0.0d0) then
                                call s_compute_stddsv(cell, kernel, q_particle(4)%sf(i,j,k), stddsv)
                                call smoothfunction( q_particle(2), nodecoord, cell , q_particle(3)%sf(i,j,k), kernel, stddsv)
                            end if
                        end do
                    end do
                end do
            end if  
        end if
        subrange(1)%beg=max(subrange(1)%beg,0*1.0d0)
        subrange(1)%end=min(subrange(1)%end,m*1.0d0)
        subrange(2)%beg=max(subrange(2)%beg,0*1.0d0)
        subrange(2)%end=min(subrange(2)%end,n*1.0d0)
        if (p > 0) then
          subrange(3)%beg=max(subrange(3)%beg,0*1.0d0)
          subrange(3)%end=min(subrange(3)%end,p*1.0d0)
        end if

        ! Limiting void fraction given max value
        do k=0,p
            do j=0,n
                do i=0,m
                    q_particle(1)%sf(i,j,k) = max(q_particle(1)%sf(i,j,k),1.d0-valmaxvoid)
                end do
            end do
        end do

    end subroutine s_smear_voidfraction

    !>  The purpose of this subroutine is to add the particle source terms following the formulation of Kazuki and Colonius (2018)
        !! @param q Conservative variables
        !! @param dq Calculated change of conservative variables
        !! @param q_prim Conservative variables
    subroutine s_add_sources(q,dq,q_prim)

        type(scalar_field), dimension(sys_size), intent(IN) :: q
        type(scalar_field), dimension(sys_size), intent(IN) :: dq
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim
        integer :: i,j,k,l

        do k=0,p
            do j=0,n
                do i=0,m
                    if (q_particle(1)%sf(i,j,k).gt.(1.0d0-valmaxvoid)) then
                        do l=1,E_idx
                            if(clusterflag.GE.4) then
                                dq(l)%sf(i,j,k) = dq(l)%sf(i,j,k) + q(l)%sf(i,j,k)*(q_particle(2)%sf(i,j,k)+q_particle(5)%sf(i,j,k))
                            else
                                dq(l)%sf(i,j,k) = dq(l)%sf(i,j,k) + q(l)%sf(i,j,k)/q_particle(1)%sf(i,j,k)*q_particle(2)%sf(i,j,k)
                            end if
                        end do
                    end if
                end do 
            end do 
        end do

        do l=1,DIM
            call gradient_dir(q_prim(E_idx), q_particle(3),l)

            do k=0,p
                do j=0,n
                    do i=0,m
                        if (q_particle(1)%sf(i,j,k).gt.(1.0d0-valmaxvoid)) then
                            dq(mom_idx%beg+l-1)%sf(i,j,k) = dq(mom_idx%beg+l-1)%sf(i,j,k) - (1.0d0-q_particle(1)%sf(i,j,k))/q_particle(1)%sf(i,j,k)*q_particle(3)%sf(i,j,k)
                        end if
                    end do
                end do
            end do

            !source in energy
            q_particle(3)%sf = q_prim(E_idx)%sf * q_prim(mom_idx%beg+l-1)%sf
            call gradient_dir(q_particle(3), q_particle(4),l)

            do k=0,p
                do j=0,n
                    do i=0,m
                        if (q_particle(1)%sf(i,j,k).gt.(1.0d0-valmaxvoid)) then
                            dq(E_idx)%sf(i,j,k) = dq(E_idx)%sf(i,j,k) - q_particle(4)%sf(i,j,k)*(1.0d0-q_particle(1)%sf(i,j,k))/q_particle(1)%sf(i,j,k)
                        end if
                    end do
                end do
            end do
        end do

    end subroutine s_add_sources

    !>  Contains the bubble dynamics subroutines
        !! @param qtime Current time from the adaptative Runge-Kutta time stepper
        !! @param step Current time step in the adaptative stepper
        !! @param q Conservative variables
        !! @param t_step Current global time step
        !! @param q_prim Primitive variables
        !! @param dq Calculated change of conservative variables
        !! @param largestep Logical variable to determine if the adaptative time step is too large
    subroutine s_RK_particle_dynamics (qtime,step,q,t_step,q_prim,dq,largestep)

        type(particlenode),pointer                  :: particle
        type(scalar_field), dimension(sys_size)   :: q
        type(scalar_field), dimension(sys_size)   :: q_prim
        type(scalar_field), dimension(sys_size), OPTIONAL :: dq
        real(kind(0.d0)),dimension(5)             :: intvalues
        real(kind(0.d0)),dimension(3)             :: totalforce, DupDt
        real(kind(0.d0))                          :: gammaparticle,vaporflux,heatflux,qtime
        integer,dimension(3)                      :: cell
        integer                                   :: i,j,k,l
        integer,intent(IN)                        :: step, t_step
        logical,OPTIONAL                          :: largestep
        real(kind(0.d0)) :: time_avg

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================


        if (avgdensFlag) then
            do i=1, sys_size
                dq(i)%sf=0.0d0
            end do
        end if
        
        !update vbles
        if (avgdensFlag) call s_smear_voidfraction (q)
        do i = 1, cont_idx%end
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        q_prim(i)%sf(j, k, l) = q(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        do i = adv_idx%beg, sys_size
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        q_prim(i)%sf(j, k, l) = q(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        
        if (coupledFlag) then
            call s_compute_rhs(q, q_prim, dq, t_step=t_step,time_avg=time_avg, qtime=qtime)
            if (num_procs > 1) then
                call bcst_largestep(largestep)
            end if
        end if

        if (avgdensFlag) then
            call s_populate_primitive_variables_buffers(q,q_particle=q_particle)
            if(model_eqns == 2 .and. (adv_alphan .neqv. .true.)) then
                q(sys_size)%sf = 1d0
                do i = adv_idx%beg, adv_idx%end
                    q(sys_size)%sf = &
                    q(sys_size)%sf - &
                    q(i)%sf
                end do
            end if
            ix%beg = -buff_size; ix%end = m + buff_size
            iy%beg = -buff_size; iy%end = n + buff_size
            if(p > 0) iz%beg = -buff_size; iz%end = p + buff_size
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
            if ((clusterflag.gt.0).or.(correctpresFlag)) call potentials (q_prim,q,qtime)

        end if

        if (bubblesources) call s_add_sources(q, dq, q_prim)

        particle  => particlesubList%List%next
        do while(Associated(particle))
            if (.NOT.particle%data%equilibrium) then
                call s_compute_interface_fluxes(particle%data, vaporflux, heatflux, gammaparticle)
                particle%data%dbdt(step)%dpbdt = deriv_gaspressure(particle%data, vaporflux, heatflux, gammaparticle )
                particle%data%dbdt(step)%dmvdt = 4.0d0*pi*particle%data%tmp%y(1)**2*vaporflux 
            else
                particle%data%dbdt(step)%dpbdt = 0.0d0
                particle%data%dbdt(step)%dmvdt = 0.0d0
            end if
            particle%data%dbdt(step)%dxdt(:) = 0.0d0
            particle%data%dbdt(step)%dudt(:) = 0.0d0
            particle => particle%next
        end do

        !Radial motion
        if (RPflag) then
            particle  => particlesubList%List%next
            do while(Associated(particle))
                if (.NOT.particle%data%equilibrium) then
                    call s_compute_RP ( particle%data, step, q, q_prim, qtime)
                else
                    particle%data%dbdt(step)%dydt(2) = 0.
                end if
                particle%data%dbdt(step)%dydt(1) = particle%data%tmp%y(2) !the derivative of the radius is the velocity
                particle => particle%next
            end do
        else
            particle  => particlesubList%List%next
            do while(Associated(particle))
                particle%data%dbdt(step)%dydt(2) = 0.0d0
                particle%data%dbdt(step)%dydt(1) = 0.0d0
                particle => particle%next
            end do
        end if

    end subroutine s_RK_particle_dynamics

    !>  This subroutine solves the Rayleight--Plesset equation
        !! @param bubbletmp Current particle information
        !! @param step Current time step in the adaptative stepper
        !! @param q Conservative variables
        !! @param q_prim Primitive variables
        !! @param qtime Current time from the adaptative Runge-Kutte time stepper
    subroutine s_compute_RP (bubbletmp, step, q, q_prim, qtime)

        type(scalar_field), dimension(sys_size)  :: q
        type(scalar_field), dimension(sys_size)  :: q_prim
        type(particledata)                 :: bubbletmp
        type(particlederivative)           :: dbdt
        real(kind(0.d0))                   :: pliqint, pbubble, deltaP, pinf, termI, aux1, &
                                              aux2, velint, rhol, cson
        integer, dimension(3)              :: cell
        integer                            :: step
        real(kind(0d0)), dimension(2)      :: Re
        real(kind(0d0)), dimension( num_fluids, num_fluids ) :: We
        real(kind(0d0)), OPTIONAL          :: qtime
        real(kind(0d0)) :: temp
            
        pbubble = min(bubbletmp%tmp%p, 1.d0) ! pres in the bubble
        pliqint = pressureliq_int(pbubble, bubbletmp%tmp%y(1), bubbletmp%tmp%y(2)) ! pres outside the bubble on the surface
        cell = get_cell_from_s(bubbletmp%tmp%s)
        pinf = get_pinf1(bubbletmp%tmp%s, q_prim(e_idx), 1, aux1, aux2) ! getting p_inf
        call get_mixture_variables(q_prim, q, pinf, cell(1), cell(2), cell(3), rhol, cson, re, we, q_particle(1))
        deltaP = pliqint - pinf
        termI = 0.0d0
        velint = bubbletmp%tmp%y(2) - bubbletmp%dbdt(step)%dmvdt/(4.0d0*pi*bubbletmp%tmp%y(1)**2*rhol)
        bubbletmp%dbdt(step)%dydt(2) =  ((1.0d0+velint/cson)*deltaP/rhol + termI &
                                        + bubbletmp%dbdt(step)%dpbdt*bubbletmp%tmp%y(1)/rhol/cson  &
                                        - velint**2*3.0d0/2.0d0*(1.0d0-velint/3.0d0/cson))         &
                                        / (bubbletmp%tmp%y(1)*(1.0d0-velint/cson))

    end subroutine s_compute_RP

    !>  This subroutine computes the fluxes at the bubbles' interface
        !! @param bubbletmp Current particle information
        !! @param vaporflux Mass flux
        !! @param heatflux Heat flux
        !! @param gammabubble Specific heat of the vapor-gas mixture in the bubble
    subroutine s_compute_interface_fluxes(bubbletmp, vaporflux, heatflux, gammabubble)

        real(kind(0.d0)) :: vaporflux, heatflux, concvapint,bubbleTemp,volbubble,kbubble,&
                            avgconc, Rmixt,gammabubble,rhogas
        type(particledata) :: bubbletmp

        if (massflag.eq.0) then
            concvapint = 0.d0
        else
            concvapint   = MWgas/MWvap*(bubbletmp%tmp%p/pvap-1.0d0)
            concvapint   = 1.0d0/(1.0d0+concvapint)
        end if

        bubbleTemp   = (bubbletmp%mg/MWgas + bubbletmp%tmp%mv/MWvap)*Runiv
        volbubble    = 4.0d0/3.0d0*pi*bubbletmp%tmp%y(1)**3
        bubbleTemp   = bubbletmp%tmp%p*volbubble/bubbleTemp

        gammabubble  = concvapint*gammavapor + (1.0d0 - concvapint)*gammagas !needed later (deriv_gaspressure)
        heatflux     = -(gammabubble-1.0d0)/gammabubble*bubbletmp%betaT*(bubbleTemp-1.0d0)/bubbletmp%tmp%y(1)

        avgconc      = bubbletmp%tmp%mv/(bubbletmp%mg + bubbletmp%tmp%mv)
        Rmixt        = (concvapint/MWvap + (1.0d0 - concvapint)/MWgas)*Runiv

        concvapint   = min(concvapint, 0.99d0)
        vaporflux    = (1.0d0- concvapint)*bubbletmp%tmp%y(1)
        rhogas       = (bubbletmp%mg+bubbletmp%tmp%mv)/(4.0d0/3.0d0*pi*bubbletmp%tmp%y(1)**3)
        vaporflux    = -diffcoefvap*bubbletmp%betaC*(avgconc -concvapint)*rhogas/vaporflux

    end subroutine s_compute_interface_fluxes

    !>  This subroutine updates the particle variables
        !! @param dt Time step
        !! @param RKstep Current time step in the adaptative stepper
        !! @param RKcoef Runge-Kutta coefficient
        !! @param largestep Logical variable to determine if the adaptative time step is too large
        !! @param q Conservative variables
        !! @param dq Calculated change of conservative variables
        !! @param q_prim Primitive variables
    subroutine s_update_particle(dt,RKstep,RKcoef, largestep,q,dq, q_prim )

        type(particlenode),pointer                :: particle
        real(kind(0.d0))                          :: dt
        real(kind(0.d0)),dimension(6),INTENT (IN) :: RKcoef
        integer, INTENT (IN)                      :: RKstep
        integer, dimension(3)                     :: oldcell, newcell
        integer                                   :: i,j
        logical                                   :: largestep,change,indomain
        type(vector_field), dimension(:), OPTIONAL :: q
        type(vector_field), dimension(:), OPTIONAL :: dq
        type(scalar_field), dimension(:), OPTIONAL :: q_prim
        type (cellwb), pointer        :: cellwbaux
        type(particlenode),pointer    :: nodeaux
        integer,dimension(3)          :: cell

        particle  => particlesubList%List%next

        do while(Associated(particle))

            oldcell =  get_cell_from_s(particle%data%tmp%s) 
            call transfertotmp (particle%data)

            do i=1,RKstep
            particle%data%tmp%y(1:2) = particle%data%tmp%y(1:2) + dt*RKcoef(i)*particle%data%dbdt(i)%dydt(1:2)
            particle%data%tmp%x(1:3) = particle%data%tmp%x(1:3) + dt*RKcoef(i)*particle%data%dbdt(i)%dxdt(1:3)
            particle%data%tmp%u(1:3) = particle%data%tmp%u(1:3) + dt*RKcoef(i)*particle%data%dbdt(i)%dudt(1:3)
            particle%data%tmp%p      = particle%data%tmp%p      + dt*RKcoef(i)*particle%data%dbdt(i)%dpbdt
            particle%data%tmp%mv     = particle%data%tmp%mv     + dt*RKcoef(i)*particle%data%dbdt(i)%dmvdt
            end do

            if ((particle%data%tmp%y(1).LE.0.0d0).or.(particle%data%tmp%x(1).NE.particle%data%tmp%x(1))) then
                if (dt.LT.2.d-15) then
                    print *, 'warning large step',dt, particle%data%id
                    call Remove_particle (particle,1)
                    call Remove_particle (particle,2)
                    goto 710
                end if
                if ((particle%data%tmp%x(1).NE.particle%data%tmp%x(1))) goto 711
            end if

            indomain = particle_in_domain(particle%data%tmp%x)

            if (.NOT.indomain) then 
            print *, 'not in domain', particle%data%id,particle%data%tmp%x(1),particle%data%xprev(1),x_cb(-buff_size-1),x_cb(n+buff_size)
            call Remove_particle (particle,2)
            goto 710
            end if

            if (particle%data%equilibrium) call s_equilibrium_state ( particle%data ,largestep, q_prim(E_idx) )

            particle => particle%next
            710 continue
        end do

711     if (num_procs > 1) then
            call bcst_largestep(largestep)
        end if

        !   !update fluid variables
        if (PRESENT(q)) then
            do i = 1, sys_size
                q(2)%vf(i)%sf(0:m,0:n,0:p) = q(1)%vf(i)%sf(0:m,0:n,0:p) 
                do j=1,RKstep
                    q(2)%vf(i)%sf(0:m,0:n,0:p) = q(2)%vf(i)%sf(0:m,0:n,0:p) & 
                                                            + dt*RKcoef(j)*dq(j)%vf(i)%sf(0:m,0:n,0:p) 
                end do
            end do
        end if

    end subroutine s_update_particle

    !>  This subroutine calculates the equilibrium state of the lagrangian bubbles
        !! @param particle Variables of the particle
        !! @param largestep Logical variable to determine if the adaptative time step is too large
        !! @param pres Pressure surrounding the bubble
    subroutine s_equilibrium_state ( particle, largestep, pres )

        type(particledata)     :: particle
        real(kind(0.d0))       :: req(1)
        real(kind(0.d0))       :: pinf, aux1, aux2
        integer, dimension(3)  :: cell
        logical                :: largestep,cond
        type(scalar_field)     :: pres

        cell = get_cell_from_s(particle%tmp%s)
        pinf = Interpolate( particle%tmp%s, pres )
        
        req(1) = particle%tmp%y(1)
        cond = .false.
        call get_equilibrium_radius (10000, req, 1.0d-10,1.0d-10, particle, pinf,cond)

        particle%tmp%y(1) = req(1)
        if ((req(1).LE.0.0d0).or.cond) then
            print *, 'released', particle%id
            particle%equilibrium = .false.
            RETURN
        end if
        particle%tmp%p    = pinf - pvap + 2.0d0*sigmabubble/req(1)
        particle%tmp%mv   = pvap*4.0d0/3.0d0*pi*req(1)**3*MWvap/Runiv

    end subroutine s_equilibrium_state

    !>  This subroutine calculates the maximum error of the RK step
        !! @param timetmp Current time
        !! @param dt Time discretization
        !! @param RKcoef Runge-Kutta coefficient
        !! @param errmax Maximum error
        !! @param largestep Logical variable to determine if the adaptative time step is too large
        !! @param t_step Time step
        !! @param q Conservative variables
        !! @param q_prim Primitive variables
        !! @param dq Calculated change of conservative variables
    subroutine s_calculate_RKerror (timetmp, dt, RKcoef, errmax, largestep, t_step, q, q_prim, dq)

        type(particlenode),pointer                  :: particle
        real(kind(0.d0)),INTENT (IN)              :: dt,timetmp
        real(kind(0.d0)),dimension(6),INTENT (IN) :: RKcoef
        real(kind(0.d0))                          :: errmax,erraux,errb
        type(vector_field), dimension(:), OPTIONAL :: q
        type(scalar_field), dimension(:), OPTIONAL :: q_prim
        type(vector_field), dimension(:), OPTIONAL :: dq
        integer :: i,j,k,l,l1,nb
        logical :: largestep
        integer, intent(IN) :: t_step
        real(kind(0.d0)) :: time_avg

        errmax = 0.0d0
        erraux = 0.0d0

        particle  => particlesubList%List%next
        do while(Associated(particle))
            errb   = 0.0d0
            if (.NOT.particle%data%equilibrium) then
                !particle radius error
                do i=1,6
                    erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dydt(1)
                end do
                errb=max(errb,abs(erraux)*dt/particle%data%R0)

                !interface velocity error
                erraux = 0.0d0
                do i=1,6
                    erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dydt(2)
                end do
                errb=max(errb,abs(erraux)*dt)

                !particle velocity error
                do j=1,3
                    erraux = 0.0d0
                    do i=1,6
                        erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dxdt(j)
                    end do
                    errb=max(errb,abs(erraux)*dt/(abs(particle%data%tmp%u(j))+1.0d-4))
                end do
            end if
            errmax=max(errmax,errb)
            particle => particle%next
        end do

        largestep = .false.
        if (PRESENT(q)) then 
            do l1=1, cont_idx%end 
                do k=0,p
                    do j=0,n
                        do i=0,m
                            erraux = q(1)%vf(l1)%sf(i,j,k)
                            do l=1,6
                                erraux = erraux + dt*RKcoef(l)*dq(l)%vf(l1)%sf(i,j,k)
                            end do
                            erraux = max(errmax,erraux)
                        end do
                    end do
                end do
            end do
            do l1=mom_idx%beg, mom_idx%beg+DIM-1
                do k=0,p
                    do j=0,n
                        do i=0,m
                            erraux = q(1)%vf(l1)%sf(i,j,k)
                            do l=1,6
                                erraux = erraux + dt*RKcoef(l)*dq(l)%vf(l1)%sf(i,j,k)
                            end do
                            erraux = max(errmax,erraux)
                        end do
                    end do
                end do
            end do
            call s_compute_rhs(q(2)%vf,q_prim,dq(6)%vf,t_step=t_step,time_avg=time_avg, qtime=timetmp)
            if (num_procs > 1) call bcst_largestep(largestep)
        end if

    end subroutine s_calculate_RKerror
 
    !>  This subroutine updates the conservative fields after performint the adaptative Runge-Kutta time stepper
        !! @param q Conservative variables
        !! @param update_fields Flag to update the fields
        !! @param q_prim Primitive variables
    subroutine s_update_RK (q, update_fields, q_prim)

        type(particlenode), pointer          :: particle
        type(vector_field), dimension(:), OPTIONAL :: q
        type(scalar_field), dimension(:), OPTIONAL :: q_prim
        integer                              :: i
        real(kind(0.d0))                     :: pinf, pcrit
        logical                              :: update_fields, release_part

        particle  => particlesubList%List%next
        release_part = .false.
        do while(Associated(particle)) 
            if (.NOT.particle%data%equilibrium) then
                pinf = get_pinf1(particle%data%tmp%s, q_prim(E_idx), 1)
                pcrit = pvap - &
                        4.0d0*sigmabubble/(3.d0*sqrt(3.0d0*(pref+.0*sigmabubble/particle%data%R0)*particle%data%R0**3/(2.0d0*sigmabubble))) 
                pcrit = min(pcrit,-pref)
                if (ABS((pcrit-pinf)/pcrit).LT.0.5d0) release_part = .true.
            end if
            call transfertodata (particle%data)
            particle => particle%next
        end do

        ! releasing particles
        particle  => particlesubList%List%next
        if (release_part) then
            do while(Associated(particle))
                particle%data%equilibrium = .false.
                particle => particle%next
            end do 
        end if

        if ((coupledflag.or.bubblesources).and.(update_fields)) then
            do i=1,sys_size; q(1)%vf(i)%sf = q(2)%vf(i)%sf ;  end do;
        end if

        if (avgdensFlag) call s_smear_voidfraction(q(1)%vf)

    end subroutine s_update_RK

    subroutine potentials (q_prim,q,qtime)

        type(scalar_field), dimension(sys_size)   :: q_prim
        type(scalar_field), dimension(sys_size), OPTIONAL   :: q
        type (cellwb), pointer            :: cellwbaux
        integer, dimension(3)             :: cell 
        integer                           :: info,i
        real(kind(0.d0)),OPTIONAL         :: qtime

        cellwbaux => cellwbList%List%next
        do while(Associated(cellwbaux))
            cell = cellwbaux%data%coord
            if(PRESENT(q)) then
                call solve_cell (cell,q_prim,q,qtime)
            else
                call solve_cell (cell,q_prim)
            end if
            cellwbaux => cellwbaux%next
        end do

    end subroutine potentials

    subroutine solve_cell (cell,q_prim,q,qtime)

        type(scalar_field), dimension(sys_size)   :: q_prim
        type(scalar_field), dimension(sys_size), OPTIONAL   :: q
        integer, dimension(3)            :: cell
        real(kind(0.d0))                 :: preterm1, term2, paux, pint, Romega, term1_fac, Rb
        real(kind(0.d0)), dimension(3)   :: scoord
        type(particlenode), pointer      :: particle
        real(kind(0.d0)), OPTIONAL       :: qtime
        real(kind(0.d0))                 :: rhol,cson
        real(kind(0d0)), dimension(2)              :: Re
        real(kind(0d0)), dimension( num_fluids, num_fluids ) :: We
    
        scoord(:) = cell(:) + 0.5d0
        paux = get_pinf1(scoord, q_prim(E_idx),2,preterm1,term2,Romega)
        particle => qbl%fp(cell(1),cell(2),cell(3))%List%next
    
        do while(Associated(particle))
            pint = pressureliq_int(particle%data%tmp%p,particle%data%tmp%y(1),particle%data%tmp%y(2)) + 0.5d0*particle%data%tmp%y(2)**2
            if(clusterflag.eq.2) then
                particle%Data%dphidt = (paux - pint) + term2
                ! Accouting for the potential induced by the bubble averaged over the control volume
                ! Note that this is based on the incompressible flow assumption near the bubble.
                Rb = particle%data%tmp%y(1)
                term1_fac=3.0d0/2.0d0*(Rb*(Romega**2d0-Rb**2d0))/(Romega**3d0-Rb**3d0)
                particle%Data%dphidt = particle%Data%dphidt/(1-term1_fac)
            end if
            particle => particle%next
        end do  

    end subroutine solve_cell

    function transfercoeff (Pe,omegaN) 

        real(kind(0.d0)) :: transfercoeff , Pe, omegaN
        complex          :: transferfunc, auxc

        transferfunc = csqrt(cmplx(0.0d0,Pe*omegaN))
        auxc = (cexp(-cmplx(2.0d0,0.0d0)*transferfunc)+1.0d0)/(-cexp(-cmplx(2.0d0,0.0d0)*transferfunc)+1.0d0)
        transferfunc = transferfunc*auxc-1.0d0
        transferfunc = 1.0d0/transferfunc - 3.0d0/cmplx(0.0d0,Pe*omegaN)
        transfercoeff = real(1.0d0/transferfunc)

    end function transfercoeff

    function pressureliq_int( pbubble, radius, bubblevel )

        real(kind(0.d0)) :: pressureliq_int, radius,bubblevel, pbubble

        pressureliq_int=  pbubble -  2.0d0*sigmabubble/radius - 4.*viscref*bubblevel/radius

    end function pressureliq_int

    function get_pinf1(scoord, pres, ptype, preterm1, term2, Romega)

        type(particlenode),pointer         :: bubble
        type(scalar_field)                 :: pres
        real(kind(0.d0)), dimension(3)     :: distance,center
        real(kind(0.d0))                   :: get_pinf1, jac, dij, dpotjdt, dc, vol, aux,&
                                            volgas, term1, Rbeq, denom, Rmax, stddsv,&
                                            charvol, charpres, charvol2, charpres2
        real(kind(0.d0)), dimension(3)     :: scoord
        real(kind(0.d0)), OPTIONAL         :: preterm1, term2, Romega
        integer, dimension(3)              :: cell, cellaux
        integer, dimension(3)              :: epsilonbaux
        integer                            :: ptype !1=p at infinity, 2= averaged P at the bubble location
        integer                            :: i, j, k, dir
        logical                            :: celloutside

        get_pinf1 = 0.0d0
        cell(:) = get_cell_from_s (scoord)
        
        if((clusterflag.eq.0)) then
        !getting p_cell in terms of only the current cell by interpolation

            call get_char_vol (cell, vol)
            bubble  => qbl%fp(cell(1),cell(2),cell(3))%List%next
            Rmax=0.0d0
            do while(Associated(bubble))
                Rmax  = Rmax+bubble%data%tmp%y(1)**3
                bubble  => bubble%next
            end do
        
            ! Surrogate bubble radius
            Rmax = Rmax**(1.0d0/3.0d0)
        
            ! Getting the cell volulme as Omega
            call get_char_dist(cell,stddsv)
            
            !p_cell (interpolated)
            get_pinf1 = Interpolate( scoord, pres )
        
            !R_Omega
            dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)
    
        else if(clusterflag .eq. 1) then
        ! Just making the characteristic volume 3^3 times bigger than the cell volume
        ! Not recomended to use for production

            ! Getting Omega (not Omega_L)
            call get_char_vol (cell, vol)
            vol=vol*2.7d1 !Just multiplying by 3^3=27
            bubble  => qbl%fp(cell(1),cell(2),cell(3))%List%next
            Rmax=0.0d0
            do while(Associated(bubble))
                Rmax  = Rmax+bubble%data%tmp%y(1)**3
                bubble  => bubble%next
            end do
        
            ! Surrogate bubble radius
            Rmax = Rmax**(1.0d0/3.0d0)
        
            ! Getting the cell volulme as Omega
            call get_char_dist(cell,stddsv)
            
            !p_cell to get interpolate
            get_pinf1 = Interpolate( scoord, pres )

            !R_Omega
            dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)
    
        else if(clusterflag .GE. 2 ) then
        ! Bubble dynamic closure from Kazuki and Colonius (2018)

            ! Range of cells included in Omega
            if(smoothtype.eq.1) then
                epsilonbaux(:) = 3
            else if(projectiontype.eq.1) then
                epsilonbaux(:) = 3
            end if
        
            charvol = 0.d0
            charpres = 0.d0
            charvol2 = 0.d0
            charpres2 = 0.d0
            vol = 0.d0
            if (DIM.eq.3) then
                k = -epsilonbaux(3)
            else
                k = 0
            end if
            i=-epsilonbaux(1);j=-epsilonbaux(2)
        
    3001    if ((i.LE.epsilonbaux(1)).and.(j.LE.epsilonbaux(2))) then
                celloutside = .false.
                cellaux(1) = cell(1) + i
                cellaux(2) = cell(2) + j
                cellaux(3) = cell(3) + k

                !Check ghost part in x-direction
                if (cellaux(1).LT.-buff_size) then
                    celloutside = .true.
                    i = i+1
                end if

                !Check ghost part in y-direction
                if (cellaux(2).LT.-buff_size) then
                    celloutside = .true.
                    j = j+1
                end if
                if(cyl_coord.and.(DIM.ne.3)) then
                    if (y_cc_lp(cellaux(2)).LT.0d0) then
                        celloutside = .true.
                        j = j+1
                    end if
                end if
        
                !Check ghost part in z-direction
                if(DIM.eq.3) then
                    if (cellaux(3).LT.-buff_size) then
                        celloutside = .true.
                        k = k+1
                    end if
                end if

                if (cellaux(1).gt.m+buff_size) celloutside =.true.
                if (cellaux(2).gt.n+buff_size) celloutside =.true.
                if (cellaux(3).gt.p+buff_size) celloutside =.true. 
                    
                if (.not. celloutside) then
                    call get_char_vol(cellaux,vol)
                    charvol  = charvol + vol
                    charpres = charpres + pres%sf(cellaux(1),cellaux(2),cellaux(3)) * vol
                    charvol2  = charvol2 + vol*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                    charpres2 = charpres2 + pres%sf(cellaux(1),cellaux(2),cellaux(3)) &
                            * vol*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                end if

                if (j.LT.epsilonbaux(2)) then
                    j = j+1
                    goto 3001
                end if
        
    3002        j=-epsilonbaux(2)
                i = i+1
                goto 3001
    
            end if
        
    3003    if ((DIM.eq.3).and.(k.LT.epsilonbaux(3))) then
                k = k+1
                i=-epsilonbaux(1);j=-epsilonbaux(2)
                goto 3001
            end if
            get_pinf1 = charpres2/charvol2
            vol=charvol
            dc = (3.0d0*abs(vol)/(4.0d0*pi))**(1.0d0/3.0d0)  !positive volume
        else

            print '(A)', 'Check cluterflag. Exiting ...'
            call s_mpi_abort()

        end if

        if (correctpresFlag.and.PRESENT(preterm1)) then
            dpotjdt = 0.0d0 !potential derivative contribution from other bubbles
            volgas  = 0.0d0
            term1 = 0.0d0
            term2 = 0.0d0
            denom = 0.0d0
            bubble  => qbl%fp(cell(1),cell(2),cell(3))%List%next

            do while(Associated(bubble))
                volgas  = volgas + bubble%data%tmp%y(1)**3 !surrogate bubble volume
                denom   = denom + bubble%data%tmp%y(1)**2
                term1   = term1 + bubble%Data%dphidt*bubble%data%tmp%y(1)**2
                term2   = term2 + bubble%data%tmp%y(2)*bubble%data%tmp%y(1)**2
                bubble  => bubble%next
            end do

            Rbeq = volgas**(1.0d0/3.0d0) !surrogate bubble radius
            aux = dc**3 - Rbeq**3
            term2 = term2/denom
            term2 = 3.0d0/2.0d0*term2**2*Rbeq**3*(1.0d0-Rbeq/dc)/aux
            preterm1  = 3.0d0/2.0d0*Rbeq*(dc**2 - Rbeq**2)/(aux*denom)

            !Control volume radius
            if(PRESENT(Romega)) Romega = dc

            ! Getting p_inf
            if (ptype.eq.1) then
                get_pinf1 =  get_pinf1 + preterm1*term1 + term2
            end if
    
        end if

        ! Test if pinf is valid
        if ((get_pinf1>10) .or. (get_pinf1<-10) .or. (ieee_is_nan(get_pinf1))) then
            print*, 'Extreme pinf value of', get_pinf1,'charvol2', charvol2,&
                    'charpres2',charpres2,'dc',dc,'Rbeq',Rbeq,'vol', vol,&
                                    'charvol',charvol,'in rank', proc_rank,&
                                    'cell m n p', cell(1), cell(2), cell(3)
            if (charvol2==0) print*, 'charvol2 is equal to zero, q_particle%sf is ', &
                q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3)), cellaux(1), cellaux(2), cellaux(3)
            if (charpres2==0) print*, 'charpres2 is equal to zero, pres%sf is ', pres%sf(cellaux(1),cellaux(2),cellaux(3))
            call s_mpi_abort()
        end if

    end function get_pinf1

    function deriv_gaspressure( bubbletmp, vaporflux, heatflux, gammabubble ) 

        type(particledata) :: bubbletmp
        real(kind(0.d0)) :: deriv_gaspressure, vaporflux,gammabubble, heatflux

        deriv_gaspressure = bubbletmp%tmp%p*bubbletmp%tmp%y(2) - heatflux - Runiv/MWvap*vaporflux
        deriv_gaspressure = -3.0d0*gammabubble/bubbletmp%tmp%y(1)*deriv_gaspressure 

    end function deriv_gaspressure

    subroutine s_deallocate_particles()

        integer :: i,j,k,imax
        type(particlenode),pointer :: particle
  
        if ((solverapproach.eq.2).and.avgdensflag) then
            imax = 4
            if(clusterflag.ge.4) imax=10 !subgrid noise model
        else if ((solverapproach.eq.0).and.avgdensflag) then
            imax = 3
        else
            imax = 2
        end if
        do i=1,imax
            deallocate(q_particle(i)%sf)
        end do
        deallocate(q_particle)
        particle  => particlesublist%list%next
        do while(associated(particle))
          call remove_particle (particle,2)
        end do
        do k=-buff_size,p+buff_size
            do j=-buff_size,n+buff_size
                do i=-buff_size,m+buff_size
                    deallocate ( qbl%fp(i,j,k)%list )
                end do
            end do
        end do
        deallocate( qbl%fp )
        deallocate(cellwblist%list)
        deallocate(cellwblist)
  
    end subroutine s_deallocate_particles  

end MODULE m_particles
