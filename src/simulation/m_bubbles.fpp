!>
!! @file m_bubbles.f90
!! @brief Contains module m_bubbles

#:include 'macros.fpp'

!> @brief This module contains the procedures shared by the ensemble-averaged and volume-averaged bubble models.
module m_bubbles

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    real(wp) :: chi_vw  !< Bubble wall properties (Ando 2010)
    real(wp) :: k_mw    !< Bubble wall properties (Ando 2010)
    real(wp) :: rho_mw  !< Bubble wall properties (Ando 2010)
    $:GPU_DECLARE(create='[chi_vw,k_mw,rho_mw]')

contains

    !> Function that computes the bubble radial acceleration based on bubble models
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param fCson Speed of sound from fP (EL)
    pure elemental function f_rddot(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, fntait, fBtait, f_bub_adv_src, f_divu, fCson, fInt, fshell, fRbuck, fRcell)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(in) :: fntait, fBtait, f_bub_adv_src, f_divu
        real(wp), intent(in) :: fshell, fRbuck, fCson, fInt, fRcell

        real(wp) :: fCpbw, fCpinf, fCpinf_dot, fH, fHdot, c_gas, c_liquid
        real(wp) :: pout
        real(wp) :: f_rddot

        if (bubble_model == 1) then
            ! Gilmore bubbles
            fCpinf = fP - pref
            fCpbw = f_cpbw(fR0, fR, fV, fpb)
            fH = f_H(fCpbw, fCpinf, fntait, fBtait)
            c_gas = f_cgas(fCpinf, fntait, fBtait, fH)
            fCpinf_dot = f_cpinfdot(fRho, fP, alf, fntait, fBtait, f_bub_adv_src, f_divu)
            fHdot = f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
            f_rddot = f_rddot_G(fCpbw, fR, fV, fH, fHdot, c_gas, fntait, fBtait)
        else if (bubble_model == 2) then
            ! Keller-Miksis bubbles
            fCpinf = fP
            fCpbw = f_cpbw_KM(fR0, fR, fV, fpb, fshell, fRbuck)
            if (bubbles_euler) then
                c_liquid = sqrt(fntait*(fP + fBtait)/(fRho*(1._wp - alf)))
            else
                if (lag_params%pressure_corrector .and. any(lag_params%interaction_model == (/1, 3/)) &
                    .and. adap_dt) then
                    pout = f_pout(fpbdot, fP, fCpbw, fRho, fR, fV, fshell, fRcell)
                    fCpinf = fCpinf + pout
                end if
                c_liquid = fCson
            end if
            f_rddot = f_rddot_KM(fpbdot, fCpinf, fCpbw, fRho, fR, fV, fR0, c_liquid, fInt, fshell, fRbuck, fpb)
        else if (bubble_model == 3) then
            ! Rayleigh-Plesset bubbles
            fCpbw = f_cpbw_KM(fR0, fR, fV, fpb)
            f_rddot = f_rddot_RP(fP, fRho, fR, fV, fCpbw)
        end if

    end function f_rddot

    !>  Function that computes that bubble wall pressure for Gilmore bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    pure elemental function f_cpbw(fR0, fR, fV, fpb)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fR0, fR, fV, fpb

        real(wp) :: f_cpbw

        if (polytropic) then
            f_cpbw = (Ca + 2._wp/Web/fR0)*((fR0/fR)**(3._wp*gam)) - Ca - 4._wp*Re_inv*fV/fR - 2._wp/(fR*Web)
        else
            f_cpbw = fpb - 1._wp - 4._wp*Re_inv*fV/fR - 2._wp/(fR*Web)
        end if

    end function f_cpbw

    !>  Function that computes the bubble enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
    pure elemental function f_H(fCpbw, fCpinf, fntait, fBtait)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpbw, fCpinf, fntait, fBtait

        real(wp) :: tmp1, tmp2, tmp3
        real(wp) :: f_H

        tmp1 = (fntait - 1._wp)/fntait
        tmp2 = (fCpbw/(1._wp + fBtait) + 1._wp)**tmp1
        tmp3 = (fCpinf/(1._wp + fBtait) + 1._wp)**tmp1

        f_H = (tmp2 - tmp3)*fntait*(1._wp + fBtait)/(fntait - 1._wp)

    end function f_H

    !> Function that computes the sound speed for the bubble
        !! @param fCpinf Driving bubble pressure
        !! @param fntait Tait EOS parameter
        !! @param fBtait Tait EOS parameter
        !! @param fH Bubble enthalpy
    pure elemental function f_cgas(fCpinf, fntait, fBtait, fH)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpinf, fntait, fBtait, fH

        real(wp) :: tmp
        real(wp) :: f_cgas

        ! get sound speed for Gilmore equations "C" -> c_gas
        tmp = (fCpinf/(1._wp + fBtait) + 1._wp)**((fntait - 1._wp)/fntait)
        tmp = fntait*(1._wp + fBtait)*tmp

        f_cgas = sqrt(tmp + (fntait - 1._wp)*fH)

    end function f_cgas

    !>  Function that computes the time derivative of the driving pressure
        !!  @param fRho Local liquid density
        !!  @param fP Local pressure
        !!  @param falf Local void fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param advsrc Advection equation source term
        !!  @param divu Divergence of velocity
    pure elemental function f_cpinfdot(fRho, fP, falf, fntait, fBtait, advsrc, divu)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fRho, fP, falf, fntait, fBtait, advsrc, divu

        real(wp) :: c2_liquid
        real(wp) :: f_cpinfdot

        ! get sound speed squared for liquid (only needed for pbdot)
        ! c_l^2 = gam (p+B) / (rho*(1-alf))
        if (mpp_lim) then
            c2_liquid = fntait*(fP + fBtait)/fRho
        else
            c2_liquid = fntait*(fP + fBtait)/(fRho*(1._wp - falf))
        end if

        ! \dot{Cp_inf} = rho sound^2 (alf_src - divu)
        f_cpinfdot = fRho*c2_liquid*(advsrc - divu)

    end function f_cpinfdot

    !>  Function that computes the time derivative of the enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fCpinf_dot Time derivative of the driving pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpbdot Time derivative of the internal bubble pressure
    pure elemental function f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpbw, fCpinf, fCpinf_dot, fntait, fBtait
        real(wp), intent(in) :: fR, fV, fR0, fpbdot

        real(wp) :: tmp1, tmp2
        real(wp) :: f_Hdot

        if (polytropic) then
            tmp1 = (fR0/fR)**(3._wp*gam)
            tmp1 = -3._wp*gam*(Ca + 2._wp/Web/fR0)*tmp1*fV/fR
        else
            tmp1 = fpbdot
        end if
        tmp2 = (2._wp/Web + 4._wp*Re_inv*fV)*fV/(fR**2._wp)

        f_Hdot = &
            (fCpbw/(1._wp + fBtait) + 1._wp)**(-1._wp/fntait)*(tmp1 + tmp2) &
            - (fCpinf/(1._wp + fBtait) + 1._wp)**(-1._wp/fntait)*fCpinf_dot

        ! Hdot = (Cpbw/(1+B) + 1)^(-1/n_tait)*(-3 gam)*(R0/R)^(3gam) V/R
        !f_Hdot = ((fCpbw/(1._wp+fBtait)+1._wp)**(-1._wp/fntait))*(-3._wp)*gam * &
        !            ( (fR0/fR)**(3._wp*gam ))*(fV/fR)

        ! Hdot = Hdot - (Cpinf/(1+B) + 1)^(-1/n_tait) Cpinfdot
        !f_Hdot = f_Hdot - ((fCpinf/(1._wp+fBtait)+1._wp)**(-1._wp/fntait))*fCpinf_dot

    end function f_Hdot

    !>  Function that computes the bubble radial acceleration for Rayleigh-Plesset bubbles
        !!  @param fCp Driving pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fCpbw Boundary wall pressure
    pure elemental function f_rddot_RP(fCp, fRho, fR, fV, fCpbw)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCp, fRho, fR, fV, fCpbw

        real(wp) :: f_rddot_RP

            !! rddot = (1/r) (  -3/2 rdot^2 + ((r0/r)^3\gamma - Cp)/rho )
            !! rddot = (1/r) (  -3/2 rdot^2 + (tmp1 - Cp)/rho )
            !! rddot = (1/r) (  tmp2 )

        f_rddot_RP = (-1.5_wp*(fV**2._wp) + (fCpbw - fCp)/fRho)/fR

    end function f_rddot_RP

    !>  Function that computes the bubble radial acceleration
        !!  @param fCpbw Bubble wall pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fH Current enthalpy
        !!  @param fHdot Current time derivative of the enthalpy
        !!  @param fcgas Current gas sound speed
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
    pure elemental function f_rddot_G(fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpbw, fR, fV, fH, fHdot
        real(wp), intent(in) :: fcgas, fntait, fBtait

        real(wp) :: tmp1, tmp2, tmp3
        real(wp) :: f_rddot_G

        tmp1 = fV/fcgas
        tmp2 = 1._wp + 4._wp*Re_inv/fcgas/fR*(fCpbw/(1._wp + fBtait) + 1._wp) &
               **(-1._wp/fntait)
        tmp3 = 1.5_wp*fV**2._wp*(tmp1/3._wp - 1._wp) + fH*(1._wp + tmp1) &
               + fR*fHdot*(1._wp - tmp1)/fcgas

        f_rddot_G = tmp3/(fR*(1._wp - tmp1)*tmp2)

    end function f_rddot_G

    !>  Function that computes the bubble wall pressure for Keller--Miksis bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure (EL polytropc: initial internal pressure)
    pure elemental function f_cpbw_KM(fR0, fR, fV, fpb, fshell, fRbuck)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fR0, fR, fV, fpb
        real(wp), intent(in), optional :: fshell, fRbuck
        real(wp) :: f_cpbw_KM
        real(wp) :: ss_mod

        if (polytropic) then
            if (bubbles_lagrange) then
                f_cpbw_KM = pv + (fpb-pv)*((fR0/fR)**(3._wp*gamma_m))
            else
                f_cpbw_KM = Ca*((fR0/fR)**(3._wp*gam)) - Ca + 1._wp
                if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM + &
                                                     (2._wp/(Web*fR0))*((fR0/fR)**(3._wp*gam))
            end if
        else
            f_cpbw_KM = fpb
        end if

        if (.not. f_is_default(Re_inv)) f_cpbw_KM = f_cpbw_KM - 4._wp*Re_inv*fV/fR

        ! Coated bubbles (buckling and elastic regimes)
        if (bubbles_lagrange .and. fshell == 1._wp) then
            ss_mod = 0._wp
            if (fR > fRbuck) ss_mod = lag_params%srfElast_ctdBub*((fR/fRbuck)**2._wp - 1._wp)
            f_cpbw_KM = f_cpbw_KM - 2._wp*ss_mod/fR - 4._wp*lag_params%srfDilVsc_ctdBub*fV/(fR**2._wp)
        else
            if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM - 2._wp/(fR*Web)
        end if

    end function f_cpbw_KM

    !>  Function that computes the bubble radial acceleration for Keller--Miksis bubbles
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param fCp Driving pressure
        !!  @param fCpbw Bubble wall pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fC Current sound speed
    pure elemental function f_rddot_KM(fpbdot, fCp, fCpbw, fRho, fR, fV, fR0, fC, fInt, fshell, fRbuck, fpb)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fpbdot, fCp, fCpbw
        real(wp), intent(in) :: fRho, fR, fV, fR0, fC, fInt, fshell, fRbuck, fpb

        real(wp) :: tmp1, tmp2, denom, cdot_star, ss_mod
        real(wp) :: f_rddot_KM
        if (polytropic) then
            if (bubbles_lagrange) then
                cdot_star = -(3._wp*gamma_m/fR)*fV*(fpb-pv)*((fR0/fR)**(3._wp*gamma_m))
            else
                cdot_star = -3._wp*gam*Ca*((fR0/fR)**(3._wp*gam))*fV/fR
                if (.not. f_is_default(Web)) cdot_star = cdot_star - &
                                                     3._wp*gam*(2._wp/(Web*fR0))*((fR0/fR)**(3._wp*gam))*fV/fR
            end if
        else
            cdot_star = fpbdot
        end if

        tmp1 = fV/fC
        denom = fR*(1._wp - tmp1)

        ! Coated bubbles
        if (bubbles_lagrange .and. fshell == 1._wp) then
            ss_mod = 0._wp
            if (fR > fRbuck) then ! Elastic regime
                ss_mod = lag_params%srfElast_ctdBub*((fR/fRbuck)**2._wp - 1._wp)
                cdot_star = cdot_star - 4._wp*lag_params%srfElast_ctdBub*fV/fRbuck**2._wp
            end if
            denom = denom + 4._wp*lag_params%srfDilVsc_ctdBub/(fRho*fC*fR)
            cdot_star = cdot_star + 2._wp*ss_mod*fV/(fR**2._wp) + &
                        8._wp*lag_params%srfDilVsc_ctdBub*(fV**2._wp)/fR**3._wp
        else
            if (.not. f_is_default(Web)) cdot_star = cdot_star + (2._wp/Web)*fV/(fR**2._wp)
        end if

        if (.not. f_is_default(Re_inv)) cdot_star = cdot_star + 4._wp*Re_inv*((fV/fR)**2._wp)

        tmp2 = 1.5_wp*(fV**2._wp)*(tmp1/3._wp - 1._wp) + &
               (1._wp + tmp1)*(fCpbw - fCp)/fRho + &
               cdot_star*fR/(fRho*fC)

        if (lag_params%pressure_corrector .and. lag_params%interaction_model == 2) tmp2 = tmp2 + fInt

        if (.not. f_is_default(Re_inv)) denom = denom + 4._wp*Re_inv/(fRho*fC)

        f_rddot_KM = tmp2/denom

    end function f_rddot_KM

    !>  Subroutine that computes bubble wall properties for vapor bubbles
        !!  @param pb Internal bubble pressure
        !!  @param iR0 Current bubble size index
    pure elemental subroutine s_bwproperty(pb_in, iR0, chi_vw_out, k_mw_out, rho_mw_out)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: pb_in
        integer, intent(in) :: iR0
        real(wp), intent(out) :: chi_vw_out  !< Bubble wall properties (Ando 2010)
        real(wp), intent(out) :: k_mw_out    !< Bubble wall properties (Ando 2010)
        real(wp), intent(out) :: rho_mw_out  !< Bubble wall properties (Ando 2010)
        real(wp) :: x_vw

        ! mass fraction of vapor
        chi_vw_out = 1._wp/(1._wp + R_v/R_n*(pb_in/pv - 1._wp))
        ! mole fraction of vapor & thermal conductivity of gas mixture
        x_vw = M_n*chi_vw_out/(M_v + (M_n - M_v)*chi_vw_out)
        k_mw_out = x_vw*k_v(iR0)/(x_vw + (1._wp - x_vw)*phi_vn) &
                   + (1._wp - x_vw)*k_n(iR0)/(x_vw*phi_nv + 1._wp - x_vw)
        ! gas mixture density
        rho_mw_out = pv/(chi_vw_out*R_v*Tw)

    end subroutine s_bwproperty

    !>  Function that computes the vapour flux
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index (EE) or bubble identifier (EL)
        !!  @param fmass_n Current gas mass (EL)
        !!  @param fbeta_c Mass transfer coefficient (EL)
        !!  @param fR_m Mixture gas constant (EL)
        !!  @param fgamma_m Mixture gamma (EL)
    pure elemental subroutine s_vflux(fR, fV, fpb, fmass_v, iR0, vflux, fmass_n, fbeta_c, fR_m, fgamma_m, fshell)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fR
        real(wp), intent(in) :: fV
        real(wp), intent(in) :: fpb, fmass_v
        integer, intent(in) :: iR0
        real(wp), intent(out) :: vflux
        real(wp), intent(out), optional :: fR_m, fgamma_m
        real(wp), intent(in), optional :: fmass_n, fbeta_c, fshell

        real(wp) :: chi_bar
        real(wp) :: rho_mw_lag
        real(wp) :: grad_chi
        real(wp) :: conc_v

        if (thermal == 3) then !transfer
            ! constant transfer model
            if (bubbles_lagrange) then
                ! Mixture properties (gas+vapor) in the bubble
                conc_v = fmass_v/(fmass_v + fmass_n)
                if (lag_params%massTransfer_model .and. (fshell == 0._wp)) then
                    conc_v = 1._wp/(1._wp + (R_v/R_n)*(fpb/pv - 1._wp))
                end if
                fR_m = (fmass_n*R_n + fmass_v*R_v)
                fgamma_m = conc_v*gamma_v + (1._wp - conc_v)*gamma_n

                ! Vapor flux
                chi_bar = fmass_v/(fmass_v + fmass_n)
                grad_chi = (chi_bar - conc_v)
                rho_mw_lag = (fmass_n + fmass_v)/(4._wp/3._wp*pi*fR**3._wp)
                vflux = 0._wp
                if (lag_params%massTransfer_model .and. (fshell == 0._wp)) then
                    vflux = -fbeta_c*rho_mw_lag*grad_chi/(1._wp - conc_v)/fR
                end if
            else
                chi_bar = fmass_v/(fmass_v + mass_n0(iR0))
                grad_chi = -Re_trans_c(iR0)*(chi_bar - chi_vw)
                vflux = rho_mw*grad_chi/Pe_c/(1._wp - chi_vw)/fR
            end if
        else
            ! polytropic
            vflux = pv*fV/(R_v*Tw)
        end if

    end subroutine s_vflux

    !>  Function that computes the time derivative of
        !!  the internal bubble pressure
        !!  @param fvflux Vapour flux
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Current internal bubble pressure
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index (EE) or bubble identifier (EL)
        !!  @param fbeta_t Heat transfer coefficient (EL)
        !!  @param fR_m Mixture gas constant (EL)
        !!  @param fgamma_m Mixture gamma (EL)
    pure elemental function f_bpres_dot(fvflux, fR, fV, fpb, fmass_v, iR0, fbeta_t, fR_m, fgamma_m, fshell)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fvflux
        real(wp), intent(in) :: fR
        real(wp), intent(in) :: fV
        real(wp), intent(in) :: fpb
        real(wp), intent(in) :: fmass_v
        integer, intent(in) :: iR0
        real(wp), intent(in), optional :: fbeta_t, fR_m, fgamma_m, fshell

        real(wp) :: T_bar
        real(wp) :: grad_T
        real(wp) :: f_bpres_dot
        real(wp) :: heatflux

        if (thermal == 3) then
            if (bubbles_lagrange) then
                T_bar = fpb*(4._wp/3._wp*pi*fR**3._wp)/fR_m
                grad_T = -fbeta_t*(T_bar - Tw)
                heatflux = 0._wp
                if (lag_params%heatTransfer_model .and. (fshell == 0._wp)) then
                    heatflux = (fgamma_m - 1._wp)/fgamma_m*grad_T/fR
                end if
                f_bpres_dot = 3._wp*fgamma_m*(-fV*fpb + fvflux*R_v*Tw &
                                              + heatflux)/fR
                return
            end if
            T_bar = Tw*(fpb/pb0(iR0))*(fR/R0(iR0))**3 &
                    *(mass_n0(iR0) + mass_v0(iR0))/(mass_n0(iR0) + fmass_v)
            grad_T = -Re_trans_T(iR0)*(T_bar - Tw)
            f_bpres_dot = 3._wp*gamma_m*(-fV*fpb + fvflux*R_v*Tw &
                                         + pb0(iR0)*k_mw*grad_T/Pe_T(iR0)/fR)/fR
        else
            f_bpres_dot = -3._wp*gamma_m*fV/fR*(fpb - pv)
        end if

    end function f_bpres_dot

        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param fCp Driving pressure
        !!  @param fCpbw Bubble wall pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fC Current sound speed
    pure elemental function f_pout(fpbdot, fCp, fCpbw, fRho, fR, fV, fshell, fRcell)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fpbdot, fCp, fCpbw, fRho, fR, fV, fRcell, fshell
        real(wp) :: f_pout, f_pout_2

        real(wp) :: c1, c2
        real(wp) :: c1_dot, c2_dot
        real(wp) :: aux, denom, ks, dphidt
        real(wp) :: a1, a2, a3

        f_pout = 0._wp

        if (p == 0) return

        aux = fRcell**3._wp - fR**3._wp
        c2 = 1.5_wp*(fR**3._wp)*(1._wp - fR/fRcell)/aux
        c1 = 1.5_wp*(fR*(fRcell**2._wp - fR**2._wp))/aux

        dphidt = fCpbw + 0.5_wp*fV**2._wp
        dphidt = (fCp - dphidt) + c2*fV**2._wp
        dphidt = dphidt/(1._wp - c1)

        f_pout = c1*dphidt + c2*fV**2._wp

        !!!!! Errors with MFC
        !Find Pout to modif Pinf Errors with MFC
        ! aux = fRcell**3._wp - fR**3._wp
        ! c2 = 1.5_wp*(fR**3._wp)*(1._wp - fR/fRcell)/aux
        ! c1 = 1.5_wp*(fR*(fRcell**2._wp - fR**2._wp))/aux

        ! dphidt = fCpbw/fRho - 0.5_wp*fV**2._wp
        ! dphidt = (fCp/fRho - dphidt) - c2*fV**2._wp
        ! dphidt = dphidt/(1._wp - c1)

        ! f_pout = fRho*(c1*dphidt - c2*fV**2._wp)    ! p_inf = pcell - pout

        !!!! Errors with MFC

        ! f_pout_2 = fCp - fCpbw - (c2 - 0.5_wp)*fRho*fV**2._wp
        ! f_pout_2 = f_pout_2/(1._wp - c1)

        ! f_pout = -f_pout_2 !+ 0.5_wp*fRho*fV**2._wp

        ! !Find Pinf_dot
        ! c1_dot = 2._wp*c1 - 3._wp + (fRcell/fR)**2._wp
        ! c1_dot = c1_dot*1.5_wp*(fV*fR**2._wp)/aux
        ! c2_dot = fV/fR + (fV*fR**2._wp)/aux
        ! c2_dot = c2_dot*3._wp*c2 - 1.5_wp*(fV*fR**3._wp)/(fRcell*aux)

        ! ks = lag_params%srfDilVsc_ctdBub

        ! a1 = fV*(1._wp - 2._wp*c2)
        ! if (.not. f_is_default(Re_inv)) a1 = a1 + 4._wp*Re_inv/(fRho*fR)
        ! if (fshell == 1._wp) a1 = a1 + 4._wp*ks/(fRho*fR**2._wp)
        ! a1 = a1/(1._wp - c1)

        ! a2 = -fpbdot
        ! if (.not. f_is_default(Re_inv)) a2 = a2 - 4._wp*Re_inv*(fV/fR)**2._wp
        ! if (.not. f_is_default(Web)) a2 = a2 - (2._wp/Web)*fV/fR**2._wp
        ! if (fshell == 1._wp) a2 = a2 - 4._wp*ks*(fV/(fR**2._wp))**2._wp
        ! a2 = a2/fRho
        ! a2 = a2 - c2_dot*fV**2._wp + dphidt*c1_dot
        ! a2 = a2/(1._wp - c1)

        ! a3 = dphidt*c1_dot - 2._wp*c2_dot*fV**2._wp

        ! fPinfdot_star = (c1*a2 + a3)*fRho
        ! fPinfdot_denom = (-c1*a1 + 2._wp*c2*fV)*fRho

    end function f_pout

    ! function f_pres_stochastic(fTzPcell, fnoise_constant, flambda_c, fdk, floc, ftime, fCson)!, fPhase_rn)
    !     !$acc routine seq
    !     real(wp), intent(in) :: fTzPcell, fnoise_constant, flambda_c, fdk, floc, ftime, fCson
    !     !real(wp), dimension(num_noise), intent(in) :: fPhase_rn

    !     real(wp) :: f_pres_stochastic
    !     real(wp) :: constant_term, k, angFreq, rndPhase, A_k_sqrd
    !     integer :: i

    !     f_pres_stochastic = 0._wp

    !     constant_term = (fnoise_constant/(0.5_wp*flambda_c*sqrt(2_wp*pi)))

    !     if (constant_term <= 0._wp) return ! Avoid complex numbers when taking squared root of negative A_k_sqrd

    !     k = 0._wp
    !     do i = 1, num_noise
    !         rndPhase = f_random_normal(0.5_wp*pi, lag_params%pnoise_dev, 0._wp, 2._wp*pi) ! mean, dev, min, max
    !         A_k_sqrd = constant_term * exp(-0.5_wp*((2._wp*pi/k - flambda_c)/(0.5_wp*flambda_c))**2._wp)
    !         f_pres_stochastic = f_pres_stochastic + sqrt(A_k_sqrd) * fdk * cos(k*floc - k*fCson*ftime + rndPhase)!+ fPhase_rn(i))
    !         if (f_pres_stochastic /= f_pres_stochastic) then
    !             print*, i, k, A_k_sqrd, sqrt(A_k_sqrd), f_pres_stochastic
    !             stop "f_pres_stochastic is NaN"
    !         end if
    !         k = k + fdk
    !     end do

    ! end function f_pres_stochastic

    ! function f_random_normal(fmean, fdev, fmin, fmax)
    !     !$acc routine seq
    !     real(wp), intent(in) :: fmean, fdev, fmin, fmax

    !     real(wp) :: f_random_normal
    !     real(wp) :: num_rn1, num_rn2

    !     do while (.true.)

    !         call random_number(num_rn1)
    !         num_rn1 = 1._wp - num_rn1
    !         call random_number(num_rn2)
    !         num_rn2 = 1._wp - num_rn2

    !         f_random_normal = fdev*sqrt(-2._wp*log(num_rn1))*cos(2._wp*pi*num_rn2) + fmean

    !         if (f_random_normal >= fmin .and. f_random_normal <= fmax) exit

    !     end do

    ! end function f_random_normal

    !> Adaptive time stepping routine for subgrid bubbles
        !!  (See Heirer, E. Hairer S.P.Nørsett G. Wanner, Solving Ordinary
        !!  Differential Equations I, Chapter II.4)
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param bub_id Bubble identifier (EL)
        !!  @param fmass_v Current mass of vapour (EL)
        !!  @param fmass_n Current mass of gas (EL)
        !!  @param fbeta_c Mass transfer coefficient (EL)
        !!  @param fbeta_t Heat transfer coefficient (EL)
        !!  @param fCson Speed of sound (EL)
        !!  @param adap_dt_stop Fail-safe exit if max iteration count reached
    impure subroutine s_advance_step(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                   fntait, fBtait, f_bub_adv_src, f_divu, &
                                   bub_id, fmass_v, fmass_n, fbeta_c, &
                                   fbeta_t, fCson, fInt, fshell, fRbuck, fRrupt, fRcell, &
                                   fnoise_constant, flambda_c, fdk, floc, ftime, fAc, &!fPhase_rn, &
                                   fQvis, fQth, fRmean, adap_dt_stop)
        $:GPU_ROUTINE(function_name='s_advance_step',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(inout) :: fR, fV, fpb, fmass_v, fshell
        real(wp), intent(in) :: fRho, fP, fR0, fpbdot, alf
        real(wp), intent(in) :: fntait, fBtait, f_bub_adv_src, f_divu
        integer, intent(in) :: bub_id
        real(wp), intent(out) :: fAc
        real(wp), intent(in) :: fmass_n, fbeta_c, fbeta_t, fCson, fInt, fRbuck, fRrupt, fRcell
        real(wp), intent(in) :: fnoise_constant, flambda_c, fdk, floc, ftime
        !real(wp), dimension(num_noise), intent(in) :: fPhase_rn
        real(wp), intent(out) :: fQvis, fQth, fRmean
        integer, intent(inout) :: adap_dt_stop

        real(wp), dimension(5) :: err !< Error estimates for adaptive time stepping
        real(wp) :: t_new !< Updated time step size
        real(wp) :: h !< Time step size
        real(wp), dimension(4) :: myR_tmp1, myV_tmp1, myR_tmp2, myV_tmp2 !< Bubble radius, radial velocity, and radial acceleration for the inner loop
        real(wp), dimension(4) :: myPb_tmp1, myMv_tmp1, myPb_tmp2, myMv_tmp2 !< Gas pressure and vapor mass for the inner loop (EL)

        real(wp) :: fR2, fV2, fpb2, fmass_v2
        integer :: iter_count
        real(wp) :: conc_v_h, R_m_h, gamma_m_h, T_bar_h, grad_T_h, heatflux_h
        real(wp) :: fAc1, fAc21, fAc22, fvis_inst, fth_inst

        call s_initial_substep_h(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                 fntait, fBtait, f_bub_adv_src, f_divu, fCson, fInt, fshell, fRbuck, fRcell, h)

        if (h /= h) then
            ! print *, 'Altering initial h (from/to)', bub_id, h, 0.1_wp*0.5_wp*dt
            ! print *, fRho, fP, fR, fV, fR0, fpb, fpbdot, fCson, fInt, fshell, fRbuck, fRcell
            h = 0.1_wp*0.5_wp*dt
        end if

        ! Advancing one step
        t_new = 0._wp
        fQvis = 0._wp
        fQth = 0._wp
        fRmean = 0._wp
        fAc = 0._wp
        iter_count = 0
        adap_dt_stop = 0

        !print*, fRho, fP, fR, fV, fR0, fpb, gamma_m, gam

        do
            if (t_new + h > 0.5_wp*dt) then
                h = 0.5_wp*dt - t_new
            end if

            ! Advancing one sub-step
            do while (iter_count < adap_dt_max_iters)

                iter_count = iter_count + 1

                ! Advance one sub-step
                call s_advance_substep(err(1), &
                                       fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                       fntait, fBtait, f_bub_adv_src, f_divu, &
                                       bub_id, fmass_v, fmass_n, fbeta_c, fbeta_t, &
                                       fnoise_constant, flambda_c, fdk, floc, ftime + t_new, fAc1, & !fPhase_rn&
                                       fCson, fInt, fshell, fRbuck, fRcell, h, &
                                       myR_tmp1, myV_tmp1, myPb_tmp1, myMv_tmp1)

                ! Advance one sub-step by advancing two half steps
                call s_advance_substep(err(2), &
                                       fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                       fntait, fBtait, f_bub_adv_src, f_divu, &
                                       bub_id, fmass_v, fmass_n, fbeta_c, fbeta_t, &
                                       fnoise_constant, flambda_c, fdk, floc, ftime + t_new, fAc21, & !fPhase_rn&
                                       fCson, fInt, fshell, fRbuck, fRcell, 0.5_wp*h, &
                                       myR_tmp2, myV_tmp2, myPb_tmp2, myMv_tmp2)

                fR2 = myR_tmp2(4); fV2 = myV_tmp2(4)
                fpb2 = myPb_tmp2(4); fmass_v2 = myMv_tmp2(4)

                call s_advance_substep(err(3), &
                                       fRho, fP, fR2, fV2, fR0, fpb2, fpbdot, alf, &
                                       fntait, fBtait, f_bub_adv_src, f_divu, &
                                       bub_id, fmass_v2, fmass_n, fbeta_c, fbeta_t, &
                                       fnoise_constant, flambda_c, fdk, floc, ftime + t_new + 0.5_wp*h, fAc22, & !fPhase_rn&
                                       fCson, fInt, fshell, fRbuck, fRcell, 0.5_wp*h, &
                                       myR_tmp2, myV_tmp2, myPb_tmp2, myMv_tmp2)

                err(4) = abs((myR_tmp1(4) - myR_tmp2(4))/myR_tmp1(4))
                err(5) = abs((myV_tmp1(4) - myV_tmp2(4))/myV_tmp1(4))
                if (abs(myV_tmp1(4)) < verysmall) err(5) = 0._wp

                ! Determine acceptance/rejection and update step size
                !   Rule 1: err1, err2, err3 < tol
                !   Rule 2: myR_tmp1(4) > 0._wp
                !   Rule 3: abs((myR_tmp1(4) - myR_tmp2(4))/fR) < tol
                !   Rule 4: abs((myV_tmp1(4) - myV_tmp2(4))/fV) < tol
                if ((err(1) <= adap_dt_tol) .and. (err(2) <= adap_dt_tol) .and. &
                    (err(3) <= adap_dt_tol) .and. (err(4) < adap_dt_tol) .and. &
                    (err(5) < adap_dt_tol) .and. myR_tmp1(4) > 0._wp) then

                    ! Accepted. Finalize the sub-step
                    t_new = t_new + h

                    ! Update R and V
                    fR = myR_tmp1(4)
                    fV = myV_tmp1(4)
                    fAc = 0.5_wp*fAc21 + 0.5_wp*fAc22

                    !if (bub_id == 50) then
                    !    print*, 'ss:', fR, fV, fAc, fAc21, fAc22, t_new
                    !end if

                    if (bubbles_lagrange) then
                        ! Update pb and mass_v
                        fpb = myPb_tmp1(4)
                        if (polytropic) then 
                          fpb = pv + (fpb - pv)*(fR0/fR)**(3._wp*gamma_m)
                        end if
                        fmass_v = myMv_tmp1(4)
                        if (fR > fRrupt) fshell = 0._wp

                        if (hifu_params%sampling) then

                            !> Mixture properties in the bubble
                            conc_v_h = 0._wp
                            if (lag_params%massTransfer_model .and. (fshell == 0._wp)) then
                                conc_v_h = 1._wp/(1._wp + (R_v/R_n)*(fpb/pv - 1._wp))
                            end if
                            R_m_h = fmass_n*R_n + fmass_v*R_v
                            gamma_m_h = conc_v_h*gamma_v + (1._wp - conc_v_h)*gamma_n

                            !> Viscous damping of the bubble (Watts)
                            fvis_inst = (4._wp*pi*fR**2._wp)*(4._wp*mul0*(fV**2._wp)/(fR))
                            fQvis = fQvis + h*fvis_inst

                            !> Thermal damping of the bubble (Watts)
                            if (.not. polytropic) then
                                T_bar_h = fpb*(4._wp/3._wp*pi*fR**3._wp)/R_m_h
                                grad_T_h = -fbeta_t*(T_bar_h - Tw)
                                if (lag_params%heatTransfer_model .and. (fshell == 0._wp)) then
                                    heatflux_h = (gamma_m_h - 1._wp)/gamma_m_h*grad_T_h/fR
                                end if
                            else
                                T_bar_h = Tw * (fR0/fR)**(3._wp*(gamma_m-1._wp)) ! Polytropic temp
                                heatflux_h = conc_v_h*k_vl + (1._wp - conc_v_h)*k_nl 
                                heatflux_h = 3._wp*heatflux_h*(1._wp-gamma_m)*T_bar_h/fR
                            end if
                            fth_inst = heatflux_h*4._wp*pi*fR**2._wp
                            fQth = fQth + h*fth_inst

                            !> Mean radius
                            fRmean = fRmean + h*fR

                            ! Checking for NaNs and negative qvis
                            if (fQvis /= fQvis .or. fQth /= fQth .or. fQvis < 0._wp) then
                                iter_count = adap_dt_max_iters
                            end if
                        end if
                    end if

                    ! Update step size for the next sub-step
                    h = h*min(2._wp, max(0.5_wp, (1e-4_wp/err(1))**(1._wp/3._wp)))

                    exit
                else
                    ! Rejected. Update step size for the next try on sub-step
                    if (err(2) <= adap_dt_tol) then
                        h = 0.5_wp*h
                    else
                        h = 0.25_wp*h
                    end if

                    ! if (iter_count >= adap_dt_max_iters) then
                    !     print *, 'h small', h, t_new, bub_id, fR, fV, fshell
                    !     print *, 'errs', err(1), err(2), err(3), err(4), err(5)
                    !     print *, 'tmp R', myR_tmp1(1), myR_tmp1(2), myR_tmp1(3), myR_tmp1(4)
                    !     print *, 'tmp V', myV_tmp1(1), myV_tmp1(2), myV_tmp1(3), myV_tmp1(4)
                    !     print *, 'tmp Pb', myPb_tmp1(1), myPb_tmp1(2), myPb_tmp1(3), myPb_tmp1(4)
                    !     print *, 'tmp mass_v', myMv_tmp1(1), myMv_tmp1(2), myMv_tmp1(3), myMv_tmp1(4)
                    !     print *, 'old vals', fR, fV, fpb, fmass_v, fP, t_new, t_new/(0.5_wp*dt)
                    !     print *, 'otherVars', fmass_n, fbeta_c, fbeta_t, fCson, fshell, fRbuck, fRrupt
                    ! end if

                end if
            end do

            ! Exit the loop if the final time reached dt
            if (f_approx_equal(t_new, 0.5_wp*dt) .or. iter_count >= adap_dt_max_iters) exit

        end do

        if (iter_count >= adap_dt_max_iters) adap_dt_stop = 1

        if (adap_dt_stop == 1) print*, iter_count, fR, fV, fR0, fP, fpb, gamma_m, err(1), err(2), err(3), err(4), err(5), h

    end subroutine s_advance_step

    !> Choose the initial time step size for the adaptive time stepping routine
        !!  (See Heirer, E. Hairer S.P.Nørsett G. Wanner, Solving Ordinary
        !!  Differential Equations I, Chapter II.4)
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param fCson Speed of sound (EL)
        !!  @param h Time step size
    pure subroutine s_initial_substep_h(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                        fntait, fBtait, f_bub_adv_src, f_divu, &
                                        fCson, fInt, fshell, fRbuck, fRcell, h)
        $:GPU_ROUTINE(function_name='s_initial_substep_h',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu
        real(wp), intent(IN) :: fCson, fshell, fRbuck, fInt, fRcell
        real(wp), intent(OUT) :: h

        real(wp), dimension(2) :: h_size !< Time step size (h0, h1)
        real(wp), dimension(3) :: d_norms !< norms (d_0, d_1, d_2)
        real(wp), dimension(2) :: myR_tmp, myV_tmp, myA_tmp !< Bubble radius, radial velocity, and radial acceleration

        ! Determine the starting time step
        ! Evaluate f(x0,y0)
        myR_tmp(1) = fR
        myV_tmp(1) = fV
        myA_tmp(1) = f_rddot(fRho, fP, myR_tmp(1), myV_tmp(1), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson, fInt, fshell, fRbuck, fRcell)

        ! Compute d_0 = ||y0|| and d_1 = ||f(x0,y0)||
        d_norms(1) = sqrt((myR_tmp(1)**2._wp + myV_tmp(1)**2._wp)/2._wp)
        d_norms(2) = sqrt((myV_tmp(1)**2._wp + myA_tmp(1)**2._wp)/2._wp)
        if (d_norms(1) < threshold_first_guess .or. d_norms(2) < threshold_first_guess) then
            h_size(1) = small_guess
        else
            h_size(1) = scale_guess*(d_norms(1)/d_norms(2))
        end if

        ! Evaluate f(x0+h0,y0+h0*f(x0,y0))
        myR_tmp(2) = myR_tmp(1) + h_size(1)*myV_tmp(1)
        myV_tmp(2) = myV_tmp(1) + h_size(1)*myA_tmp(1)
        myA_tmp(2) = f_rddot(fRho, fP, myR_tmp(2), myV_tmp(2), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson, fInt, fshell, fRbuck, fRcell)

        ! Compute d_2 = ||f(x0+h0,y0+h0*f(x0,y0))-f(x0,y0)||/h0
        d_norms(3) = sqrt(((myV_tmp(2) - myV_tmp(1))**2._wp + (myA_tmp(2) - myA_tmp(1))**2._wp)/2._wp)/h_size(1)

        ! Set h1 = (0.01/max(d_1,d_2))^{1/(p+1)}
        !      if max(d_1,d_2) < 1.e-15_wp, h_size(2) = max(1.e-6_wp, h0*1.e-3_wp)
        if (max(d_norms(2), d_norms(3)) < threshold_second_guess) then
            h_size(2) = max(small_guess, h_size(1)*scale_first_guess)
        else
            h_size(2) = (scale_guess/max(d_norms(2), d_norms(3)))**(1._wp/3._wp)
        end if

        h = min(h_size(1)/scale_guess, h_size(2))

    end subroutine s_initial_substep_h

    !>  Integrate bubble variables over the given time step size, h, using a
        !!      third-order accurate embedded Runge–Kutta scheme.
        !!  @param err Estimated error
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param bub_id Bubble identifier (EL)
        !!  @param fmass_v Current mass of vapour (EL)
        !!  @param fmass_n Current mass of gas (EL)
        !!  @param fbeta_c Mass transfer coefficient (EL)
        !!  @param fbeta_t Heat transfer coefficient (EL)
        !!  @param fCson Speed of sound (EL)
        !!  @param fshell Shell switch, Marmottant model (EL)
        !!  @param fRbuck Buckling radius, Marmottant model (EL)
        !!  @param h Time step size
        !!  @param myR_tmp Bubble radius at each stage
        !!  @param myV_tmp Bubble radial velocity at each stage
        !!  @param myPb_tmp Internal bubble pressure at each stage (EL)
        !!  @param myMv_tmp Mass of vapor in the bubble at each stage (EL)
    pure subroutine s_advance_substep(err, fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                      fntait, fBtait, f_bub_adv_src, f_divu, &
                                      bub_id, fmass_v, fmass_n, fbeta_c, fbeta_t, &
                                      fnoise_constant, flambda_c, fdk, floc, ftime, fAc, & !fPhase_rn &
                                      fCson, fInt, fshell, fRbuck, fRcell, h, &
                                      myR_tmp, myV_tmp, myPb_tmp, myMv_tmp)
        $:GPU_ROUTINE(function_name='s_advance_substep',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(OUT) :: err
        real(wp), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu, h
        integer, intent(IN) :: bub_id
        real(wp), intent(out) :: fAc
        real(wp), intent(IN) :: fmass_v, fmass_n, fbeta_c, fbeta_t, fCson, fshell, fRbuck, fInt, fRcell
        real(wp), intent(in) :: fnoise_constant, flambda_c, fdk, floc, ftime
        real(wp), dimension(4), intent(OUT) :: myR_tmp, myV_tmp, myPb_tmp, myMv_tmp

        real(wp), dimension(4) :: myA_tmp, mydPbdt_tmp, mydMvdt_tmp
        real(wp) :: err_R, err_V, Pinf

        Pinf = fP
        ! if (bubbles_lagrange .and. num_dims == 2) then
        !     Pnoise_tmp = f_pres_stochastic(fP, fnoise_constant, flambda_c, fdk, floc, ftime, fCson)!, fPhase_rn)
        !     Pinf = fP + Pnoise_tmp
        ! end if

        myPb_tmp(1:4) = fpb
        mydPbdt_tmp(1:4) = fpbdot

        ! Stage 0
        myR_tmp(1) = fR
        myV_tmp(1) = fV
        if (bubbles_lagrange) then
            
            myPb_tmp(1) = fpb
            myMv_tmp(1) = fmass_v
            call s_advance_EL(myR_tmp(1), myV_tmp(1), myPb_tmp(1), myMv_tmp(1), bub_id, &
                              fmass_n, fbeta_c, fbeta_t, mydPbdt_tmp(1), mydMvdt_tmp(1), fshell)
        end if
        myA_tmp(1) = f_rddot(fRho, Pinf, myR_tmp(1), myV_tmp(1), fR0, &
                             myPb_tmp(1), mydPbdt_tmp(1), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson, fInt, fshell, fRbuck, fRcell)

        ! Stage 1
        myR_tmp(2) = myR_tmp(1) + h*myV_tmp(1)
        myV_tmp(2) = myV_tmp(1) + h*myA_tmp(1)
        if (bubbles_lagrange) then
            myPb_tmp(2) = myPb_tmp(1) + h*mydPbdt_tmp(1)
            myMv_tmp(2) = myMv_tmp(1) + h*mydMvdt_tmp(1)
            call s_advance_EL(myR_tmp(2), myV_tmp(2), myPb_tmp(2), myMv_tmp(2), &
                              bub_id, fmass_n, fbeta_c, fbeta_t, mydPbdt_tmp(2), mydMvdt_tmp(2), fshell)
        end if
        myA_tmp(2) = f_rddot(fRho, Pinf, myR_tmp(2), myV_tmp(2), fR0, &
                             myPb_tmp(2), mydPbdt_tmp(2), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson, fInt, fshell, fRbuck, fRcell)

        ! Stage 2
        myR_tmp(3) = myR_tmp(1) + (h/4._wp)*(myV_tmp(1) + myV_tmp(2))
        myV_tmp(3) = myV_tmp(1) + (h/4._wp)*(myA_tmp(1) + myA_tmp(2))
        if (bubbles_lagrange) then
            myPb_tmp(3) = myPb_tmp(1) + (h/4._wp)*(mydPbdt_tmp(1) + mydPbdt_tmp(2))
            myMv_tmp(3) = myMv_tmp(1) + (h/4._wp)*(mydMvdt_tmp(1) + mydMvdt_tmp(2))
            call s_advance_EL(myR_tmp(3), myV_tmp(3), myPb_tmp(3), myMv_tmp(3), &
                              bub_id, fmass_n, fbeta_c, fbeta_t, mydPbdt_tmp(3), mydMvdt_tmp(3), fshell)
        end if
        myA_tmp(3) = f_rddot(fRho, Pinf, myR_tmp(3), myV_tmp(3), fR0, &
                             myPb_tmp(3), mydPbdt_tmp(3), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson, fInt, fshell, fRbuck, fRcell)

        ! Stage 3
        myR_tmp(4) = myR_tmp(1) + (h/6._wp)*(myV_tmp(1) + myV_tmp(2) + 4._wp*myV_tmp(3))
        myV_tmp(4) = myV_tmp(1) + (h/6._wp)*(myA_tmp(1) + myA_tmp(2) + 4._wp*myA_tmp(3))
        if (bubbles_lagrange) then
            myPb_tmp(4) = myPb_tmp(1) + (h/6._wp)*(mydPbdt_tmp(1) + mydPbdt_tmp(2) + 4._wp*mydPbdt_tmp(3))
            myMv_tmp(4) = myMv_tmp(1) + (h/6._wp)*(mydMvdt_tmp(1) + mydMvdt_tmp(2) + 4._wp*mydMvdt_tmp(3))
            call s_advance_EL(myR_tmp(4), myV_tmp(4), myPb_tmp(4), myMv_tmp(4), &
                              bub_id, fmass_n, fbeta_c, fbeta_t, mydPbdt_tmp(4), mydMvdt_tmp(4), fshell)
        end if
        myA_tmp(4) = f_rddot(fRho, Pinf, myR_tmp(4), myV_tmp(4), fR0, &
                             myPb_tmp(4), mydPbdt_tmp(4), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson, fInt, fshell, fRbuck, fRcell)
        fAc = myA_tmp(4)
        ! Estimate error
        err_R = (-5._wp*h/24._wp)*(myV_tmp(2) + myV_tmp(3) - 2._wp*myV_tmp(4)) &
                /max(abs(myR_tmp(1)), abs(myR_tmp(4)))
        err_V = (-5._wp*h/24._wp)*(myA_tmp(2) + myA_tmp(3) - 2._wp*myA_tmp(4)) &
                /max(abs(myV_tmp(1)), abs(myV_tmp(4)))
        ! Error correction for non-oscillating bubbles
        if (bubbles_lagrange .and. f_approx_equal(myA_tmp(1), 0._wp) .and. f_approx_equal(myA_tmp(2), 0._wp) .and. &
            f_approx_equal(myA_tmp(3), 0._wp) .and. f_approx_equal(myA_tmp(4), 0._wp)) then
            err_V = 0._wp
        end if
        err = sqrt((err_R**2._wp + err_V**2._wp)/2._wp)

    end subroutine s_advance_substep

    !>  Changes of pressure and vapor mass in the lagrange bubbles.
        !!  @param bub_id Bubble identifier
        !!  @param fmass_n Current mass of gas
        !!  @param fbeta_c Mass transfer coefficient
        !!  @param fbeta_t Heat transfer coefficient
        !!  @param fR_tmp Bubble radius
        !!  @param fV_tmp Bubble radial velocity
        !!  @param fPb_tmp Internal bubble pressure
        !!  @param fMv_tmp Mass of vapor in the bubble
        !!  @param fdPbdt_tmp Rate of change of the internal bubble pressure
        !!  @param fdMvdt_tmp Rate of change of the mass of vapor in the bubble
    pure elemental subroutine s_advance_EL(fR_tmp, fV_tmp, fPb_tmp, fMv_tmp, bub_id, &
                                           fmass_n, fbeta_c, fbeta_t, fdPbdt_tmp, advance_EL, fshell)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(IN) :: fR_tmp, fV_tmp, fPb_tmp, fMv_tmp
        real(wp), intent(IN) :: fmass_n, fbeta_c, fbeta_t, fshell
        integer, intent(IN) :: bub_id
        real(wp), intent(INOUT) :: fdPbdt_tmp
        real(wp), intent(out) :: advance_EL
        real(wp) :: fVapFlux, myR_m, mygamma_m

        call s_vflux(fR_tmp, fV_tmp, fPb_tmp, fMv_tmp, bub_id, fVapFlux, fmass_n, fbeta_c, myR_m, mygamma_m, fshell)
        fdPbdt_tmp = f_bpres_dot(fVapFlux, fR_tmp, fV_tmp, fPb_tmp, fMv_tmp, bub_id, fbeta_t, myR_m, mygamma_m, fshell)
        advance_EL = 4._wp*pi*fR_tmp**2._wp*fVapFlux

        if (polytropic) then
            fdPbdt_tmp = 0._wp; advance_EL = 0._wp
        end if

    end subroutine s_advance_EL

end module m_bubbles
