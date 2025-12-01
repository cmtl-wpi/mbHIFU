#!/usr/bin/env python3

import math
import json

#################################################################
##### Evaporated milk phantom (with or w/o microbubbles)    #####
##### Tries to reproduce Kaleb's experiment (2024-10-10)    #####
##### Author: Diego Vaca Revelo                             #####
#################################################################

# Function to calculate the EOS parameters
# specific heat ratio and stiffness - Pa
def calculat_eos_param(rhoFun, cFun, cpFun, pFun, tFun):
    gammaFun = cFun**2/(cpFun*tFun) + 1
    pi_infFun = rhoFun * cpFun * tFun * (gammaFun-1)/gammaFun - pFun
    return gammaFun, pi_infFun

# Function to find dynamic viscosity from abs coefficient 
# (assume bulk viscosity to be thrice the dynamic viscosity)
def absCoef_to_mu(absFun, freqFun, rhoFun, cFun):
    omega = 2 * math.pi * freqFun
    A = (absFun * rhoFun * cFun**3)/(omega**2)
    muFun=(3/11)*(A)
    return muFun

# Define reference values for nondimensionalization (ok)
x0 = 1.e-03         # length - m
rho0 = 1.e+03       # density - kg/m3
c0 = 1475.          # speed of sound - m/s
p0 = rho0*c0*c0     # pressure - Pa
T0 = 298            # temperature - K

# Define acoustic source properties (ok)
# (Gain = 23.2817 and Pfocus = 1.99 MPa 'assume pure water')
patm = 101325.          # Atmospheric pressure - Pa
pamplitude = 2.0e6   # Amplitud of the acoustic source - Pa
freq = 1.11e+06         # Source frequency - Hz
focLen  = 46.0e-03      # Focal length - m
aperture= 41.5e-03      # Transducer aperture - m
waveLen = c0/freq       # wave length - m

# Define water properties (ok)
c_water = 1475.                 # speed of sound - m/s
rho_water = 1000                # density - kg/m3
T_water = 298                   # temperature - K
cp_water = 4180                 # specific heat - J/kg*K
tdiff_water = 1.46e-7           # thermal diffusivity - m2/s
abs_coef_water = 0.025          # attenuation of water -  Np/m at 1MHz
abs_coef_water = abs_coef_water * ((freq / 1.0e6)**2)                    # power 1 linear, 2 quadratic
mu_water = absCoef_to_mu(abs_coef_water, freq, rho_water, c_water)    # Dynamic viscosity - Pa.s
[gamma_water, pi_inf_water]  = calculat_eos_param(rho_water, c_water, cp_water, 1e5, T_water)   # specific heat ratio and stiffness - Pa

# Define host properties (EMP) (ok)
c_host = 1570.                  # speed of sound - m/s
rho_host = 1040                 # density kg/m3
T_host = 298                    # temperature K
cp_host = 3850                  # specific heat - J/kg*K
tdiff_host = 1.32e-7            # thermal diffusivity - m2/s
abs_coef_host = 4.20            # attenuation of water -  Np/m at 1MHz
abs_coef_host = abs_coef_host * ((freq / 1.0e6)**1)                  # power 1 linear, 2 quadratic
mu_host = absCoef_to_mu(abs_coef_host, freq, rho_host, c_host)    # Dynamic viscosity - Pa.s
[gamma_host, pi_inf_host]  = calculat_eos_param(rho_host, c_host, cp_host, 1e5, T_host)   # specific heat ratio and stiffness - Pa

# Lagrangian bubble's properties (ok)
    # Sonazoid contrast agent: 
    # gas in the core: Perfluorobutane (C4F10)
    # shell: monomolecular membrane of hydrogenated egg phosphatidylserine
    # vapor addition through the interface (if mass transfer is TRUE) is considered to be water vapor only (NIST properties)
R_uni = 8314                        # Universal gas constant - J/kmol/K
MW_g = 238.027                      # Molar weigth of the gas - kg/kmol (https://pubs.acs.org/doi/10.1021/acs.iecr.1c02969)
MW_v = 18.0                         # Molar weigth of the vapor - kg/kmol
gamma_g = 1.0699                    # Specific heat ratio of the gas (https://doi.org/10.3390/pharmaceutics14010098)
gamma_v = 1.333                     # Specific heat ratio of the vapor
pv = 2350                           # Vapor pressure of the host - Pa
cp_g = 0.809e+03                    # Specific heat of the gas - J/kg/K (https://www.f2chemicals.com/perfluorobutane.html)
cp_v = 2.1e+03                      # Specific heat of the vapor - J/kg/K
k_g = 0.00674                       # Thermal conductivity of the gas - W/m/K (https://doi.org/10.2172/6986083)
k_v = 0.02                          # Thermal conductivity of the vapor - W/m/K
diffVapor = 2.178e-05               # Diffusivity coefficient of the vapor in air - m2/s (https://doi.org/10.1016/S1352-2310(97)00391-9)
sigBubble = 0.069                   # Clean gas-water surface tension of the bubble (assume water) - N/m
# Marmotant model inputs:
sigmaInit = 0.0                     # Initial surface tension of the coated bubble - N/m (Assume initially R=Rbuck)
elasticity = 0.53                   # Elasticity of the shell - N/m (https://doi.org/10.3390/pharmaceutics14010098)
dilatationalViscosity = 1.2e-8      # Dilatiational viscosity of the shell - Pa.s.m (https://doi.org/10.1121/1.3418685)
mu_g = 1.48e-5

# Domain stgs 1 - 2
xb = 0.0e-03       # Domain boundaries - m (x direction)
xe = 6.0e-03
yb = 0.e-3         # Domain boundaries - m (y direction)
ye = 3.0e-03
zb = 0.e-3         # Domain boundaries - m (z direction)
ze = 3.0e-03

mltplo  = 10
waveSol = 15
Nx      = int(round(((xe-xb)/(waveLen/waveSol))/mltplo)*mltplo)      # number of elements into x direction (wavelength/35)
Ny      = int(round(((ye-yb)/(waveLen/waveSol))/mltplo)*mltplo)
Nz      = int(round(((ze-zb)/(waveLen/waveSol))/mltplo)*mltplo)

# Domain stg 3 - Heat solver cartesian coords
#xb_ht = 21.e-3 # Domain boundaries - m (x direction)
#xe_ht = 60.e-3
#ye_ht = 6.e-3          # Domain boundaries - m (y direction)
#ze_ht = 6.e-3          # Domain boundaries - m (z direction)

#Nx_ht = int(round(((xe_ht-xb_ht)/(waveLen/35))/50)*50)      # number of elements into x, y, z, dir
#Ny_ht = int(round(((2*ye_ht)/(waveLen/35))/50)*50)
#Nz_ht = int(round(((2*ze_ht)/(waveLen/35))/50)*50)

#scaleMesh = 1
#Nx_ht = Nx_ht / scaleMesh
#Ny_ht = Ny_ht / scaleMesh
#Nz_ht = Nz_ht / scaleMesh

########### Time set up #########

# E-L solver
dt_hyd = (1/freq)/(4*25)    # time-step - sec
t_save_hyd  = 2.e-06       # save time - sec0
t_stop_stg1 = 6.e-06       # stop time stg1 - sec 
t_stop_stg2 = 8.e-06       # stop time stg2 - sec

# Heat solver
cfl_heat = 0.1            # Courant number for diffusion equation cartesian (even dx, dy, dz)
dt_heat = cfl_heat*((xe-xb)/Nx)**2 / tdiff_host    # time-step - sec
t_save_heat = 1.           # save time - sec
t_stop_stg3 = 2.           # stop time - sec
t_stop_hifu_source = 1.    # stop hifu time - sec

t_step_stop_stg1 = int(t_stop_stg1*(c0/x0)/round(dt_hyd*c0/x0,6))
t_step_stop_stg2 = int(t_stop_stg2*(c0/x0)/round(dt_hyd*c0/x0,6))
t_step_save_hyd = int(t_save_hyd*(c0/x0)/round(dt_hyd*c0/x0,6))

t_step_stop_stg3 = int(t_stop_stg3*(c0/x0)/round(dt_heat*c0/x0,6))
t_step_save_heat = int(t_save_heat*(c0/x0)/round(dt_heat*c0/x0,6))
t_step_stop_source = int(t_stop_hifu_source*(c0/x0)/round(dt_heat*c0/x0,6))

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'cyl_coord'                    : 'F',
    'x_domain%beg'                 : xb/x0,
    'x_domain%end'                 : xe/x0,
    'y_domain%beg'                 : yb/x0,
    'y_domain%end'                 : ye/x0,
    'z_domain%beg'                 : zb/x0,
    'z_domain%end'                 : ze/x0,
    'stretch_x'                    : 'T',
    'a_x'                          : 40,
    'x_a'                          : -65e-03/x0,
    'x_b'                          : 65e-03/x0,
    'stretch_y'                    : 'T',
    'a_y'                          : 0.3,
    'y_a'                          : -16e-03/x0,
    'y_b'                          : 16e-03/x0,
    'stretch_z'                    : 'T',
    'a_z'                          : 0.3,
    'z_a'                          : -16e-03/x0,
    'z_b'                          : 16e-03/x0,
    'm'                            : Nx,
    'n'                            : Ny,
    'p'                            : Nz,
    'adap_dt'                      : 'T',       #Strang splitting
    'cfl_adap_dt'                  : 'T',
    'cfl_target'                   : 0.7,
    'dt'                           : round(dt_hyd*c0/x0,6),
    't_step_start'                 : 0,     ############### also modify stg1, stg2 or stg3 flags
    't_step_save'                  : t_step_save_hyd,       # Always stg1
    't_step_stop'                  : t_step_stop_stg1,      # Always stg1
    'n_start'                      : 0,     ## ADAPTIVE ### also modify stg1, stg2 or stg3 flags
    't_save'                       : t_save_hyd*(c0/x0),    # Always stg1
    't_stop'                       : t_stop_stg1*(c0/x0),   # Always stg1
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_fluids'                   : 3,         # Water/Phantom/BubbleGas
    'num_patches'                  : 2,
    'viscous'                      : 'T',
    'model_eqns'                   : 2,         # 5 model eqns
    'alt_soundspeed'               : 'F',       # Alternate sound speed (5 eqn model only)
    'mpp_lim'                      : 'T',       # Mixture physical parameters limits
    'mixture_err'                  : 'T',       # Mixture properties correction
    'time_stepper'                 : 3,         # O(3) TVD RK
    'weno_order'                   : 5,
    'weno_eps'                     : 1.0e-16,
    'mapped_weno'                  :'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     :-20,        # Acoustic input BC
    'bc_x%end'                     :-6,         # Nonreflective subsonic BC
    'bc_y%beg'                     :-2,         # Symmetric BC
    'bc_y%end'                     :-6,         # Nonreflective subsonic BC
    'bc_z%beg'                     :-2,         # Symmetric BC
    'bc_z%end'                     :-6,         # Nonreflective subsonic BC
    # ==========================================================

    # Acoustic source (bc == -20)===============================
    'acoustic_bc_params%iwave'      : 2, # transducer
    'acoustic_bc_params%ncycles'    : int(1e+6),
    'acoustic_bc_params%Pbase'      : patm/p0,
    'acoustic_bc_params%rho'        : rho_water/rho0,
    'acoustic_bc_params%cson'       : c_water/c0,
    'acoustic_bc_params%Pamp'       : pamplitude/p0,
    'acoustic_bc_params%freq'       : freq*x0/c0,
    'acoustic_bc_params%focLen'     : focLen/x0, 
    'acoustic_bc_params%focCal'     : 3.3e-3/x0,    # add calibration parameter bulb1: 3.3 mm -> focLen + focCal
    'acoustic_bc_params%apert'      : aperture/x0,
    # ==========================================================

    # HIFU parameters ==========================================
    'hifu'                          : 'T',
    'hifu_params%atmPres'           : patm/p0,
    'hifu_params%Tref'              : T_host/T0,
    # Automatize going from one stage to another
    'hifu_params%automatic_stages'  : 'T',
    # STG1: Develop hydrodynamic field
    'hifu_params%stg1'              : 'T',
    'hifu_params%t_stop_stg1'       : t_stop_stg1*(c0/x0),
    'hifu_params%t_step_stop_stg1'  : t_step_stop_stg1,
    # STG2: Sampling
    'hifu_params%stg2'              : 'T',
    'hifu_params%dt_stg2'           : round(dt_hyd*c0/x0,6),
    'hifu_params%t_stop_stg2'       : t_stop_stg2*(c0/x0),
    'hifu_params%t_step_stop_stg2'  : t_step_stop_stg2,
    # STG3: Solving heat equation
    'hifu_params%stg3'              : 'T',
    'hifu_params%intPrms'           : 'F', # True: Utilize qus from Prms
    'hifu_params%streaming'         : 'F', # True: Consider streaming velocity

    #'hifu_params%stg3_3d'           : 'F', # True: Transform 2D axisymmetric plane to 3D cylinder or 3D cartesian
    #'hifu_params%cartesian'         : 'F', # True: 3D cartesian via interpolation of the 2D plane
    # 'hifu_params%z_max'             : 2.0*math.pi, # Cylinder sector
    # 'hifu_params%p_cyl'             : 127,         # No. of cells in the azimuthal dir.
    # 'hifu_params%xb'                : xb_ht/x0, # Domain reduction x_beg
    # 'hifu_params%xe'                : xe_ht/x0, # Domain reduction x_end
    # 'hifu_params%ye'                : ye_ht/x0, # Domain reduction y_end & z_end
    # 'hifu_params%m'                 : int(Nx_ht), # No. cells x dir
    # 'hifu_params%n'                 : int(Ny_ht), # No. cells y dir
    # 'hifu_params%p'                 : int(Nz_ht), # No. cells z dir

    'hifu_params%dt_stg3'           : round(dt_heat*c0/x0,6),
    'hifu_params%t_step_save_stg3'  : t_step_save_heat,
    'hifu_params%t_step_stop_stg3'  : t_step_stop_stg3,
    'hifu_params%stepStopSource'    : t_step_stop_source,
    # ==========================================================

    # Lagrangian Bubbles ===========================
     'bubbles_lagrange'                 : 'T',
     'bubble_model'                     : 2,    # Keller-Miksis model
     'lag_params%nBubs_glb'             : 5,  # Number of bubbles
     'lag_params%solver_approach'       : 2,    # Two-way coupled
     'lag_params%cluster_type'          : 2,    # 1: p_inf from intepolation, 2: p_inf avg surrounding cells
     'lag_params%pressure_corrector'    : 'T',
     'lag_params%interaction_model'     : 1, # Interaction model: 1 -> kazuki & 2 -> Aditya's model
    #  'lag_params%influence'             : 3, # Number of surrounding cells to define influence volume
     'lag_params%smooth_type'           : 1,
     'lag_params%coatedBub_model'       : 'T',  # Marmmotant model
     'lag_params%heatTransfer_model'    : 'T',
     'lag_params%massTransfer_model'    : 'F',
     'lag_params%epsilonb'              : 1.0,
     'lag_params%valmaxvoid'            : 0.9,
     'lag_params%write_bubbles'         : 'T',
     'lag_params%write_bubbles_stats'   : 'F',
     'lag_params%c0'                    : c0,
     'lag_params%rho0'                  : rho0,
     'lag_params%T0'                    : T0,
     'lag_params%x0'                    : x0,
     'lag_params%diffcoefvap'           : diffVapor,
     'lag_params%Thost'                 : T_host,
     'lag_params%ss0_ctdBub'            : sigmaInit,
     'lag_params%srfElast_ctdBub'       : elasticity,
     'lag_params%srfDilVsc_ctdBub'      : dilatationalViscosity,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    'probe_wrt'                    :'F',
    # 'num_probes'                   : 3,
    # 'probe(1)%x'                   : focLen/x0,
    # 'probe(1)%y'                   : 0.,
    # 'probe(1)%z'                   : 0.,
    # 'probe(2)%x'                   : (focLen - 20.e-3)/x0,
    # 'probe(2)%y'                   : 0.,
    # 'probe(2)%z'                   : 0.,
    # 'probe(3)%x'                   : (1.e-3)/x0,
    # 'probe(3)%y'                   : 0.,
    # 'probe(3)%z'                   : 0.,
    # ==========================================================

    # Patch 1: Water (left) ====================================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0.5*(xe+xb)/x0,
    'patch_icpp(1)%y_centroid'     : 0.5*(ye+yb)/x0,
    'patch_icpp(1)%z_centroid'     : 0.5*(ze+zb)/x0,
    'patch_icpp(1)%length_x'       : 4*(xe-xb)/x0,
    'patch_icpp(1)%length_y'       : 4*(ye-yb)/x0,
    'patch_icpp(1)%length_z'       : 4*(ze-zb)/x0,
    'patch_icpp(1)%vel(1)'         : 0.,
    'patch_icpp(1)%vel(2)'         : 0.,
    'patch_icpp(1)%vel(3)'         : 0.,
    'patch_icpp(1)%pres'           : patm/p0,
    'patch_icpp(1)%alpha_rho(1)'   : rho_water/rho0,
    'patch_icpp(1)%alpha_rho(2)'   : 0.,
    'patch_icpp(1)%alpha_rho(3)'   : 0.,
    'patch_icpp(1)%alpha(1)'       : 1.,
    'patch_icpp(1)%alpha(2)'       : 0.,
    'patch_icpp(1)%alpha(3)'       : 0.,
    # ==========================================================

    # Patch 2: EMP (right) ====================================
    'patch_icpp(2)%geometry'        : 11,
    'patch_icpp(2)%alter_patch(1)'  : 'T',
    'patch_icpp(2)%smoothen'        : 'T',
    'patch_icpp(2)%smooth_patch_id' : 1,
    'patch_icpp(2)%smooth_coeff'    : 0.4,
    'patch_icpp(2)%x_centroid'      : 2.e-3/x0,
    'patch_icpp(2)%y_centroid'      : 1.e-3/x0,
    'patch_icpp(2)%z_centroid'      : 1.e-3/x0,
    'patch_icpp(2)%normal(1)'       : 1.0,
    'patch_icpp(2)%normal(2)'       : 0.0,
    'patch_icpp(2)%normal(3)'       : 0.0,
    'patch_icpp(2)%vel(1)'          : 0.,
    'patch_icpp(2)%vel(2)'          : 0.,
    'patch_icpp(2)%vel(3)'          : 0.,
    'patch_icpp(2)%pres'            : patm/p0,
    'patch_icpp(2)%alpha_rho(1)'    : 0.,
    'patch_icpp(2)%alpha_rho(2)'    : rho_host/rho0,
    'patch_icpp(2)%alpha_rho(3)'    : 0.,
    'patch_icpp(2)%alpha(1)'        : 0.,
    'patch_icpp(2)%alpha(2)'        : 1.,
    'patch_icpp(2)%alpha(3)'        : 0.,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    # Water
    'fluid_pp(1)%gamma'            : 1.0/(gamma_water-1.0),
    'fluid_pp(1)%pi_inf'           : gamma_water*(pi_inf_water/p0)/(gamma_water-1.0),
    'fluid_pp(1)%Re(1)'            : 1.0/(mu_water/(rho0*c0*x0)),
    'fluid_pp(1)%Re(2)'            : 1.0/(3*mu_water/(rho0*c0*x0)),
    'fluid_pp(1)%rho_cp'           : (rho_water/rho0)*(cp_water*(T0/(c0*c0))),
    'fluid_pp(1)%tdiff'            : tdiff_water/(x0*c0),
    'fluid_pp(1)%absCoef'          : abs_coef_water*x0,

    # EMP (host medium)
    'fluid_pp(2)%gamma'            : 1.0/(gamma_host-1.0),
    'fluid_pp(2)%pi_inf'           : gamma_host*(pi_inf_host/p0)/(gamma_host-1.0),
    'fluid_pp(2)%Re(1)'            : 1.0/(mu_host/(rho0*c0*x0)),
    'fluid_pp(2)%Re(2)'            : 1.0/(3*mu_host/(rho0*c0*x0)),
    'fluid_pp(2)%rho_cp'           : (rho_host/rho0)*(cp_host*(T0/(c0*c0))),
    'fluid_pp(2)%tdiff'            : tdiff_host/(x0*c0),
    'fluid_pp(2)%absCoef'          : abs_coef_host*x0,
    'fluid_pp(2)%mul0'             : mu_host,
    'fluid_pp(2)%ss'               : sigBubble,
    'fluid_pp(2)%pv'               : pv,
    'fluid_pp(2)%gamma_v'          : gamma_v,
    'fluid_pp(2)%M_v'              : MW_v,
    'fluid_pp(2)%k_v'              : k_v,
    'fluid_pp(2)%cp_v'             : cp_v,

    # Bubble gas state
    'fluid_pp(3)%gamma'            : 1./(gamma_g-1.),
    'fluid_pp(3)%pi_inf'           : 0.0E+00,
    'fluid_pp(3)%Re(1)'            : 1.0/(mu_g/(rho0*c0*x0)),
    'fluid_pp(3)%Re(2)'            : 1.0/(mu_g/(rho0*c0*x0)),
    'fluid_pp(3)%gamma_v'          : gamma_g,
    'fluid_pp(3)%M_v'              : MW_g,
    'fluid_pp(3)%k_v'              : k_g,
    'fluid_pp(3)%cp_v'             : cp_g,
    # ==========================================================
 }))

# ==============================================================================
