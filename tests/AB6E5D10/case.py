#!/usr/bin/env python3
#
# /home/dgvacarevelo/MFC/tests/AB6E5D10/case.py:
# 3D -> mixlayer_perturb -> bubbles

import json
import argparse

parser = argparse.ArgumentParser(
    prog="/home/dgvacarevelo/MFC/tests/AB6E5D10/case.py",
    description="/home/dgvacarevelo/MFC/tests/AB6E5D10/case.py: 3D -> mixlayer_perturb -> bubbles",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("dict", type=str, metavar="DICT", help=argparse.SUPPRESS)

ARGS = vars(parser.parse_args())

ARGS["dict"] = json.loads(ARGS["dict"])

case = {
    "run_time_info": "T",
    "m": 24,
    "n": 35,
    "p": 24,
    "dt": 1e-06,
    "t_step_start": 0,
    "t_step_stop": 50,
    "t_step_save": 50,
    "num_patches": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "F",
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 0.99999,
    "patch_icpp(1)%alpha(1)": 1e-05,
    "patch_icpp(2)%pres": 0.5,
    "patch_icpp(2)%alpha_rho(1)": 0.5,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(3)%pres": 0.1,
    "patch_icpp(3)%alpha_rho(1)": 0.125,
    "patch_icpp(3)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 0.16393442623,
    "fluid_pp(1)%pi_inf": 22.312399959394575,
    "fluid_pp(1)%cv": 0.0,
    "fluid_pp(1)%qv": 0.0,
    "fluid_pp(1)%qvp": 0.0,
    "bubbles": "T",
    "Ca": 0.7160271976687712,
    "Web": 5.660481099656358,
    "Re_inv": 0.0069829599021229965,
    "pref": 101325.0,
    "rhoref": 1000.0,
    "bubble_model": 3,
    "polytropic": "T",
    "polydisperse": "F",
    "thermal": 1,
    "R0ref": 1e-05,
    "patch_icpp(1)%r0": 1,
    "patch_icpp(1)%v0": 0,
    "patch_icpp(2)%r0": -1000000.0,
    "patch_icpp(2)%v0": -1000000.0,
    "patch_icpp(3)%r0": -1000000.0,
    "patch_icpp(3)%v0": -1000000.0,
    "qbmm": "F",
    "dist_type": 2,
    "poly_sigma": 0.3,
    "R0_type": 1,
    "sigR": 0.1,
    "sigV": 0.1,
    "rhoRV": 0.0,
    "acoustic_source": "F",
    "num_source": 1,
    "acoustic(1)%loc(1)": 0.5,
    "acoustic(1)%mag": 0.2,
    "acoustic(1)%length": 0.25,
    "acoustic(1)%dir": 1.0,
    "acoustic(1)%npulse": 1,
    "acoustic(1)%pulse": 1,
    "rdma_mpi": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": 360.0,
    "y_domain%beg": -180.0,
    "y_domain%end": 180.0,
    "z_domain%beg": 0.0,
    "z_domain%end": 360.0,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -6,
    "bc_y%end": -6,
    "bc_z%beg": -3,
    "bc_z%end": -3,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%z_centroid": 180.0,
    "patch_icpp(1)%length_z": 360.0,
    "patch_icpp(2)%z_centroid": -1000000.0,
    "patch_icpp(2)%length_z": -1000000.0,
    "patch_icpp(3)%z_centroid": -1000000.0,
    "patch_icpp(3)%length_z": -1000000.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_y": 360.0,
    "patch_icpp(1)%x_centroid": 180.0,
    "patch_icpp(1)%length_x": 360.0,
    "patch_icpp(1)%vel(1)": 1.1966855884162177,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(2)%geometry": -100,
    "patch_icpp(2)%y_centroid": -1000000.0,
    "patch_icpp(2)%length_y": -1000000.0,
    "patch_icpp(2)%x_centroid": -1000000.0,
    "patch_icpp(2)%length_x": -1000000.0,
    "patch_icpp(2)%vel(1)": -1000000.0,
    "patch_icpp(2)%vel(2)": -1000000.0,
    "patch_icpp(2)%vel(3)": -1000000.0,
    "patch_icpp(3)%geometry": -100,
    "patch_icpp(3)%y_centroid": -1000000.0,
    "patch_icpp(3)%length_y": -1000000.0,
    "patch_icpp(3)%x_centroid": -1000000.0,
    "patch_icpp(3)%length_x": -1000000.0,
    "patch_icpp(3)%vel(1)": -1000000.0,
    "patch_icpp(3)%vel(2)": -1000000.0,
    "patch_icpp(3)%vel(3)": -1000000.0,
    "mixlayer_vel_profile": "T",
    "mixlayer_domain": 1.475,
    "mixlayer_vel_coef": 0.6,
    "mixlayer_perturb": "T",
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "fluid_pp(1)%Re(1)": 1.6881644098979287,
    "fluid_pp(2)%gamma": 2.5,
    "fluid_pp(2)%pi_inf": 0.0,
    "adv_n": "T",
    "nb": 1
}
mods = {}

if "post_process" in ARGS["dict"]["targets"]:
    mods = {
        'parallel_io'  : 'T', 'cons_vars_wrt'   : 'T',
        'prim_vars_wrt': 'T', 'alpha_rho_wrt(1)': 'T',
        'rho_wrt'      : 'T', 'mom_wrt(1)'      : 'T',
        'vel_wrt(1)'   : 'T', 'E_wrt'           : 'T',
        'pres_wrt'     : 'T', 'alpha_wrt(1)'    : 'T',
        'gamma_wrt'    : 'T', 'heat_ratio_wrt'  : 'T',
        'pi_inf_wrt'   : 'T', 'pres_inf_wrt'    : 'T',
        'c_wrt'        : 'T',
    }

    if case['p'] != 0:
        mods['fd_order']  = 1
        mods['omega_wrt(1)'] = 'T'
        mods['omega_wrt(2)'] = 'T'
        mods['omega_wrt(3)'] = 'T'

print(json.dumps({**case, **mods}))
