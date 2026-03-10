#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --job-name="${name}"
##SBATCH --time=${walltime}
#SBATCH --mail-user=dgvacarevelo@wpi.edu
#SBATCH --mail-type=BEGIN
% if partition:
##SBATCH --partition=${partition}
#SBATCH --time=1-00:00:00
% if gpu:
#SBATCH --partition=gpuA100x4
% else:
#SBATCH --partition=cpu
% endif
% else:
#SBATCH --time=0-01:00:00
% if gpu:
#SBATCH --partition=gpuA100x4-interactive
% else:
#SBATCH --partition=cpu-interactive
% endif
% endif
% if account:
#SBATCH --account="${account}"
% endif
% if gpu:
#SBATCH --gpus-per-node=4
#SBATCH --mem=208G
#SBATCH --gpu-bind=closest
#SBATCH --account=bgko-delta-gpu
###SBATCH --partition=gpuA100x4-interactive
% else:
#SBATCH --account=bgko-delta-cpu
###SBATCH --partition=cpu-interactive
% endif
#SBATCH --output="${name}.out"
#SBATCH --error="${name}.err"
#SBATCH --export=ALL
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type="BEGIN, END, FAIL"
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
. ./mfc.sh load -c d -m ${'g' if gpu else 'c'}
cd - > /dev/null
echo

# % if gpu:
#     export MPICH_GPU_SUPPORT_ENABLED=1 # Disable GPU-Direct MPI
# % endif

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
    	% if gpu:
            (set -x; ${profiler}                         \
                    mpirun  -x UCX_NET_DEVICES=all       \
                            -mca coll_hcoll_enable 0     \
                            -np ${nodes*tasks_per_node}  \
                            "${target.get_install_binpath(case)}")
	    % else:
            % if partition:
                (set -x; ${profiler}                            \
                    srun    --account=bgko-delta-cpu            \
                            --partition=cpu                     \
                            --ntasks=${nodes*tasks_per_node}    \
                            "${target.get_install_binpath(case)}")
            % else:
                (set -x; ${profiler}                            \
                    srun    --account=bgko-delta-cpu            \
                            --partition=cpu-interactive         \
                            --ntasks=${nodes*tasks_per_node}    \
                            "${target.get_install_binpath(case)}")
            % endif
	    % endif
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
