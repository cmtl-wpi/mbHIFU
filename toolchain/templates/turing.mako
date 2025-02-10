#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --job-name="${name}"
#SBATCH --time=${walltime}
#SBATCH --error="${name}.err"
#SBATCH --output="${name}.out"
% if account:
#SBATCH --account=${account}
% endif
% if partition:
#SBATCH --partition=${partition}
% endif
% if quality_of_service:
#SBATCH --qos=${quality_of_service}
% endif
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type="BEGIN, END, FAIL"	
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
% if engine == 'batch':
	. ./mfc.sh load -c t -m ${'g' if gpu else 'c'}
% endif
cd - > /dev/null
echo

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        (set -x; ${profiler}                              \
	    srun --mpi=pmi2     \
                     ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                     "${target.get_install_binpath(case)}")
    % endif


    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
