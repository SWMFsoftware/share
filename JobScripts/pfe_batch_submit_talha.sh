#!/bin/bash

# This script generates a PBS job file 'tmp.pfe', then submits it using a custom command.
#
# Usage:
#   ./generate_job.sh <NAME> <CPUs> [options]
#
# Options:
#   -long           : Adds a directive for the 'long' queue.
#   -t=<hours>      : Sets the walltime in hours (e.g., -t=24). Default is 8.
#   -no<model>      : Skips generating a line for a specific model (e.g., -noivy).

# --- Validation for Positional Arguments ---
if [ "$#" -lt 2 ] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
    echo "❌ Error: Invalid arguments."
    echo "Usage: $0 <NAME> <CPUs> [options]"
    exit 1
fi

# --- Positional Argument Assignment ---
JOB_NAME=$1
TARGET_CPUS=$2
shift 2 # Removes the first two arguments, so only flags remain.

# --- Default values & Initialization ---
WALLTIME_HOURS=8
QUEUE_LINE=""
declare -A EXCLUDE_MODELS

# --- Argument Parsing Loop for Optional Flags ---
for arg in "$@"; do
  case $arg in
    -long)
      QUEUE_LINE="#PBS -q long"
      ;;
    -t=*)
      val="${arg#*=}"
      if [[ "$val" =~ ^[0-9]+$ ]]; then
        WALLTIME_HOURS="$val"
      fi
      ;;
    -no*)
      model_to_exclude="${arg#-no}"
      EXCLUDE_MODELS["$model_to_exclude"]=1
      ;;
  esac
done

# --- Machine Configuration ---
declare -A models
models=( ["ivy"]=20 ["has"]=24 ["bro"]=28 ["bro_ele"]=28 ["sky_ele"]=40 ["cas_ait"]=40 ["rom_ait"]=128 )
MODEL_ORDER=("ivy" "has" "bro" "bro_ele" "sky_ele" "cas_ait" "rom_ait")

# --- Dynamically Build PBS Select Lines (Corrected Logic) ---
declare -a pbs_lines_array
is_first_active_line=true

for model in "${MODEL_ORDER[@]}"; do
  if [[ -v EXCLUDE_MODELS[$model] ]]; then
    continue
  fi

  ncpus=${models[$model]}
  select_val=$(( (TARGET_CPUS + ncpus - 1) / ncpus ))
  # CORRECTED: 'line' now only contains the directive itself, not the #PBS prefix.
  line="-l select=${select_val}:ncpus=${ncpus}:model=${model}"

  if $is_first_active_line; then
    # CORRECTED: Add the #PBS prefix here.
    pbs_lines_array+=("#PBS ${line}")
    is_first_active_line=false
  else
    # CORRECTED: Add the correct '### PBS' prefix here.
    pbs_lines_array+=("### PBS ${line}")
  fi
done
PBS_SELECT_LINES=$(printf "%s\n" "${pbs_lines_array[@]}")

# --- File Generation ---
cat << EOF > tmp.pfe
#!/bin/csh
########################################################
#PBS -S /bin/csh
#PBS -N SWMF
########################################################
${PBS_SELECT_LINES}
########################################################
${QUEUE_LINE}
#PBS -l walltime=${WALLTIME_HOURS}:00:00
#PBS -j oe
#PBS -m e
########################################################
cd \$PBS_O_WORKDIR
setenv MPI_TYPE_DEPTH 20
mpiexec ./SWMF.exe > runlog_\`date +%y%m%d%H%M\`
exit
if(! -f SWMF.SUCCESS) exit
if(-f SWMF.DONE) exit
./Restart.pl
qsub tmp.pfe
EOF
./qsub.pfe.pbspl.pl tmp.pfe "${JOB_NAME}"
rm -rf tmp.pfe
echo "✅ Requested ${TARGET_CPUS} CPUs on all queues"
