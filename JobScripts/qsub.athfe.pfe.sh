#!/usr/bin/env bash
CURRENT_HOST="$(hostname -s)"
echo
echo
echo '########### NAS Job Script Called ###########'
echo
case "$CURRENT_HOST" in
    pfe*)
      	echo '################ PFE Detected ###############'
      	echo
        if [[ " $* " =~ [[:space:]]-i[[:space:]] ]]; then
            ############################ PFE INTERACTIVE STARTS ##############################
            orig_dir=$(pwd)
            cd "$HOME" || exit 1
            table=$(node_stats.sh | awk '
                /^Nodes used\/free by hardware type:/ { in_hw=1; next }
                /^GPUs used\/free/ { in_hw=0 }
                in_hw && /Cores/ {
                    nodename=$1
                    for (i=1; i<=NF; i++) {
                        if ($i=="Cores") {
                            cores[nodename]=$(i-1)
                            break
                        }
                    }
                }
                /^Nodes currently allocated to the devel queue:/ {
                    in_block=1
                    next
                }
                in_block==1 && NF==0 { exit }
                in_block==1 {
                    nodename=$1
                    free=$NF
                    data[nodename, "cores"]=cores[nodename]
                    data[nodename, "free"]=free
                    data[nodename, "total"]=cores[nodename]*free
                    nodes[++count]=nodename
                }
                END {
                    max_node_len=length("NODETYPE")
                    max_cpu_len=length("CPUs / NODE")
                    max_free_len=length("FREE NODES")
                    max_total_len=length("TOTAL CPUs")
                    for (i=1; i<=count; i++) {
                        n = nodes[i]
                        if (length(n) > max_node_len) max_node_len=length(n)
                        if (length(data[n,"cores"]) > max_cpu_len) max_cpu_len=length(data[n,"cores"])
                        if (length(data[n,"free"]) > max_free_len) max_free_len=length(data[n,"free"])
                        if (length(data[n,"total"]) > max_total_len) max_total_len=length(data[n,"total"])
                    }
                    printf "%-*s | %-*s | %-*s | %-*s\n", max_node_len, "NODETYPE", max_cpu_len, "CPUs / NODE", max_free_len, "FREE NODES", max_total_len, "FREE CPUs"
                    printf "%s\n", gensub(/./,"-","g",sprintf("%*s", max_node_len+max_cpu_len+max_free_len+max_total_len+9,""))
                    for (i=1; i<=count; i++) {
                        n = nodes[i]
                        printf "%-*s | %-*d | %-*s | %-*d\n", max_node_len, n, max_cpu_len, data[n,"cores"], max_free_len, data[n,"free"], max_total_len, data[n,"total"]
                    }
                }')
            echo "$table"
            read -p "Please enter the NODETYPE: " nodetype
            read -p "How many nodes would you like to request? " nodereq
            ncpus=$(echo "$table" | awk -F'|' -v nt="$nodetype" '{
                gsub(/^[ \t]+|[ \t]+$/, "", $1)
                gsub(/^[ \t]+|[ \t]+$/, "", $2)
                if ($1==nt) print $2
            }')
            qsub -I -V -X -q devel -lselect=${nodereq}:ncpus=${ncpus}:model=${nodetype},walltime=02:00:00
            ############################# PFE INTERACTIVE ENDS ###############################
        else
            ############################ PFE NON-INTERACTIVE STARTS ##############################
            rm -f job.pfe.new
            cancel_check() {
                if [ $? -ne 0 ]; then
                    echo "Cancelled by user."
                    exit 1
                fi
            }
            declare -A models=(
                ["bro_ele"]=28
                ["cas_ait"]=40
                ["mil_ait"]=128
                ["rom_ait"]=128
                ["sky_ele"]=40
            )
            MODEL_ORDER=("bro_ele" "cas_ait" "mil_ait" "rom_ait" "sky_ele")
            declare -A previous_selected
            if [ -f job.pfe ]; then
                while read -r line; do
                    if [[ "$line" =~ select=.*:model= ]]; then
                        model=$(echo "$line" | sed -E 's/.*model=([^ :]+).*/\1/')
                        previous_selected[$model]=1
                    fi
                done < job.pfe
            fi
            JOB_NAME=$(whiptail --inputbox "Enter job name (leave empty for SWMF):" 8 60 "" 3>&1 1>&2 2>&3)
            cancel_check
            [ -z "$JOB_NAME" ] && JOB_NAME="SWMF"
            while true; do
                TARGET_CPUS=$(whiptail --inputbox "Enter total number of CPUs to request:" 8 60 "" 3>&1 1>&2 2>&3)
                cancel_check
                if [[ "$TARGET_CPUS" =~ ^[0-9]+$ ]] && [ "$TARGET_CPUS" -gt 0 ]; then
                    break
                fi
                whiptail --msgbox "Please enter a valid positive integer for CPUs." 8 50
            done
            WALLTIME_HOURS=$(whiptail --inputbox "Enter walltime (hours, default 8):" 8 60 "" 3>&1 1>&2 2>&3)
            cancel_check
            if ! [[ "$WALLTIME_HOURS" =~ ^[0-9]+$ ]] || [ -z "$WALLTIME_HOURS" ]; then
                WALLTIME_HOURS=8
            fi
            if [ "$WALLTIME_HOURS" -gt 8 ]; then
                QUEUE_LINE="#PBS -q long"
            else
                QUEUE_LINE="#PBS -q normal"
            fi
            checklist_items=()
            max_len=0
            for model in "${MODEL_ORDER[@]}"; do
                desc="${models[$model]} CPUs/node"
                if [ ${#previous_selected[@]} -eq 0 ] || [ -n "${previous_selected[$model]}" ]; then
                    checklist_items+=("$model" "$desc" "ON")
                else
                    checklist_items+=("$model" "$desc" "OFF")
                fi
                item_len=$(( ${#model} + ${#desc} + 1 ))
                (( item_len > max_len )) && max_len=$item_len
            done
            whiptail_width=$((max_len + 10))
            [ "$whiptail_width" -lt 40 ] && whiptail_width=40
            [ "$whiptail_width" -gt 100 ] && whiptail_width=100
            num_models=${#MODEL_ORDER[@]}
            whiptail_height=$((num_models + 5))
            MODEL_SELECTION=$(whiptail --title "Select Models" --checklist \
                "Select which models to include (Space = toggle, Enter = confirm):" \
                $whiptail_height $whiptail_width $num_models \
                "${checklist_items[@]}" \
                3>&1 1>&2 2>&3)
            cancel_check
            MODEL_SELECTION=$(echo "$MODEL_SELECTION" | tr -d '"')
            declare -a pbs_lines_array
            is_first_active_line=true
            for model in ${MODEL_SELECTION}; do
                ncpus=${models[$model]}
                select_val=$(( (TARGET_CPUS + ncpus - 1) / ncpus ))
                line="-l select=${select_val}:ncpus=${ncpus}:model=${model}"
                if $is_first_active_line; then
                    pbs_lines_array+=("#PBS ${line}")
                    is_first_active_line=false
                else
                    pbs_lines_array+=("### PBS ${line}")
                fi
            done
            PBS_SELECT_LINES=$(printf "%s\n" "${pbs_lines_array[@]}")
            cat << EOF > job.pfe
#!/bin/csh
########################################################
#PBS -S /bin/csh
#PBS -N ${JOB_NAME}
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
qsub job.pfe
EOF
            if [ -x ./qsub.pfe.pbspl.pl ]; then
                ./qsub.pfe.pbspl.pl job.pfe "${JOB_NAME}"
            fi
            whiptail --msgbox "Job '${JOB_NAME}' submitted requesting ${TARGET_CPUS} CPUs (Walltime:${WALLTIME_HOURS}h)" 8 70
            ############################# PFE NON-INTERACTIVE ENDS ###############################
        fi
        ;;
    ath*)
    	echo '############### ATHENA Detected #############'
    	echo
        if [[ " $* " =~ [[:space:]]-i[[:space:]] ]]; then
            ############################ ATHENA INTERACTIVE STARTS ##############################
            read -rp "How many interactive nodes? Enter 1 or 2: " NODES
            if [[ "$NODES" != "1" && "$NODES" != "2" ]]; then
                echo "Invalid input: You entered '$NODES'. Please enter only 1 or 2. Cancelling."
                exit 1
            fi
            qsub -I -q devel -l select="${NODES}":ncpus=256:model=tur_ath -l walltime=2:00:00
            ############################# ATHENA INTERACTIVE ENDS ###############################
        else
            ############################ ATHENA NON-INTERACTIVE STARTS ##############################
            read -rp "Enter job name: " JOB_NAME
            read -rp "Enter total number of cores: " CORES
            read -rp "Enter walltime (in hours): " HOURS
            if ! [[ "$HOURS" =~ ^[0-9]+$ ]]; then
                echo "Error: Please enter walltime as an integer number of hours."
                exit 1
            fi
            WALLTIME=$(printf "%02d:00:00" "$HOURS")
            if [ "$HOURS" -gt 8 ]; then
                QUEUE="long"
            else
                QUEUE="normal"
            fi
            NODES=$(( (CORES + 255) / 256 ))
            CORES=$(( NODES * 256 ))
            SCRIPT_FILE="tjob.athfe"
            cat << EOF > "$SCRIPT_FILE"
#!/bin/csh
#PBS -N ${JOB_NAME}
#PBS -l select=${NODES}:ncpus=256:mpiprocs=256:model=tur_ath
#PBS -q ${QUEUE}
#PBS -l walltime=${WALLTIME}
#PBS -j oe
#PBS -m e
module purge
module load PrgEnv-intel cray-pals idl
setenv FI_PROVIDER cxi
cd \$PBS_O_WORKDIR
mpiexec -n ${CORES} ./SWMF.exe > runlog_\`date +%y%m%d%H%M\`
EOF
            echo "Created ${SCRIPT_FILE} (Nodes:${NODES}, Queue: ${QUEUE}, Walltime:${WALLTIME})"
            qsub "$SCRIPT_FILE"
            ############################# ATHENA NON-INTERACTIVE ENDS ###############################
        fi
        ;;
    *)
        echo "Host not recognized" >&2
        exit 1
        ;;
esac
