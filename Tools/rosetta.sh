#!/bin/bash
# ==============================================================================
#                 Rosetta Batch Design Script (Parallel Bash Version)
# ==============================================================================

# --- 用户配置区：请在这里修改您的路径和参数 ---
PDB_INPUT_DIRECTORY="/Rosetta/rosetta/pratice/TNfa/to3rosetta_inputs_pdbs_high/process" # 包含蛋白和多肽复合物的pdb文件
RESFILE_INPUT_DIRECTORY="/Rosetta/rosetta/pratice/TNfa/to3rosetta_resfiles_highMPNN"  # resfile文件的路径

OUTPUT_BASE_DIRECTORY="/Rosetta/rosetta/pratice/TNfa/rosetta3rd_outputs_high"   # resfile文件的路径

ROSETTA_XML_PROTOCOL="/Rosetta/rosetta/pratice/design_relax.xml"    # xml文件路径
ROSETTA_EXECUTABLE="/usr/local/src/rosetta_src_2021.16.61629_bundle/main/source/build/src/release/linux/4.15/64/x86/gcc/7/mpi/rosetta_scripts.mpi.linuxgccrelease" # 无需更改
ROSETTA_DATABASE="/usr/local/src/rosetta_src_2021.16.61629_bundle/main/database"  # 无需更改

NUM_RESFILES_PER_PDB=5 # 告知脚本每个复合物一共有几个.resfile文件，例如上图，我们可以输入5

NUM_PARALLEL_PROCESSES=6 # 告知脚本使用几个cpu

PRINT_ROSETTA_VERBOSE_OUTPUT="true" 
TEMP_COMMAND_LIST="rosetta_commands_to_run.tmp"
LOG_FILE="rosetta_batch_log_$(date +%Y%m%d_%H%M%S).log"
ERROR_LOG_FILE="rosetta_batch_errors_$(date +%Y%m%d_%H%M%S).log"

log_message() {
    local message="$1"
    echo -e "$message" | tee -a "$LOG_FILE"
}

log_error() {
    local message="$1"
    echo -e "ERROR: $message" | tee -a "$LOG_FILE" "$ERROR_LOG_FILE" >&2
}

log_message "--- Starting Rosetta Batch Design Script (Parallel Bash Version) ---"
log_message "PDB Input Directory: $PDB_INPUT_DIRECTORY"
log_message "Resfile Input Directory: $RESFILE_INPUT_DIRECTORY"
log_message "Output Base Directory: $OUTPUT_BASE_DIRECTORY"
log_message "Rosetta XML Protocol: $ROSETTA_XML_PROTOCOL"
log_message "Rosetta Executable: $ROSETTA_EXECUTABLE"
log_message "Rosetta Database: $ROSETTA_DATABASE"
log_message "Number of Resfiles per PDB: $NUM_RESFILES_PER_PDB"
log_message "Parallel Processes: $NUM_PARALLEL_PROCESSES"
log_message "Print Rosetta Verbose Output: $PRINT_ROSETTA_VERBOSE_OUTPUT"
log_message "Log File: $LOG_FILE"
log_message "Error Log File: $ERROR_LOG_FILE"
log_message "--------------------------------------------------"

log_message "--- Running Pre-checks ---"

if [ ! -f "$ROSETTA_EXECUTABLE" ]; then
    log_error "Rosetta Executable '$ROSETTA_EXECUTABLE' not found! Please check path."
    exit 1
fi
log_message "Rosetta Executable: OK"

if [ ! -d "$ROSETTA_DATABASE" ]; then
    log_error "Rosetta Database '$ROSETTA_DATABASE' not found! Please check path."
    exit 1
fi
log_message "Rosetta Database: OK"

if [ ! -f "$ROSETTA_XML_PROTOCOL" ]; then
    log_error "Rosetta XML Protocol File '$ROSETTA_XML_PROTOCOL' not found! Please check path."
    exit 1
fi
log_message "Rosetta XML Protocol: OK"

if [ ! -d "$PDB_INPUT_DIRECTORY" ]; then
    log_error "PDB Input Directory '$PDB_INPUT_DIRECTORY' not found! Please check path."
    exit 1
fi
log_message "PDB Input Directory: OK"

if [ ! -d "$RESFILE_INPUT_DIRECTORY" ]; then
    log_error "Resfile Input Directory '$RESFILE_INPUT_DIRECTORY' not found! Please check path."
    exit 1
fi
log_message "Resfile Input Directory: OK"

mkdir -p "$OUTPUT_BASE_DIRECTORY"
log_message "Output Base Directory '$OUTPUT_BASE_DIRECTORY' ready."
log_message "--- Checks Complete, Starting Process ---"
log_message "--------------------------------------------------"

log_message "Building task list..."
> "$TEMP_COMMAND_LIST"

PDB_FILES=$(find "$PDB_INPUT_DIRECTORY" -maxdepth 1 -type f -name "*.pdb")
if [ -z "$PDB_FILES" ]; then
    log_message "Warning: No .pdb files found in '$PDB_INPUT_DIRECTORY'."
    exit 0
fi

TOTAL_PDB_FILES=$(echo "$PDB_FILES" | wc -l)
log_message "Found $TOTAL_PDB_FILES PDB files. Each PDB will run $NUM_RESFILES_PER_PDB designs."

TASK_COUNTER=0
for PDB_PATH in $PDB_FILES; do
    PDB_FILENAME=$(basename "$PDB_PATH")
    PDB_BASENAME="${PDB_FILENAME%.*}" 
    
    CURRENT_OUTPUT_DIR="${OUTPUT_BASE_DIRECTORY}/${PDB_BASENAME}_designs"
    mkdir -p "$CURRENT_OUTPUT_DIR"

    for ((i=1; i<=NUM_RESFILES_PER_PDB; i++)); do
        RESFILE_NAME="${PDB_BASENAME}_${i}.resfile"
        RESFILE_PATH="${RESFILE_INPUT_DIRECTORY}/${RESFILE_NAME}"

        if [ ! -f "$RESFILE_PATH" ]; then
            log_message "  --> Warning (Task Build Phase): Resfile '$RESFILE_NAME' not found. Skipping this combination."
            continue
        fi

        OUTPUT_PDB_PREFIX="${PDB_BASENAME}_design_${i}_"
        ROSETTA_COMMAND=(
            "$ROSETTA_EXECUTABLE"
            "-database" "$ROSETTA_DATABASE"
            "-in:file:s" "$PDB_PATH"
            "-parser:protocol" "$ROSETTA_XML_PROTOCOL"
            "-parser:script_vars" "resfile_name=${RESFILE_PATH}"
            "-nstruct" "1"
            "-beta"
            "-ex1" "-ex2aro"
            "-out:path:all" "$CURRENT_OUTPUT_DIR"
            "-out:prefix" "$OUTPUT_PDB_PREFIX"
            "-ignore_zero_occupancy" "false"
        )
        


        if [ "$PRINT_ROSETTA_VERBOSE_OUTPUT" = "true" ]; then
            CMD_STRING="${ROSETTA_COMMAND[@]}"
            TEMP_TASK_LOG="${CURRENT_OUTPUT_DIR}/${PDB_BASENAME}_${i}_rosetta.log"
            echo "(time ${CMD_STRING} > ${TEMP_TASK_LOG} 2>&1) && echo \"  --> Design Success: ${PDB_FILENAME} using ${RESFILE_NAME}\" || echo \"  --> Design Failed: ${PDB_FILENAME} using ${RESFILE_NAME}\" && cat ${TEMP_TASK_LOG} >> \"$LOG_FILE\"" >> "$TEMP_COMMAND_LIST"
        else
            CMD_STRING="${ROSETTA_COMMAND[@]} > /dev/null 2>&1"
            echo "(${CMD_STRING}) && echo \"  --> Design Success: ${PDB_FILENAME} using ${RESFILE_NAME}\" || echo \"  --> Design Failed: ${PDB_FILENAME} using ${RESFILE_NAME}\"" >> "$TEMP_COMMAND_LIST"
        fi

        TASK_COUNTER=$((TASK_COUNTER + 1))
    done
done

TOTAL_EXECUTABLE_TASKS=$(wc -l < "$TEMP_COMMAND_LIST")

if [ "$TOTAL_EXECUTABLE_TASKS" -eq 0 ]; then
    log_message "No executable tasks found (please check PDBs and Resfiles and naming conventions). Script exiting."
    rm -f "$TEMP_COMMAND_LIST"
    exit 0
fi

log_message "Task list built. Total executable tasks: $TOTAL_EXECUTABLE_TASKS. Will run in $NUM_PARALLEL_PROCESSES parallel processes."
log_message "--------------------------------------------------"

OVERALL_START_TIME=$(date +%s)

log_message "Starting parallel execution..."

cat "$TEMP_COMMAND_LIST" | xargs -P "$NUM_PARALLEL_PROCESSES" -I {} bash -c '{}'

OVERALL_END_TIME=$(date +%s)
TOTAL_RUN_TIME=$((OVERALL_END_TIME - OVERALL_START_TIME))

HOURS=$((TOTAL_RUN_TIME / 3600))
MINUTES=$(( (TOTAL_RUN_TIME % 3600) / 60 ))
SECONDS_FLOAT=$((TOTAL_RUN_TIME % 60))

log_message "\n======================================================================"
log_message "All PDB and Resfile design tasks completed!"
log_message "Total elapsed time: ${HOURS}h ${MINUTES}min ${SECONDS_FLOAT}.0s" 
log_message "All results are saved in '$OUTPUT_BASE_DIRECTORY' directory."
log_message "======================================================================"

rm -f "$TEMP_COMMAND_LIST"

log_message "Script execution finished."