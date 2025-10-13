#!/bin/bash
#PBS -N pancreatic_tumor_sweep
#PBS -t 0-9
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -l walltime=02:00:00
#PBS -o logs/job_${PBS_JOBID}.out
#PBS -e logs/job_${PBS_JOBID}.err
#PBS -j oe

# PBS Job Array Script for Pancreatic Tumor Model Parameter Sweep
#
# Usage:
#   1. Generate parameter files first:
#      python3 generate_param_sweep.py --sweep-type tumor_growth --num-jobs 10
#   2. Submit job array:
#      qsub submit_param_sweep_pbs.sh
#
# This script will:
# - Run one simulation per array task
# - Use different parameter files for each task
# - Save outputs in separate directories

# Load required modules (adjust for your HPC system)
# module load gcc/9.3.0
# module load cmake/3.18.0
# module load python/3.8

# Set up environment
export OMP_NUM_THREADS=1

# Base directory (adjust to your path)
BASE_DIR="${PBS_O_WORKDIR}"
PARAM_SWEEP_DIR="${BASE_DIR}/param_sweep"
RESULTS_DIR="${BASE_DIR}/results"

# Create necessary directories
mkdir -p logs
mkdir -p "${RESULTS_DIR}"

echo "============================================"
echo "PBS Job Array Task: ${PBS_ARRAYID}"
echo "Job ID: ${PBS_JOBID}"
echo "Working Directory: ${PBS_O_WORKDIR}"
echo "Node: ${PBS_O_HOST}"
echo "Start Time: $(date)"
echo "============================================"

# Change to working directory
cd "${BASE_DIR}"

# Get job configuration
if [ ! -f "${PARAM_SWEEP_DIR}/job_configs.json" ]; then
    echo "Error: job_configs.json not found. Run generate_param_sweep.py first."
    exit 1
fi

# Extract job information for this array task
PYTHON_CMD=$(cat << 'EOF'
import json
import sys
import os

array_id = int(os.environ.get('PBS_ARRAYID', '0'))
with open('param_sweep/job_configs.json', 'r') as f:
    config = json.load(f)

if array_id >= len(config['jobs']):
    print(f"Error: Array ID {array_id} exceeds number of jobs {len(config['jobs'])}")
    sys.exit(1)

job = config['jobs'][array_id]
print(f"JOB_NAME={job['job_name']}")
print(f"PARAM_FILE={job['param_file']}")
print(f"DESCRIPTION={job['description']}")
EOF
)

# Get job parameters
eval $(python3 -c "$PYTHON_CMD")

# Check if job parameters were set
if [ -z "$JOB_NAME" ]; then
    echo "Error: Failed to extract job parameters"
    exit 1
fi

# Create output directory for this job
OUTPUT_DIR="${RESULTS_DIR}/${JOB_NAME}"
mkdir -p "${OUTPUT_DIR}"

# Copy parameter file to output directory for reference
cp "${PARAM_FILE}" "${OUTPUT_DIR}/params.txt"

echo "Job Name: ${JOB_NAME}"
echo "Parameter File: ${PARAM_FILE}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Description: ${DESCRIPTION}"
echo "============================================"

# Run the simulation
echo "Starting simulation..."
export PARAM_FILE="${PARAM_FILE}"
export OUTPUT_DIR="${OUTPUT_DIR}"

# Use timeout to limit simulation time (PBS walltime should handle this, but just in case)
timeout 1h ./build/pancreatic_tumor_new --param-file "${PARAM_FILE}" --output-dir "${OUTPUT_DIR}"
EXIT_CODE=$?

echo "============================================"
echo "Simulation completed with exit code: ${EXIT_CODE}"
echo "End Time: $(date)"

# Check if simulation was successful
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Simulation completed successfully"
    
    # Create a summary file
    echo "Job: ${JOB_NAME}" > "${OUTPUT_DIR}/job_summary.txt"
    echo "Description: ${DESCRIPTION}" >> "${OUTPUT_DIR}/job_summary.txt"
    echo "Parameter File: ${PARAM_FILE}" >> "${OUTPUT_DIR}/job_summary.txt"
    echo "PBS Job ID: ${PBS_JOBID}" >> "${OUTPUT_DIR}/job_summary.txt"
    echo "Start Time: $(date)" >> "${OUTPUT_DIR}/job_summary.txt"
    echo "End Time: $(date)" >> "${OUTPUT_DIR}/job_summary.txt"
    echo "Exit Code: ${EXIT_CODE}" >> "${OUTPUT_DIR}/job_summary.txt"
    echo "Node: ${PBS_O_HOST}" >> "${OUTPUT_DIR}/job_summary.txt"
    
    # Check output files
    if [ -f "${OUTPUT_DIR}/populations.csv" ]; then
        wc -l "${OUTPUT_DIR}/populations.csv" >> "${OUTPUT_DIR}/job_summary.txt"
        echo "Output file created: populations.csv" >> "${OUTPUT_DIR}/job_summary.txt"
    else
        echo "WARNING: populations.csv not found" >> "${OUTPUT_DIR}/job_summary.txt"
    fi
    
    # List all generated files
    echo "Generated files:" >> "${OUTPUT_DIR}/job_summary.txt"
    ls -la "${OUTPUT_DIR}" >> "${OUTPUT_DIR}/job_summary.txt"
    
elif [ $EXIT_CODE -eq 124 ]; then
    echo "TIMEOUT: Simulation exceeded time limit"
    echo "Status: TIMEOUT" >> "${OUTPUT_DIR}/job_summary.txt"
else
    echo "ERROR: Simulation failed with exit code ${EXIT_CODE}"
    echo "Status: FAILED" >> "${OUTPUT_DIR}/job_summary.txt"
fi

echo "============================================"