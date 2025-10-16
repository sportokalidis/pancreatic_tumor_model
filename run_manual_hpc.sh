#!/bin/bash
# Manual HPC Execution Script
# For HPC systems without job schedulers or when you have direct access to compute nodes
#
# Usage: ./run_manual_hpc.sh [sweep_type] [num_jobs] [nodes_file]

set -e

SWEEP_TYPE=${1:-"tumor_growth"}
NUM_JOBS=${2:-10}
NODES_FILE=${3:-"nodes.txt"}

echo "============================================"
echo "Manual HPC Parameter Sweep"
echo "============================================"
echo "Sweep Type: $SWEEP_TYPE"
echo "Number of Jobs: $NUM_JOBS"
echo "Nodes File: $NODES_FILE"
echo "============================================"

# Step 1: Generate parameter files
echo "Step 1: Generating parameter files..."
python3 generate_param_sweep.py \
    --sweep-type "$SWEEP_TYPE" \
    --num-jobs "$NUM_JOBS" \
    --output-dir param_sweep

if [ ! -f "param_sweep/job_configs.json" ]; then
    echo "Error: Parameter generation failed"
    exit 1
fi

echo "✓ Generated $NUM_JOBS parameter files"

# Step 2: Check nodes file
if [ ! -f "$NODES_FILE" ]; then
    echo "Creating example nodes file: $NODES_FILE"
    cat > "$NODES_FILE" << EOF
# List your compute nodes here (one per line)
# Examples:
# node001
# node002  
# node003
# compute-0-1
# compute-0-2
localhost
EOF
    echo "Please edit $NODES_FILE with your actual compute node names"
    echo "Then run this script again"
    exit 1
fi

# Read available nodes
NODES=($(grep -v '^#' "$NODES_FILE" | grep -v '^$'))
NUM_NODES=${#NODES[@]}

if [ $NUM_NODES -eq 0 ]; then
    echo "Error: No nodes found in $NODES_FILE"
    exit 1
fi

echo "Available nodes: ${NODES[@]}"
echo "Number of nodes: $NUM_NODES"

# Step 3: Create directories
mkdir -p results logs

# Step 4: Distribute jobs across nodes
echo ""
echo "Step 2: Distributing $NUM_JOBS jobs across $NUM_NODES nodes..."

# Get job list
PYTHON_CMD='
import json
with open("param_sweep/job_configs.json", "r") as f:
    config = json.load(f)
for i, job in enumerate(config["jobs"]):
    print(f"{i} {job['job_name']} {job['param_file']}")
'

# Create job distribution
JOB_LIST=$(python3 -c "$PYTHON_CMD")

# Function to run job on remote node
run_remote_job() {
    local node=$1
    local job_id=$2
    local job_name=$3
    local param_file=$4
    
    echo "Starting job $job_id ($job_name) on $node..."
    
    # Create remote command
    REMOTE_CMD="
cd $(pwd) && 
mkdir -p results/$job_name logs &&
export PARAM_FILE=$param_file &&
export OUTPUT_DIR=results/$job_name &&
timeout 3600s ./build/pancreatic_tumor_new --param-file $param_file --output-dir results/$job_name > logs/${job_name}.out 2> logs/${job_name}.err &&
echo 'Job: $job_name' > results/$job_name/job_summary.txt &&
echo 'Status: SUCCESS' >> results/$job_name/job_summary.txt &&
echo 'Node: $node' >> results/$job_name/job_summary.txt &&
echo 'End Time: \$(date)' >> results/$job_name/job_summary.txt
"
    
    if [ "$node" = "localhost" ]; then
        # Run locally
        eval "$REMOTE_CMD" &
    else
        # Run on remote node
        ssh "$node" "$REMOTE_CMD" &
    fi
    
    echo "Job $job_id submitted to $node (PID: $!)"
}

# Distribute jobs
JOB_COUNTER=0
echo "$JOB_LIST" | while read -r job_id job_name param_file; do
    node_index=$((JOB_COUNTER % NUM_NODES))
    node=${NODES[$node_index]}
    
    run_remote_job "$node" "$job_id" "$job_name" "$param_file"
    
    JOB_COUNTER=$((JOB_COUNTER + 1))
    
    # Optional: Add delay between job submissions
    sleep 1
done

echo ""
echo "All jobs submitted!"
echo ""
echo "Monitor progress with:"
echo "  watch 'ls results/*/job_summary.txt 2>/dev/null | wc -l'"
echo ""
echo "Check logs:"
echo "  ls logs/"
echo ""
echo "When complete, analyze results:"
echo "  python3 analyze_results.py --results-dir results --output analysis_output"

echo ""
echo "============================================"