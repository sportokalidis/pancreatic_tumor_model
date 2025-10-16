#!/bin/bash
# Local Parameter Sweep Runner
# Simulates HPC job array execution on local machine
#
# Usage: ./run_local_sweep.sh [sweep_type] [num_jobs]

set -e

SWEEP_TYPE=${1:-"tumor_growth"}
NUM_JOBS=${2:-5}

echo "============================================"
echo "Local Parameter Sweep Runner"
echo "============================================"
echo "Sweep Type: $SWEEP_TYPE"
echo "Number of Jobs: $NUM_JOBS" 
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

# Step 2: Create results and logs directories
mkdir -p results logs

# Step 3: Run all simulations locally
echo ""
echo "Step 2: Running simulations locally..."
echo "This may take a while for large parameter sweeps..."

# Read job configurations
PYTHON_CMD='
import json
import sys

with open("param_sweep/job_configs.json", "r") as f:
    config = json.load(f)

for i, job in enumerate(config["jobs"]):
    print(f"{i} {job['job_name']} {job['param_file']}")
'

# Run each simulation
SUCCESS_COUNT=0
FAILED_COUNT=0

python3 -c "$PYTHON_CMD" | while read -r job_id job_name param_file; do
    echo "Running simulation $job_id: $job_name..."
    
    output_dir="results/$job_name"
    mkdir -p "$output_dir"
    
    # Copy parameter file to output directory
    cp "$param_file" "$output_dir/params.txt"
    
    # Run simulation with timeout
    start_time=$(date +%s)
    if timeout 300s ./build/pancreatic_tumor_new \
        --param-file "$param_file" \
        --output-dir "$output_dir" \
        > "logs/${job_name}.out" 2> "logs/${job_name}.err"; then
        
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        
        # Create job summary
        cat > "$output_dir/job_summary.txt" << EOF
Job: $job_name
Description: Local parameter sweep
Parameter File: $param_file
Exit Code: 0
Status: SUCCESS
Duration: ${duration}s
Start Time: $(date -d @$start_time)
End Time: $(date -d @$end_time)
EOF
        
        echo "  ✓ Completed successfully in ${duration}s"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        
    else
        exit_code=$?
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        
        # Create failure summary
        cat > "$output_dir/job_summary.txt" << EOF
Job: $job_name
Description: Local parameter sweep
Parameter File: $param_file
Exit Code: $exit_code
Status: FAILED
Duration: ${duration}s
Start Time: $(date -d @$start_time)
End Time: $(date -d @$end_time)
EOF
        
        echo "  ✗ Failed with exit code $exit_code"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
done

# Wait for background processes
wait

# Count actual results
SUCCESS_COUNT=$(find results -name "job_summary.txt" -exec grep -l "Status: SUCCESS" {} \; | wc -l)
FAILED_COUNT=$(find results -name "job_summary.txt" -exec grep -l "Status: FAILED" {} \; | wc -l)

echo ""
echo "============================================"
echo "Step 3: Summary"
echo "============================================"
echo "Successful simulations: $SUCCESS_COUNT"
echo "Failed simulations: $FAILED_COUNT"
echo "Total simulations: $((SUCCESS_COUNT + FAILED_COUNT))"

if [ $SUCCESS_COUNT -gt 0 ]; then
    echo ""
    echo "Step 4: Running analysis..."
    python3 analyze_results.py \
        --results-dir results \
        --output analysis_output \
        --sweep-type "$SWEEP_TYPE"
    
    echo ""
    echo "✓ Analysis complete!"
    echo ""
    echo "Results available in:"
    echo "  - results/               # Individual simulation outputs"
    echo "  - analysis_output/       # Aggregated analysis"
    echo "  - logs/                  # Simulation logs"
else
    echo ""
    echo "⚠️  No successful simulations to analyze"
    echo "Check logs/ directory for error details"
fi

echo ""
echo "============================================"