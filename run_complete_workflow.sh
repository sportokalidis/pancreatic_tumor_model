#!/bin/bash
# Complete HPC Parameter Sweep Workflow Script
# 
# This script runs the complete workflow:
# 1. Generate parameter files
# 2. Submit jobs to HPC scheduler 
# 3. Monitor completion
# 4. Analyze results
#
# Usage: ./run_complete_workflow.sh [sweep_type] [num_jobs] [scheduler]

set -e  # Exit on any error

# Default parameters
SWEEP_TYPE=${1:-"tumor_growth"}
NUM_JOBS=${2:-10}
SCHEDULER=${3:-"slurm"}  # slurm or pbs

echo "============================================"
echo "HPC Parameter Sweep Workflow"
echo "============================================"
echo "Sweep Type: $SWEEP_TYPE"
echo "Number of Jobs: $NUM_JOBS"
echo "Scheduler: $SCHEDULER"
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

# Step 2: Submit jobs
echo ""
echo "Step 2: Submitting jobs to $SCHEDULER..."

# Create logs directory
mkdir -p logs

if [ "$SCHEDULER" = "slurm" ]; then
    # Submit SLURM job array
    JOB_ID=$(sbatch --array=0-$((NUM_JOBS-1)) --parsable submit_param_sweep_slurm.sh)
    echo "✓ Submitted SLURM job array: $JOB_ID"
    
    # Monitor job progress
    echo ""
    echo "Step 3: Monitoring job progress..."
    echo "You can monitor jobs with: squeue -j $JOB_ID"
    echo "Or check all your jobs with: squeue -u \$USER"
    
    # Wait for jobs to complete (optional)
    read -p "Wait for jobs to complete before analysis? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Waiting for jobs to complete..."
        while squeue -j "$JOB_ID" 2>/dev/null | grep -q "$JOB_ID"; do
            echo "Jobs still running... (checking every 30 seconds)"
            sleep 30
        done
        echo "✓ All jobs completed"
    fi
    
elif [ "$SCHEDULER" = "pbs" ]; then
    # Modify PBS script for correct array range
    sed -i "s/#PBS -t 0-9/#PBS -t 0-$((NUM_JOBS-1))/g" submit_param_sweep_pbs.sh
    
    # Submit PBS job array
    JOB_ID=$(qsub submit_param_sweep_pbs.sh)
    echo "✓ Submitted PBS job array: $JOB_ID"
    
    echo ""
    echo "Step 3: Monitoring job progress..."
    echo "You can monitor jobs with: qstat -t $JOB_ID"
    echo "Or check all your jobs with: qstat -u \$USER"
    
    # Wait for jobs to complete (optional)
    read -p "Wait for jobs to complete before analysis? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Waiting for jobs to complete..."
        while qstat "$JOB_ID" 2>/dev/null | grep -q "$JOB_ID"; do
            echo "Jobs still running... (checking every 30 seconds)"
            sleep 30
        done
        echo "✓ All jobs completed"
    fi
    
else
    echo "Error: Unknown scheduler '$SCHEDULER'. Use 'slurm' or 'pbs'"
    exit 1
fi

# Step 4: Analyze results (optional)
echo ""
echo "Step 4: Results analysis (optional)..."
read -p "Run results analysis now? (y/n): " -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]; then
    if [ -d "results" ]; then
        echo "Running analysis..."
        python3 analyze_results.py \
            --results-dir results \
            --output analysis_output \
            --sweep-type "$SWEEP_TYPE"
        
        echo "✓ Analysis complete!"
        echo ""
        echo "Results available in:"
        echo "  - analysis_output/simulation_summary.csv"
        echo "  - analysis_output/analysis_report.md"
        echo "  - analysis_output/plots/"
    else
        echo "Warning: Results directory not found. Jobs may still be running."
        echo "Run analysis later with:"
        echo "  python3 analyze_results.py --results-dir results --output analysis_output"
    fi
fi

echo ""
echo "============================================"
echo "Workflow Summary:"
echo "- Parameter files: param_sweep/"
echo "- Job submission: $SCHEDULER job array"
echo "- Results: results/"
echo "- Analysis: analysis_output/"
echo "============================================"
echo ""
echo "Useful commands:"
if [ "$SCHEDULER" = "slurm" ]; then
    echo "  Monitor jobs: squeue -u \$USER"
    echo "  Cancel jobs: scancel $JOB_ID"
else
    echo "  Monitor jobs: qstat -u \$USER" 
    echo "  Cancel jobs: qdel $JOB_ID"
fi
echo "  Check logs: ls logs/"
echo "  Analyze results: python3 analyze_results.py --results-dir results"
echo ""
echo "See HPC_WORKFLOW.md for detailed documentation"