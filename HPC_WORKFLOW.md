# HPC Parameter Sweep Workflow

This guide explains how to run parameter sweeps for the pancreatic tumor model on HPC systems.

## Overview

The workflow consists of three main steps:
1. **Parameter Generation**: Create multiple parameter files for different scenarios
2. **Job Submission**: Submit array jobs to HPC scheduler (SLURM or PBS)
3. **Results Analysis**: Aggregate and analyze results from all simulations

## Files Overview

- `generate_param_sweep.py`: Generate parameter files and job configurations
- `submit_param_sweep_slurm.sh`: SLURM job array script
- `submit_param_sweep_pbs.sh`: PBS job array script  
- `analyze_results.py`: Results analysis and visualization
- `run_complete_workflow.sh`: End-to-end workflow script

## Quick Start

### 1. Generate Parameter Files

```bash
# Generate tumor growth parameter sweep (10 jobs)
python3 generate_param_sweep.py --sweep-type tumor_growth --num-jobs 10

# Generate immune response parameter sweep (15 jobs)
python3 generate_param_sweep.py --sweep-type immune_response --num-jobs 15

# Generate comprehensive sweep with all parameters (20 jobs)
python3 generate_param_sweep.py --sweep-type comprehensive --num-jobs 20
```

### 2. Submit Jobs to HPC

#### For SLURM systems:
```bash
# Make script executable
chmod +x submit_param_sweep_slurm.sh

# Submit job array (adjust array range to match number of parameter files)
sbatch --array=0-9 submit_param_sweep_slurm.sh
```

#### For PBS systems:
```bash
# Make script executable
chmod +x submit_param_sweep_pbs.sh

# Submit job array
qsub submit_param_sweep_pbs.sh
```

### 3. Monitor Jobs

#### SLURM:
```bash
# Check job status
squeue -u $USER

# Check specific job
squeue -j <job_id>

# Cancel jobs if needed
scancel <job_id>
```

#### PBS:
```bash
# Check job status
qstat -u $USER

# Check specific job
qstat <job_id>

# Cancel jobs if needed
qdel <job_id>
```

### 4. Analyze Results

```bash
# Install required Python packages (if not already installed)
pip install pandas matplotlib numpy

# Run analysis (seaborn is optional for advanced plots)
python3 analyze_results.py --results-dir results --output analysis_output

# Generate specific analysis
python3 analyze_results.py --sweep-type tumor_growth --format pdf
```

## Detailed Usage

### Parameter Generation Options

The `generate_param_sweep.py` script supports several sweep types:

#### Tumor Growth Sweep
Varies PDC replication probability and max simulation steps:
```bash
python3 generate_param_sweep.py --sweep-type tumor_growth --num-jobs 10 --output-dir param_sweep
```

#### Immune Response Sweep  
Varies Treg suppression and immune cell activation:
```bash
python3 generate_param_sweep.py --sweep-type immune_response --num-jobs 15
```

#### Initial Conditions Sweep
Varies initial cell populations:
```bash
python3 generate_param_sweep.py --sweep-type initial_conditions --num-jobs 12
```

#### Treg Suppression Sweep
Focuses on regulatory T cell effects:
```bash
python3 generate_param_sweep.py --sweep-type treg_suppression --num-jobs 8
```

### Custom Parameter Ranges

You can modify the parameter ranges in `generate_param_sweep.py`:

```python
# Example: Custom tumor growth parameters
tumor_growth_params = {
    'pdc_replication_prob': np.linspace(0.05, 0.20, num_jobs),
    'max_steps': np.linspace(500, 2000, num_jobs, dtype=int),
    # Add more parameters as needed
}
```

### HPC System Configuration

#### SLURM Configuration
Edit `submit_param_sweep_slurm.sh` to match your HPC system:

```bash
#SBATCH --partition=your_partition
#SBATCH --account=your_account
#SBATCH --cpus-per-task=4        # Increase for multi-core jobs
#SBATCH --mem=8GB                # Adjust memory requirements
#SBATCH --time=04:00:00          # Adjust time limit

# Load your system's modules
module load gcc/9.3.0
module load cmake/3.18.0
```

#### PBS Configuration
Edit `submit_param_sweep_pbs.sh` similarly:

```bash
#PBS -l nodes=1:ppn=4            # Adjust CPU count
#PBS -l mem=8gb                  # Adjust memory
#PBS -l walltime=04:00:00        # Adjust time limit
```

### Results Structure

After jobs complete, results are organized as:

```
results/
├── tumor_growth_job_0/
│   ├── params.txt               # Parameter file used
│   ├── populations.csv          # Population time series
│   ├── job_summary.txt          # Job metadata
│   └── *-Cells.csv             # Individual cell type data
├── tumor_growth_job_1/
│   └── ...
└── ...
```

### Analysis Outputs

The analysis script generates:

```
analysis_output/
├── simulation_summary.csv       # Summary table of all runs
├── analysis_report.md          # Comprehensive report
└── plots/
    ├── tumor_vs_pdc_replication.png
    ├── immune_vs_treg_suppression.png
    ├── population_dynamics.png
    └── correlation_matrix.png
```

## Troubleshooting

### Common Issues

1. **Job Arrays Not Working**
   - Check that array range matches number of parameter files
   - Verify `job_configs.json` exists in `param_sweep/` directory

2. **Simulations Failing**
   - Check log files in `logs/` directory
   - Verify parameter file format is correct
   - Ensure simulation binary exists and is executable

3. **Analysis Script Errors**
   - Install required Python packages: `pip install pandas matplotlib numpy`
   - For advanced plots: `pip install seaborn`
   - Check that results directory contains expected files

4. **Memory/Time Limits**
   - Increase memory allocation in job scripts
   - Extend walltime for longer simulations
   - Consider reducing simulation complexity

### Debugging Tips

1. **Test Single Job First**
   ```bash
   # Test with one parameter file
   export PARAM_FILE="param_sweep/tumor_growth_job_0.txt"
   export OUTPUT_DIR="test_output"
   ./build/pancreatic_tumor_new --param-file $PARAM_FILE --output-dir $OUTPUT_DIR
   ```

2. **Check Job Logs**
   ```bash
   # SLURM
   less logs/job_<jobid>_<arrayid>.out
   
   # PBS  
   less logs/job_<jobid>.out
   ```

3. **Verify Parameter Files**
   ```bash
   # Check generated parameters
   head -20 param_sweep/tumor_growth_job_0.txt
   ```

## Advanced Usage

### Custom Sweep Types

Add new sweep types to `generate_param_sweep.py`:

```python
elif sweep_type == 'custom_sweep':
    # Define your custom parameter ranges
    custom_params = {
        'your_param1': np.linspace(min_val, max_val, num_jobs),
        'your_param2': np.logspace(log_min, log_max, num_jobs),
    }
    
    for i in range(num_jobs):
        params = base_params.copy()
        params.update({
            'your_param1': custom_params['your_param1'][i],
            'your_param2': custom_params['your_param2'][i],
        })
        # ... rest of parameter generation
```

### Parallel Analysis

For large result sets, consider parallel analysis:

```bash
# Split analysis by sweep type
python3 analyze_results.py --sweep-type tumor_growth &
python3 analyze_results.py --sweep-type immune_response &
wait
```

### Integration with Workflow Managers

For complex workflows, consider using workflow managers like Snakemake or Nextflow to coordinate parameter generation, job submission, and analysis steps.

## Performance Considerations

- **Memory Usage**: Monitor memory usage and adjust job requests accordingly
- **I/O Performance**: For large parameter sweeps, consider staging data on fast storage
- **CPU Utilization**: Balance single-core jobs vs. multi-core simulations based on your model
- **Queue Limits**: Check your HPC system's job array limits and split large sweeps if necessary

## Contact and Support

For issues specific to:
- BioDynaMo framework: Check BioDynaMo documentation
- HPC system: Consult your local HPC support team
- This workflow: Review logs and parameter files for debugging information