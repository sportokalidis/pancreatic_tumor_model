# HPC Automation Without SLURM: Complete Guide

This guide shows how to automate parameter sweeps on different HPC systems that don't use SLURM.

## 🔍 **Step 1: Identify Your HPC Scheduler**

Run these commands on your HPC system to identify the scheduler:

```bash
# Check for different schedulers
which qsub     # PBS/Torque or SGE
which bsub     # LSF
which condor_submit  # HTCondor

# Check system info
qstat --version 2>/dev/null || echo "No qstat"
bjobs --version 2>/dev/null || echo "No bjobs"  
condor_version 2>/dev/null || echo "No condor"
```

## 🎛️ **Step 2: Choose Your Automation Method**

### **Method A: PBS/Torque Systems**

**Commands:**
```bash
# 1. Generate parameters
python3 generate_param_sweep.py --sweep-type tumor_growth --num-jobs 50

# 2. Update PBS script array range
sed -i 's/#PBS -t 0-9/#PBS -t 0-49/' submit_param_sweep_pbs.sh

# 3. Customize for your system (edit the file)
nano submit_param_sweep_pbs.sh
# Add lines like:
# #PBS -A your_project_account
# #PBS -q your_queue_name

# 4. Submit jobs
qsub submit_param_sweep_pbs.sh

# 5. Monitor
qstat -u $USER
watch 'qstat -u $USER'

# 6. Cancel if needed
qdel JOB_ID
```

**PBS Script Customization:**
```bash
#PBS -A project_name      # Your project/account
#PBS -q gpu_queue         # Queue name
#PBS -l nodes=1:ppn=4     # 4 CPU cores
#PBS -l mem=8gb           # 8GB memory
#PBS -l walltime=04:00:00 # 4 hour limit
```

---

### **Method B: LSF Systems**

**Commands:**
```bash
# 1. Generate parameters
python3 generate_param_sweep.py --sweep-type tumor_growth --num-jobs 50

# 2. Update LSF script array range
sed -i 's/\[0-9\]/[0-49]/' submit_param_sweep_lsf.sh

# 3. Customize for your system
nano submit_param_sweep_lsf.sh

# 4. Submit jobs
bsub < submit_param_sweep_lsf.sh

# 5. Monitor
bjobs -u $USER
watch 'bjobs -u $USER'

# 6. Cancel if needed
bkill JOB_ID
```

**LSF Script Customization:**
```bash
#BSUB -P project_name         # Project name
#BSUB -q gpu_queue            # Queue name  
#BSUB -n 4                    # 4 CPU cores
#BSUB -R "rusage[mem=8000]"   # 8GB memory per core
#BSUB -W 4:00                 # 4 hour limit
```

---

### **Method C: SGE Systems**

**Commands:**
```bash
# 1. Generate parameters
python3 generate_param_sweep.py --sweep-type tumor_growth --num-jobs 50

# 2. Update SGE script (SGE uses 1-based indexing!)
sed -i 's/#$ -t 1-10/#$ -t 1-50/' submit_param_sweep_sge.sh

# 3. Customize for your system
nano submit_param_sweep_sge.sh

# 4. Submit jobs
qsub submit_param_sweep_sge.sh

# 5. Monitor
qstat -u $USER

# 6. Cancel if needed
qdel JOB_ID
```

**SGE Script Customization:**
```bash
#$ -P project_name           # Project name
#$ -q gpu.q                  # Queue name
#$ -pe smp 4                 # 4 CPU cores
#$ -l h_vmem=8G              # 8GB memory
#$ -l h_rt=04:00:00          # 4 hour limit
```

---

### **Method D: Manual/SSH Distribution**

For systems without job schedulers:

```bash
# 1. Create nodes file
echo "node001" > nodes.txt
echo "node002" >> nodes.txt
echo "node003" >> nodes.txt

# 2. Run manual distribution
./run_manual_hpc.sh tumor_growth 50 nodes.txt

# 3. Monitor progress
watch 'ls results/*/job_summary.txt 2>/dev/null | wc -l'
```

---

## 📊 **Step 3: Monitor and Analyze**

### **Monitor Job Progress:**

```bash
# Check running jobs (adjust command for your scheduler)
qstat -u $USER        # PBS/SGE
bjobs -u $USER         # LSF

# Count completed jobs
ls results/*/job_summary.txt 2>/dev/null | wc -l

# Check for failures
grep -l "FAILED" results/*/job_summary.txt 2>/dev/null

# View logs
tail -f logs/job_*.out
```

### **Results Analysis:**

```bash
# When jobs complete, analyze results
python3 analyze_results.py --results-dir results --output analysis_output

# View summary
cat analysis_output/analysis_report.md
cat analysis_output/simulation_summary.csv
```

## 🛠️ **Troubleshooting**

### **Common Issues:**

1. **Wrong Array Range:**
   ```bash
   # PBS/LSF: 0-based indexing (0-49 for 50 jobs)
   # SGE: 1-based indexing (1-50 for 50 jobs)
   ```

2. **Module Loading:**
   ```bash
   # Add to your submission script:
   module load gcc/9.3.0
   module load cmake/3.18.0
   module load python/3.8
   ```

3. **Job Limits:**
   ```bash
   # Check system limits
   qstat -Q          # PBS
   bqueues           # LSF
   qconf -sq all.q   # SGE
   ```

4. **Memory Issues:**
   ```bash
   # Increase memory allocation in job script
   #PBS -l mem=16gb        # PBS
   #BSUB -R "rusage[mem=16000]"  # LSF
   #$ -l h_vmem=16G        # SGE
   ```

## 📋 **Quick Reference**

| Scheduler | Submit Command | Monitor Command | Cancel Command |
|-----------|----------------|-----------------|----------------|
| **PBS/Torque** | `qsub script.sh` | `qstat -u $USER` | `qdel JOB_ID` |
| **LSF** | `bsub < script.sh` | `bjobs -u $USER` | `bkill JOB_ID` |
| **SGE** | `qsub script.sh` | `qstat -u $USER` | `qdel JOB_ID` |
| **Manual** | `./run_manual_hpc.sh` | `watch 'ls results/*/'` | Kill processes |

## 🚀 **Example Complete Workflow:**

```bash
# For PBS system with 100 parameter combinations:

# Step 1: Generate parameters
python3 generate_param_sweep.py --sweep-type immune_response --num-jobs 100

# Step 2: Update script
sed -i 's/#PBS -t 0-9/#PBS -t 0-99/' submit_param_sweep_pbs.sh

# Step 3: Submit
qsub submit_param_sweep_pbs.sh

# Step 4: Monitor (job will have ID like 12345[])
qstat -t 12345

# Step 5: Analyze when complete
python3 analyze_results.py --results-dir results --output final_analysis
```

This approach works on virtually any HPC system! 🎯