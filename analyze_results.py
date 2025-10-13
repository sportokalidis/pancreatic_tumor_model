#!/usr/bin/env python3
"""
Results Analysis Script for Pancreatic Tumor Model Parameter Sweep

This script analyzes the results from multiple parameter sweep simulations,
aggregates data, and generates comparison plots.

Usage:
    python3 analyze_results.py --results-dir results --output plots
    python3 analyze_results.py --help
"""

import argparse
import json
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Optional seaborn import
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze parameter sweep results for pancreatic tumor model"
    )
    parser.add_argument(
        "--results-dir", 
        type=str, 
        default="results",
        help="Directory containing simulation results (default: results)"
    )
    parser.add_argument(
        "--output", 
        type=str, 
        default="analysis_output",
        help="Output directory for analysis results (default: analysis_output)"
    )
    parser.add_argument(
        "--sweep-type",
        type=str,
        choices=['tumor_growth', 'immune_response', 'initial_conditions', 'treg_suppression'],
        help="Specific sweep type to analyze (if not specified, analyzes all)"
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=['png', 'pdf', 'svg'],
        default='png',
        help="Output format for plots (default: png)"
    )
    
    return parser.parse_args()

def load_job_configs(param_sweep_dir):
    """Load job configurations from parameter sweep."""
    config_file = param_sweep_dir / "job_configs.json"
    if not config_file.exists():
        print(f"Warning: {config_file} not found")
        return None
    
    with open(config_file, 'r') as f:
        return json.load(f)

def parse_parameter_file(param_file):
    """Parse parameter file and return dictionary."""
    params = {}
    try:
        with open(param_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        
                        # Try to convert to appropriate type
                        try:
                            # Try integer first
                            params[key] = int(value)
                        except ValueError:
                            try:
                                # Try float
                                params[key] = float(value)
                            except ValueError:
                                # Keep as string
                                params[key] = value
    except Exception as e:
        print(f"Error parsing {param_file}: {e}")
    
    return params

def load_simulation_data(results_dir):
    """Load all simulation results and metadata."""
    results_dir = Path(results_dir)
    
    if not results_dir.exists():
        print(f"Error: Results directory {results_dir} does not exist")
        return None
    
    simulations = []
    
    # Iterate through result directories
    for job_dir in results_dir.iterdir():
        if not job_dir.is_dir():
            continue
            
        print(f"Processing {job_dir.name}...")
        
        simulation_data = {
            'job_name': job_dir.name,
            'job_dir': job_dir,
            'parameters': {},
            'success': False,
            'populations': None,
            'final_populations': {},
            'summary': {}
        }
        
        # Load job summary
        summary_file = job_dir / "job_summary.txt"
        if summary_file.exists():
            with open(summary_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if ':' in line:
                        key, value = line.split(':', 1)
                        simulation_data['summary'][key.strip()] = value.strip()
            
            simulation_data['success'] = 'SUCCESS' in simulation_data['summary'].get('Status', 'FAILED') or \
                                       simulation_data['summary'].get('Exit Code', '1') == '0'
        
        # Load parameters
        param_file = job_dir / "params.txt"
        if param_file.exists():
            simulation_data['parameters'] = parse_parameter_file(param_file)
        
        # Load population data
        pop_file = job_dir / "populations.csv"
        if pop_file.exists():
            try:
                df = pd.read_csv(pop_file)
                simulation_data['populations'] = df
                
                # Get final population counts
                if not df.empty:
                    final_row = df.iloc[-1]
                    for col in df.columns:
                        if col != 'iteration':
                            simulation_data['final_populations'][col] = final_row[col]
                
            except Exception as e:
                print(f"Error loading {pop_file}: {e}")
        
        simulations.append(simulation_data)
    
    print(f"Loaded {len(simulations)} simulations")
    return simulations

def create_summary_table(simulations, output_dir):
    """Create a summary table of all simulations."""
    output_dir = Path(output_dir)
    
    summary_data = []
    
    for sim in simulations:
        row = {
            'Job Name': sim['job_name'],
            'Success': sim['success'],
            'Description': sim['summary'].get('Description', 'N/A')
        }
        
        # Add key parameters
        for key in ['max_steps', 'pdc_replication_prob', 'treg_suppression', 
                   'initial_tumor_cells', 'initial_immune_cells']:
            if key in sim['parameters']:
                row[key] = sim['parameters'][key]
        
        # Add final populations
        for key in ['total_tumor', 'total_immune', 'N-Cells', 'P-Cells', 'H-Cells']:
            if key in sim['final_populations']:
                row[f'final_{key}'] = sim['final_populations'][key]
        
        summary_data.append(row)
    
    # Create DataFrame and save
    df = pd.DataFrame(summary_data)
    
    # Save as CSV
    csv_file = output_dir / "simulation_summary.csv"
    df.to_csv(csv_file, index=False)
    print(f"Summary table saved to {csv_file}")
    
    return df

def plot_parameter_effects(simulations, output_dir, file_format='png'):
    """Create plots showing parameter effects on outcomes."""
    output_dir = Path(output_dir)
    
    # Create plots directory
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    successful_sims = [sim for sim in simulations if sim['success']]
    
    if not successful_sims:
        print("No successful simulations found for plotting")
        return
    
    # Extract data for plotting
    plot_data = []
    for sim in successful_sims:
        row = {'job_name': sim['job_name']}
        row.update(sim['parameters'])
        row.update(sim['final_populations'])
        plot_data.append(row)
    
    if not plot_data:
        print("No data available for plotting")
        return
    
    df = pd.DataFrame(plot_data)
    
    # Set style
    plt.style.use('default')
    if HAS_SEABORN:
        sns.set_palette("husl")
    
    # 1. Tumor growth vs PDC replication probability
    if 'pdc_replication_prob' in df.columns and 'total_tumor' in df.columns:
        plt.figure(figsize=(10, 6))
        plt.scatter(df['pdc_replication_prob'], df['total_tumor'], alpha=0.7, s=60)
        plt.xlabel('PDC Replication Probability')
        plt.ylabel('Final Tumor Cell Count')
        plt.title('Tumor Growth vs PDC Replication Probability')
        plt.grid(True, alpha=0.3)
        
        # Add trend line
        if len(df) > 1:
            z = np.polyfit(df['pdc_replication_prob'], df['total_tumor'], 1)
            p = np.poly1d(z)
            plt.plot(df['pdc_replication_prob'], p(df['pdc_replication_prob']), "r--", alpha=0.8)
        
        plt.tight_layout()
        plt.savefig(plots_dir / f"tumor_vs_pdc_replication.{file_format}", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 2. Immune response effectiveness
    if 'treg_suppression' in df.columns and 'total_immune' in df.columns:
        plt.figure(figsize=(10, 6))
        plt.scatter(df['treg_suppression'], df['total_immune'], alpha=0.7, s=60, c='green')
        plt.xlabel('Treg Suppression Rate')
        plt.ylabel('Final Immune Cell Count')
        plt.title('Immune Response vs Treg Suppression')
        plt.grid(True, alpha=0.3)
        
        if len(df) > 1:
            z = np.polyfit(df['treg_suppression'], df['total_immune'], 1)
            p = np.poly1d(z)
            plt.plot(df['treg_suppression'], p(df['treg_suppression']), "r--", alpha=0.8)
        
        plt.tight_layout()
        plt.savefig(plots_dir / f"immune_vs_treg_suppression.{file_format}", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Population dynamics over time (if we have time series data)
    time_series_sims = [sim for sim in successful_sims if sim['populations'] is not None]
    
    if time_series_sims:
        plt.figure(figsize=(12, 8))
        
        for i, sim in enumerate(time_series_sims[:5]):  # Plot first 5 simulations
            df_pop = sim['populations']
            if 'iteration' in df_pop.columns:
                for col in ['N-Cells', 'P-Cells', 'H-Cells']:
                    if col in df_pop.columns:
                        plt.plot(df_pop['iteration'], df_pop[col], 
                               label=f"{col} - {sim['job_name'][:15]}", alpha=0.7)
        
        plt.xlabel('Iteration')
        plt.ylabel('Cell Count')
        plt.title('Population Dynamics Over Time')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(plots_dir / f"population_dynamics.{file_format}", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 4. Parameter correlation heatmap
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    param_cols = [col for col in numeric_cols if col in ['pdc_replication_prob', 'treg_suppression', 
                                                        'initial_tumor_cells', 'initial_immune_cells', 'max_steps']]
    outcome_cols = [col for col in numeric_cols if col.startswith('total_') or col in ['N-Cells', 'P-Cells', 'H-Cells']]
    
    if len(param_cols) > 1 and len(outcome_cols) > 0:
        corr_data = df[param_cols + outcome_cols].corr()
        
        plt.figure(figsize=(10, 8))
        if HAS_SEABORN:
            sns.heatmap(corr_data, annot=True, cmap='coolwarm', center=0, 
                       square=True, fmt='.2f', cbar_kws={"shrink": .8})
        else:
            plt.imshow(corr_data, cmap='coolwarm', aspect='auto')
            plt.colorbar()
            # Add text annotations
            for i in range(len(corr_data)):
                for j in range(len(corr_data.columns)):
                    plt.text(j, i, f'{corr_data.iloc[i, j]:.2f}', 
                            ha='center', va='center')
            plt.xticks(range(len(corr_data.columns)), corr_data.columns, rotation=45)
            plt.yticks(range(len(corr_data)), corr_data.index)
        plt.title('Parameter-Outcome Correlation Matrix')
        plt.tight_layout()
        plt.savefig(plots_dir / f"correlation_matrix.{file_format}", dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Plots saved to {plots_dir}")

def generate_report(simulations, summary_df, output_dir):
    """Generate a comprehensive analysis report."""
    output_dir = Path(output_dir)
    report_file = output_dir / "analysis_report.md"
    
    with open(report_file, 'w') as f:
        f.write("# Pancreatic Tumor Model Parameter Sweep Analysis\n\n")
        f.write(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Summary statistics
        total_sims = len(simulations)
        successful_sims = len([s for s in simulations if s['success']])
        
        f.write("## Summary Statistics\n\n")
        f.write(f"- Total simulations: {total_sims}\n")
        f.write(f"- Successful simulations: {successful_sims}\n")
        f.write(f"- Success rate: {successful_sims/total_sims*100:.1f}%\n\n")
        
        # Parameter ranges
        f.write("## Parameter Ranges\n\n")
        if successful_sims > 0:
            param_summary = {}
            for sim in simulations:
                if sim['success']:
                    for key, value in sim['parameters'].items():
                        if isinstance(value, (int, float)):
                            if key not in param_summary:
                                param_summary[key] = []
                            param_summary[key].append(value)
            
            for key, values in param_summary.items():
                if values:
                    f.write(f"- {key}: {min(values):.3f} - {max(values):.3f} (mean: {np.mean(values):.3f})\n")
        
        f.write("\n## Key Findings\n\n")
        
        # Analyze outcomes
        if successful_sims > 0:
            outcomes = {}
            for sim in simulations:
                if sim['success']:
                    for key, value in sim['final_populations'].items():
                        if isinstance(value, (int, float)):
                            if key not in outcomes:
                                outcomes[key] = []
                            outcomes[key].append(value)
            
            for key, values in outcomes.items():
                if values:
                    f.write(f"- {key}: mean = {np.mean(values):.1f}, std = {np.std(values):.1f}, range = {min(values):.1f}-{max(values):.1f}\n")
        
        f.write("\n## Files Generated\n\n")
        f.write("- `simulation_summary.csv`: Detailed results table\n")
        f.write("- `plots/`: Directory containing analysis plots\n")
        f.write("- `analysis_report.md`: This report\n")
    
    print(f"Analysis report saved to {report_file}")

def main():
    """Main analysis function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Analyzing results from: {args.results_dir}")
    print(f"Output directory: {output_dir}")
    
    # Load simulation data
    simulations = load_simulation_data(args.results_dir)
    if not simulations:
        print("No simulation data found")
        return 1
    
    # Filter by sweep type if specified
    if args.sweep_type:
        simulations = [s for s in simulations if args.sweep_type in s['job_name']]
        print(f"Filtered to {len(simulations)} {args.sweep_type} simulations")
    
    # Create summary table
    summary_df = create_summary_table(simulations, output_dir)
    
    # Create plots
    plot_parameter_effects(simulations, output_dir, args.format)
    
    # Generate report
    generate_report(simulations, summary_df, output_dir)
    
    print(f"\nAnalysis complete! Results saved to: {output_dir}")
    return 0

if __name__ == "__main__":
    sys.exit(main())