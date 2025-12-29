import subprocess
import json
import os
import sys
import re
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing
import csv

# =============================================================================
# CONFIGURATION
# =============================================================================
EXECUTABLE  = "../build/main"
CONFIG_PATH = "../data/config.json"
OUTPUT_DIR  = "../output"
RESULTS_DIR = "../results"

# Base Grid Size for 1 Processor
# (e.g., 64^3 = 262k cells per core. Keep this heavy enough to mask overhead)
BASE_GRID_SIZE = 64
NUM_STEPS = 50

MAX_CORES = multiprocessing.cpu_count()
PROCS_TO_TEST = [1, 2, 4, 8]
if 16 <= MAX_CORES: PROCS_TO_TEST.append(16)
# =============================================================================

def update_config_weak(n_procs):
    """
    Updates config.json for Weak Scaling:
    1. Scale Grid dimensions so TotalCells approx = BaseCells * N_Procs
    2. Disable I/O.
    """
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    with open(CONFIG_PATH, 'r') as f:
        config = json.load(f)

    # --- WEAK SCALING LOGIC ---
    # We want Volume_N = Volume_1 * N_procs
    # L_new^3 = L_base^3 * N_procs
    # L_new = L_base * (N_procs)^(1/3)
    scale_factor = n_procs ** (1.0/3.0)
    new_N = int(round(BASE_GRID_SIZE * scale_factor))

    # Ensure dimensions are even to be friendly to decomposition
    if new_N % 2 != 0: new_N += 1

    config['mesh']['nx'] = new_N
    config['mesh']['ny'] = new_N
    config['mesh']['nz'] = new_N

    # Standard settings
    config['time']['dt'] = 0.001 # Small dt for stability on fine grids
    config['time']['t_end'] = config['time']['dt'] * NUM_STEPS
    config['output']['output_frequency'] = NUM_STEPS + 100 # Disable I/O
    config['output']['base_filename'] = "weak_test"

    with open(CONFIG_PATH, 'w') as f:
        json.dump(config, f, indent=4)

    total_cells = new_N**3
    print(f"[SETUP NP={n_procs}] Grid: {new_N}x{new_N}x{new_N} | Total Cells: {total_cells} ({total_cells/n_procs:.0f}/core)")
    return new_N, total_cells

def parse_output(stdout_str):
    metrics = {}

    # Regex patterns
    re_total = r"Total CPU Time \(solve loop\):\s+([\d\.]+)"
    re_norm  = r"CPU Time / \(Steps \* Cells\):\s+([\d\.eE\+\-]+)"

    match_total = re.search(re_total, stdout_str)
    match_norm  = re.search(re_norm, stdout_str)

    if match_total: metrics['total_time'] = float(match_total.group(1))
    if match_norm:  metrics['norm_metric']= float(match_norm.group(1))

    return metrics

def run_test(n_procs):
    # Update config and get grid info
    grid_dim, total_cells = update_config_weak(n_procs)

    cmd = [EXECUTABLE]
    if n_procs > 1:
        cmd = ["mpirun", "-np", str(n_procs)] + cmd

    print(f"  Running... ", end="", flush=True)
    try:
        # Run and capture output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        metrics = parse_output(result.stdout)

        # Inject grid info into metrics for logging later
        metrics['grid_dim'] = grid_dim
        metrics['total_cells'] = total_cells

        print(f"Done. (Time: {metrics.get('total_time', 'N/A')}s)")
        return metrics
    except subprocess.CalledProcessError as e:
        print(f"FAILED. {e}")
        return None

def plot_results(results):
    procs = [r['np'] for r in results]
    times = [r['metrics']['total_time'] for r in results]

    # Weak Scaling Efficiency: E = T_1 / T_N
    # (Since workload scales with P, T_N should ideally equal T_1)
    t_base = times[0]
    efficiency = [(t_base / t) * 100.0 for t in times]

    plt.figure(figsize=(10, 5))

    # --- Plot 1: Execution Time (Ideal is Flat) ---
    plt.subplot(1, 2, 1)
    plt.plot(procs, times, 'o-', linewidth=2, label="Measured Time")
    plt.axhline(t_base, color='k', linestyle='--', alpha=0.5, label="Ideal (Constant)")
    plt.title(f"Weak Scaling (Work/Core ~ {BASE_GRID_SIZE}^3)")
    plt.xlabel("Number of Processors")
    plt.ylabel("Execution Time (s)")
    plt.ylim(0, max(times)*1.2)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xticks(procs)

    # --- Plot 2: Efficiency ---
    plt.subplot(1, 2, 2)
    plt.plot(procs, efficiency, 'rs-', linewidth=2)
    plt.axhline(100, color='k', linestyle='--', alpha=0.5)
    plt.title("Weak Scaling Efficiency")
    plt.xlabel("Number of Processors")
    plt.ylabel("Efficiency (%)")
    plt.ylim(0, 110)
    plt.grid(True, alpha=0.3)
    plt.xticks(procs)

    plt.tight_layout()

    # Save in both formats
    base_plot_name = os.path.join(RESULTS_DIR, "weak_scalability_results")
    plt.savefig(f"{base_plot_name}.png")
    plt.savefig(f"{base_plot_name}.eps", format='eps')

    print(f"\n[OUTPUT] Plots saved to:")
    print(f"  - {base_plot_name}.png")
    print(f"  - {base_plot_name}.eps")

def main():
    print("============================================================")
    print(" WEAK SCALABILITY TEST SUITE")
    print("============================================================")

    print(f"[WARNING] Available CPU Cores: {MAX_CORES}")

    results = []

    for np in PROCS_TO_TEST:
        metrics = run_test(np)
        if metrics:
            results.append({'np': np, 'metrics': metrics})

    if not results: return

    # --- Print Table & Save CSV ---
    t_base = results[0]['metrics']['total_time']
    csv_filename = os.path.join(RESULTS_DIR, "weak_scalability_results.csv")

    print("\n" + "="*85)
    print(f"{'NP':<5} | {'Grid':<15} | {'Total Cells':<12} | {'Time (s)':<10} | {'Eff (%)':<10} | {'Metric (s/step/cell)'}")
    print("-" * 85)

    with open(csv_filename, 'w', newline='') as csvfile:
        fieldnames = ['NP', 'Grid', 'Total Cells', 'Total Time (s)', 'Efficiency (%)', 'Metric (s/step/cell)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for res in results:
            np_val = res['np']
            time = res['metrics']['total_time']
            # Default to 0 if regex failed
            metric = res['metrics'].get('norm_metric', 0.0)

            # Grid info
            N = res['metrics']['grid_dim']
            grid_str = f"{N}x{N}x{N}"
            total_cells = res['metrics']['total_cells']

            eff = (t_base / time) * 100.0

            # 1. Print to Console
            print(f"{np_val:<5} | {grid_str:<15} | {total_cells:<12} | {time:<10.4f} | {eff:<10.1f} | {metric:<.4e}")

            # 2. Write to CSV
            writer.writerow({
                'NP': np_val,
                'Grid': grid_str,
                'Total Cells': total_cells,
                'Total Time (s)': time,
                'Efficiency (%)': eff,
                'Metric (s/step/cell)': metric
            })

    print("="*85)
    print(f"[OUTPUT] Data logged to '{csv_filename}'")

    plot_results(results)

if __name__ == "__main__":
    main()