import subprocess
import json
import os
import sys
import re
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing

# =============================================================================
# CONFIGURATION
# =============================================================================
EXECUTABLE  = "../build/main"
CONFIG_PATH = "../data/config.json"
OUTPUT_DIR  = "../output"

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

    config['grid']['Nx'] = new_N
    config['grid']['Ny'] = new_N
    config['grid']['Nz'] = new_N

    # Standard settings
    config['time']['dt'] = 0.001 # Small dt for stability on fine grids
    config['time']['t_end'] = config['time']['dt'] * NUM_STEPS
    config['output']['output_frequency'] = NUM_STEPS + 100 # Disable I/O
    config['output']['base_filename'] = "weak_test"

    with open(CONFIG_PATH, 'w') as f:
        json.dump(config, f, indent=4)

    total_cells = new_N**3
    print(f"[SETUP NP={n_procs}] Grid: {new_N}x{new_N}x{new_N} | Total Cells: {total_cells} ({total_cells/n_procs:.0f}/core)")
    return total_cells

def parse_output(stdout_str):
    metrics = {}
    # Extract Total CPU Time
    match = re.search(r"Total CPU Time \(solve loop\):\s+([\d\.]+)", stdout_str)
    if match: metrics['total_time'] = float(match.group(1))
    return metrics

def run_test(n_procs):
    update_config_weak(n_procs)

    cmd = [EXECUTABLE]
    if n_procs > 1:
        cmd = ["mpirun", "-np", str(n_procs)] + cmd

    print(f"  Running... ", end="", flush=True)
    try:
        # Run and capture output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        metrics = parse_output(result.stdout)
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
    plt_filename = "../results/weak_scalability_results.png"
    plt.savefig(plt_filename)
    print(f"\n[OUTPUT] Plot saved to '{plt_filename}'")

def main():
    print("============================================================")
    print(" WEAK SCALABILITY TEST SUITE")
    print("============================================================")

    results = []

    for np in PROCS_TO_TEST:
        metrics = run_test(np)
        if metrics:
            results.append({'np': np, 'metrics': metrics})

    if not results: return

    # Print Table
    print("\n" + "="*60)
    print(f"{'NP':<5} | {'Grid':<15} | {'Time (s)':<10} | {'Efficiency (%)':<15}")
    print("-" * 60)

    t_base = results[0]['metrics']['total_time']

    for res in results:
        np_val = res['np']
        time = res['metrics']['total_time']

        # Calculate grid size roughly for display
        scale = np_val ** (1.0/3.0)
        N = int(round(BASE_GRID_SIZE * scale))
        if N%2!=0: N+=1
        grid_str = f"{N}x{N}x{N}"

        eff = (t_base / time) * 100.0

        print(f"{np_val:<5} | {grid_str:<15} | {time:<10.4f} | {eff:<15.1f}")
    print("="*60)

    plot_results(results)

if __name__ == "__main__":
    main()