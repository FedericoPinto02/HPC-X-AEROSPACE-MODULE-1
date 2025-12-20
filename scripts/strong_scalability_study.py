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

# List of processor counts to test.
# Ensure you don't exceed your physical core count.
# Example: [1, 2, 4, 8, 16]
MAX_CORES = multiprocessing.cpu_count()
PROCS_TO_TEST = [1, 2, 4, 8]
if 16 <= MAX_CORES: PROCS_TO_TEST.append(16)

# Simulation settings for the test
NUM_STEPS = 50       # Run enough steps to average out initialization jitter
GRID_SIZE = 64       # Nx=Ny=Nz. Make this large enough to justify parallelism!
# =============================================================================

def setup_config():
    """
    Modifies config.json to:
    1. Set a reasonable grid size.
    2. Disable heavy I/O (VTK writing) so we measure CPU scaling, not Disk I/O.
    """
    with open(CONFIG_PATH, 'r') as f:
        config = json.load(f)

    # Set Grid (Strong Scaling: Grid size stays constant)
    config['mesh']['nx'] = GRID_SIZE
    config['mesh']['ny'] = GRID_SIZE
    config['mesh']['nz'] = GRID_SIZE

    # Set Time
    config['time']['t_end'] = config['time']['dt'] * NUM_STEPS

    # DISABLE I/O: Writing files kills scalability metrics.
    # We set frequency higher than total steps.
    config['output']['output_frequency'] = NUM_STEPS + 100
    config['output']['base_filename'] = "perf_test"

    with open(CONFIG_PATH, 'w') as f:
        json.dump(config, f, indent=4)

    print(f"[SETUP] Grid: {GRID_SIZE}^3 | Steps: {NUM_STEPS} | I/O: Disabled")

def parse_output(stdout_str):
    """
    Parses the specific output format provided by the user.
    """
    metrics = {}

    # Regex patterns based on your example
    # Total CPU Time (solve loop):  2.641829 s
    re_total = r"Total CPU Time \(solve loop\):\s+([\d\.]+)"
    # Avg CPU Time per Timestep:    0.132091 s
    re_avg   = r"Avg CPU Time per Timestep:\s+([\d\.]+)"
    # CPU Time / (Steps * Cells):   2.1135e-06 s
    re_norm  = r"CPU Time / \(Steps \* Cells\):\s+([\d\.eE\+\-]+)"

    match_total = re.search(re_total, stdout_str)
    match_avg   = re.search(re_avg, stdout_str)
    match_norm  = re.search(re_norm, stdout_str)

    if match_total: metrics['total_time'] = float(match_total.group(1))
    if match_avg:   metrics['avg_step']   = float(match_avg.group(1))
    if match_norm:  metrics['norm_metric']= float(match_norm.group(1))

    return metrics

def run_test(n_procs):
    cmd = [EXECUTABLE]
    if n_procs > 1:
        cmd = ["mpirun", "-np", str(n_procs)] + cmd

    print(f"  Running with NP = {n_procs} ... ", end="", flush=True)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        metrics = parse_output(result.stdout)
        print(f"Done. (Time: {metrics.get('total_time', 'N/A')}s)")
        return metrics
    except subprocess.CalledProcessError as e:
        print("FAILED.")
        print(e.stderr)
        return None

def plot_results(results):
    procs = [r['np'] for r in results]
    times = [r['metrics']['total_time'] for r in results]

    # Calculate Speedup and Efficiency
    t_serial = times[0]
    speedup = [t_serial / t for t in times]
    ideal_speedup = procs
    efficiency = [s / p * 100.0 for s, p in zip(speedup, procs)]

    # --- Plot 1: Speedup ---
    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.plot(procs, speedup, 'o-', label='Measured Speedup', linewidth=2)
    plt.plot(procs, ideal_speedup, 'k--', label='Ideal Linear', alpha=0.5)
    plt.title(f"Strong Scaling (Grid {GRID_SIZE}^3)")
    plt.xlabel("Number of Processors")
    plt.ylabel("Speedup (T1 / Tn)")
    plt.grid(True, which="both", ls="-", alpha=0.4)
    plt.legend()
    plt.xticks(procs)

    # --- Plot 2: Parallel Efficiency ---
    plt.subplot(1, 2, 2)
    plt.plot(procs, efficiency, 'rs-', label='Efficiency %', linewidth=2)
    plt.axhline(100, color='k', linestyle='--', alpha=0.5)
    plt.title("Parallel Efficiency")
    plt.xlabel("Number of Processors")
    plt.ylabel("Efficiency (%)")
    plt.ylim(0, 110)
    plt.grid(True, which="both", ls="-", alpha=0.4)
    plt.xticks(procs)

    plt.tight_layout()
    plt_filename = "../results/strong_scalability_results.png"
    plt.savefig(plt_filename)
    print(f"\n[OUTPUT] Plot saved to '{plt_filename}'")

def main():
    print("============================================================")
    print(" SCALABILITY TEST SUITE")
    print("============================================================")

    setup_config()

    results = []

    print("\n[STARTING TESTS]")
    for np in PROCS_TO_TEST:
        metrics = run_test(np)
        if metrics:
            results.append({
                'np': np,
                'metrics': metrics
            })

    if not results:
        print("No results collected.")
        sys.exit(1)

    # Print Table
    print("\n" + "="*75)
    print(f"{'NP':<5} | {'Total Time (s)':<15} | {'Speedup':<10} | {'Efficiency (%)':<15} | {'Metric (s/step/cell)':<20}")
    print("-" * 75)

    t_base = results[0]['metrics']['total_time']

    for res in results:
        np_val = res['np']
        time = res['metrics']['total_time']
        metric = res['metrics']['norm_metric']

        speedup = t_base / time
        eff = (speedup / np_val) * 100.0

        print(f"{np_val:<5} | {time:<15.4f} | {speedup:<10.2f} | {eff:<15.1f} | {metric:<20.4e}")
    print("="*75)

    plot_results(results)

if __name__ == "__main__":
    main()