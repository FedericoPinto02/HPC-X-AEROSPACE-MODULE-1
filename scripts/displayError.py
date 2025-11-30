#!/usr/bin/env python3
import json
import copy
import subprocess
import os
import shutil
import re
import sys
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from scipy.stats import skew, kurtosis
import time

# --- USER CONFIGURATION ---
EXECUTABLE = "../build/main"
CONFIG_FILE = "../data/config.json"
FINAL_OUTPUT_DIR = "../output/"

# Subdirectories for error analysis
ERROR_DIR = "errorPlots"
PRESSURE_ANALYSIS_DIR = os.path.join(ERROR_DIR, "pressure_detailed") # New folder
ERROR_VTK_DIR = os.path.join(FINAL_OUTPUT_DIR, "ERROR_VTK")

# =============================================================================
# === SIMULATION PARAMETERS ===
# =============================================================================
NX_LIST    = [30]
T_END = 0.5 

NX_REF = 50.0
dx = 6.0 / (NX_REF + 0.5)
nu = 6
DT_REF = dx**2 / nu * 0.1 
DT_LIST = [DT_REF * (NX_REF / nx)**2 for nx in NX_LIST]

T_END_LIST = np.full_like(NX_LIST, T_END, dtype=float)

OUTPUT_FREQ = 1
DOMAIN_LEN_X = 6.0
ENABLE_L2_ANALYSIS = True

# =============================================================================
# === ANALYTICAL SOLUTIONS ===
# =============================================================================
def u_ana(points, t, h): 
    return np.sin(t) * np.sin(points[:,0] + h/2) * np.sin(points[:,1]) * np.sin(points[:,2])

def v_ana(points, t, h): 
    return np.sin(t) * np.cos(points[:,0]) * np.cos(points[:,1] + h/2) * np.cos(points[:,2])

def w_ana(points, t, h): 
    return np.sin(t) * np.cos(points[:,0]) * np.sin(points[:,1]) * (np.sin(points[:,2] + h/2) + np.cos(points[:,2] + h/2))

def p_ana(points, t): 
    return -3 * np.sin(t) * np.cos(points[:,0]) * np.sin(points[:,1]) * (np.sin(points[:,2]) - np.cos(points[:,2]))

def get_step_from_filename(filename):
    match = re.search(r'_(\d+)\.vtk$', filename)
    return int(match.group(1)) if match else -1

def calculate_rms_error(error_field):
    return np.sqrt(np.mean(error_field**2))

# =============================================================================
# === NEW: DETAILED PRESSURE ERROR ANALYSIS ===
# =============================================================================
def analyze_pressure_detailed(points, err_p, nx, dt, step, run_id):
    """
    Performs statistical and spatial analysis of pressure error.
    Splits data into Boundary vs Center and computes directional moments.
    """
    os.makedirs(PRESSURE_ANALYSIS_DIR, exist_ok=True)
    
    # 1. Geometry Setup to distinguish Boundary vs Center
    # We assume a cubic/rectangular domain. We define "Boundary" as the outer 10%
    min_coords = np.min(points, axis=0)
    max_coords = np.max(points, axis=0)
    lengths = max_coords - min_coords
    margin = lengths * 0.10 # 10% margin from edges
    
    # Create masks
    mask_x = (points[:,0] < (min_coords[0] + margin[0])) | (points[:,0] > (max_coords[0] - margin[0]))
    mask_y = (points[:,1] < (min_coords[1] + margin[1])) | (points[:,1] > (max_coords[1] - margin[1]))
    mask_z = (points[:,2] < (min_coords[2] + margin[2])) | (points[:,2] > (max_coords[2] - margin[2]))
    
    mask_boundary = mask_x | mask_y | mask_z
    mask_center = ~mask_boundary
    
    err_p_bnd = err_p[mask_boundary]
    err_p_ctr = err_p[mask_center]
    
    # 2. Statistical Moments
    stats = {
        'Global':   [np.mean(err_p), np.std(err_p), skew(err_p), kurtosis(err_p)],
        'Boundary': [np.mean(err_p_bnd), np.std(err_p_bnd), skew(err_p_bnd), kurtosis(err_p_bnd)],
        'Center':   [np.mean(err_p_ctr), np.std(err_p_ctr), skew(err_p_ctr), kurtosis(err_p_ctr)]
    }
    
    # Print Report
    report_file = os.path.join(PRESSURE_ANALYSIS_DIR, f"stats_report_{run_id}.txt")
    with open(report_file, 'w') as f:
        f.write(f"PRESSURE ERROR ANALYSIS | Run: {run_id} | Step: {step}\n")
        f.write("="*80 + "\n")
        f.write(f"{'Region':<12} | {'Mean':<12} | {'Std Dev':<12} | {'Skewness':<12} | {'Kurtosis':<12}\n")
        f.write("-" * 80 + "\n")
        for region, vals in stats.items():
            f.write(f"{region:<12} | {vals[0]:<12.4e} | {vals[1]:<12.4e} | {vals[2]:<12.4f} | {vals[3]:<12.4f}\n")
        f.write("-" * 80 + "\n")
        f.write("NOTE:\n")
        f.write("  Skewness != 0 implies asymmetry (bias).\n")
        f.write("  Kurtosis > 3 implies heavy tails (outliers drive the error).\n")
    
    # 3. Plotting Histograms
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.hist(err_p, bins=50, density=True, alpha=0.6, color='gray', label='Global')
    plt.title(f'Global Pressure Error PDF\nRun {run_id}')
    plt.xlabel('Error Value (Pa)')
    plt.ylabel('Density')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.hist(err_p_ctr, bins=50, density=True, alpha=0.5, color='blue', label='Center')
    plt.hist(err_p_bnd, bins=50, density=True, alpha=0.5, color='red', label='Boundary')
    plt.title('Boundary vs Center Comparison')
    plt.xlabel('Error Value (Pa)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(PRESSURE_ANALYSIS_DIR, f"hist_dist_{run_id}.png"))
    plt.close()

    # 4. Directional Analysis (Slicing/Binning)
    # Binning data along axes to see spatial trends
    n_bins = 20
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    dirs = ['X', 'Y', 'Z']
    
    for i, ax_idx in enumerate([0, 1, 2]):
        coord = points[:, ax_idx]
        # Create bins
        bins = np.linspace(coord.min(), coord.max(), n_bins + 1)
        centers = (bins[:-1] + bins[1:]) / 2
        
        # Digitize coordinates to find which bin they belong to
        bin_indices = np.digitize(coord, bins) - 1
        
        # Calculate mean error and std dev in each bin
        bin_means = []
        bin_stds = []
        
        for b in range(n_bins):
            # Select points in this bin
            mask_bin = bin_indices == b
            if np.any(mask_bin):
                bin_vals = err_p[mask_bin]
                bin_means.append(np.mean(bin_vals))
                bin_stds.append(np.std(bin_vals))
            else:
                bin_means.append(0)
                bin_stds.append(0)
        
        # Plot
        axes[i].errorbar(centers, bin_means, yerr=bin_stds, fmt='-o', capsize=3, label='Mean Error')
        axes[i].set_title(f'Mean Error along {dirs[i]}')
        axes[i].set_xlabel(f'{dirs[i]} Coordinate')
        axes[i].set_ylabel('Mean Pressure Error')
        axes[i].grid(True, alpha=0.4)
        axes[i].axhline(0, color='k', linestyle='--', linewidth=0.8)

    plt.suptitle(f"Directional Error Profiles - Run {run_id}", fontsize=14)
    plt.savefig(os.path.join(PRESSURE_ANALYSIS_DIR, f"directional_profiles_{run_id}.png"))
    plt.close()
    
    print(f"    --> Detailed pressure analysis saved to: {PRESSURE_ANALYSIS_DIR}")


# =============================================================================
# === MAIN RUN LOOP ===
# =============================================================================

def run_single_case(nx, dt, t_end, original_config_content, run_index):
    dt_str = str(dt).replace('.', 'p')
    run_id = f"n{nx}_dt{dt_str}"
    
    print(f"\n{'='*60}")
    print(f" [Run {run_index+1}] Starting simulation: Nx={nx}, dt={dt}, T_end={t_end}")
    print(f"{'='*60}")

    try:
        conf = json.loads(original_config_content)
        conf['mesh']['nx'] = nx
        conf['mesh']['input_for_manufactured_solution'] = True
        conf['time']['dt'] = dt
        conf['time']['t_end'] = t_end
        base_name = f"error_{run_id}"
        conf['output']['base_filename'] = base_name
        conf['output']['output_frequency'] = OUTPUT_FREQ
        
        with open(CONFIG_FILE, 'w') as f:
            json.dump(conf, f, indent=4)

        # Cleanup
        for f in os.listdir(FINAL_OUTPUT_DIR):
            if f.startswith(base_name) and f.endswith(".vtk"):
                try: os.remove(os.path.join(FINAL_OUTPUT_DIR, f))
                except OSError: pass

        # Run Solver
        print("  -> Running solver...")
        ret = subprocess.run([EXECUTABLE], capture_output=True, text=True)
        if ret.returncode != 0:
            print("  !!! Solver FAILED. STDERR:")
            print(ret.stderr)
            return None

        # Post-Processing
        h = DOMAIN_LEN_X / (nx + 0.5)
        all_files = os.listdir(FINAL_OUTPUT_DIR)
        vtk_files = [os.path.join(FINAL_OUTPUT_DIR, f) for f in all_files 
                     if f.startswith(base_name) and f.endswith(".vtk")]
        vtk_files.sort(key=lambda x: get_step_from_filename(os.path.basename(x)))

        print(f"  -> Processing {len(vtk_files)} timesteps...")
        
        l2_history = {'t': [], 'u': [], 'v': [], 'w': [], 'p': [], 'mag': []} 
        os.makedirs(ERROR_VTK_DIR, exist_ok=True)
        os.makedirs(ERROR_DIR, exist_ok=True)

        for i, f_path in enumerate(vtk_files):
            step = get_step_from_filename(os.path.basename(f_path))
            t = step * dt
            
            grid = pv.read(f_path)
            points = grid.points
            
            vel = grid.point_data['velocity']
            p = grid.point_data['pressure']
            
            u_ex = u_ana(points, t, h)
            v_ex = v_ana(points, t, h)
            w_ex = w_ana(points, t, h)
            p_ex = p_ana(points, t)
            
            err_u = vel[:, 0] - u_ex
            err_v = vel[:, 1] - v_ex
            err_w = vel[:, 2] - w_ex
            err_p = p - p_ex
            err_mag = np.sqrt(err_u**2 + err_v**2 + err_w**2)
            
            # Save Error VTK
            grid.point_data.clear()
            grid.point_data['error_vel_mag'] = err_mag
            grid.point_data['error_p'] = err_p
            grid.save(os.path.join(ERROR_VTK_DIR, f"error_{run_id}_step_{step:05d}.vtk"))
            
            if ENABLE_L2_ANALYSIS:
                l2_history['t'].append(t)
                l2_history['u'].append(calculate_rms_error(err_u))
                l2_history['v'].append(calculate_rms_error(err_v))
                l2_history['w'].append(calculate_rms_error(err_w))
                l2_history['p'].append(calculate_rms_error(err_p))
                l2_history['mag'].append(calculate_rms_error(err_mag))
            
            # --- TRIGGER DETAILED PRESSURE ANALYSIS ON THE FINAL TIMESTEP ---
            is_last_step = (i == len(vtk_files) - 1)
            if is_last_step:
                print(f"  -> Performing detailed pressure analysis for step {step}...")
                analyze_pressure_detailed(points, err_p, nx, dt, step, run_id)

        # Save History Data and Plot
        if l2_history['t']:
            txt_path = os.path.join(ERROR_DIR, f"l2_error_data_Nx{nx}_dt{dt}.txt")
            with open(txt_path, 'w') as f:
                f.write(f"# Nx={nx}, dt={dt}\n# Time | RMS_p\n")
                for j in range(len(l2_history['t'])):
                    f.write(f"{l2_history['t'][j]:.6f} | {l2_history['p'][j]:.6e}\n")
            
            plt.figure(figsize=(10, 6))
            plt.plot(l2_history['t'], l2_history['p'], label='RMS Error p', color='red')
            plt.yscale('log')
            plt.title(f'Pressure RMS Error (Nx={nx})')
            plt.savefig(os.path.join(ERROR_DIR, f"l2_plot_P_Nx{nx}_dt{dt}.png"))
            plt.close()

        return {'nx': nx, 'dt': dt, 't': l2_history['t'], 'mag': l2_history['mag'], 'success': True}

    except Exception as e:
        print(f"  !!! Error in Run Nx={nx}: {e}")
        return {'success': False, 'nx': nx, 'dt': dt}
    finally:
        with open(CONFIG_FILE, 'w') as f: f.write(original_config_content)

def main():
    print("--- SIMULATION & ANALYSIS ---")
    os.makedirs(FINAL_OUTPUT_DIR, exist_ok=True)
    os.makedirs(ERROR_DIR, exist_ok=True)
    
    try:
        with open(CONFIG_FILE, 'r') as f: original_config = f.read()
    except FileNotFoundError:
        sys.exit("Config file not found.")

    tasks = list(zip(NX_LIST, DT_LIST, T_END_LIST))
    for i, (nx, dt, tend) in enumerate(tasks):
        run_single_case(nx, dt, tend, original_config, i)

if __name__ == "__main__":
    main()