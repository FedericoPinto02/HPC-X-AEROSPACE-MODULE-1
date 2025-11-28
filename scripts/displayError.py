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
import time # Time is still useful for tracking total execution time

# --- USER CONFIGURATION ---
# Paths are relative to the 'scripts' directory
EXECUTABLE = "../build/main"
CONFIG_FILE = "../data/config.json"
FINAL_OUTPUT_DIR = "../output/"

# Subdirectories for error analysis
ERROR_DIR = "errorPlots"
ERROR_VTK_DIR = os.path.join(FINAL_OUTPUT_DIR, "ERROR_VTK")

# =============================================================================
# === SIMULATION PARAMETERS (VECTORS) ===
# =============================================================================
# The lists below MUST have the same length (coupled lists mode is mandatory).
# The script will execute: Run 0 (Nx[0], Dt[0]...), Run 1 (Nx[1], Dt[1]...), etc.

NX_LIST    = [20, 30, 40, 50]
T_END = 2.0  # Fixed end time for all runs

NX_REF = 50.0
dx = 6.0 / (NX_REF + 0.5)
nu = 6
DT_REF = dx**2 / nu * 0.1  # Stability condition for diffusion equation
DT_LIST = [DT_REF * (NX_REF / nx)**2 for nx in NX_LIST]

T_END_LIST = np.full_like(NX_LIST, T_END, dtype=float)

OUTPUT_FREQ = 1

# Other fixed parameters
DOMAIN_LEN_X = 6.0
ENABLE_L2_ANALYSIS = True

# Removed Multicore settings, running serially now.
# =============================================================================

# --- ANALYTICAL SOLUTIONS (Staggered Grid) ---
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
    """
    Calculates the Root Mean Square (RMS) error.
    RMS = sqrt(mean(error^2))
    """
    return np.sqrt(np.mean(error_field**2))

def run_single_case(nx, dt, t_end, original_config_content, run_index):
    """
    Runs a single simulation case serially. Modifies the config file, runs the solver, 
    and restores the config file content afterward.
    """
    
    # Unique identifier for this run
    dt_str = str(dt).replace('.', 'p')
    run_id = f"n{nx}_dt{dt_str}"
    
    print(f"\n{'='*60}")
    print(f" [Run {run_index+1}] Starting simulation: Nx={nx}, dt={dt}, T_end={t_end}")
    print(f"{'='*60}")

    # 1. Prepare Config and File Paths
    try:
        # Load the original content from string (safer than reading the file)
        conf = json.loads(original_config_content)
        
        conf['mesh']['nx'] = nx
        conf['mesh']['input_for_manufactured_solution'] = True
        conf['time']['dt'] = dt
        conf['time']['t_end'] = t_end
        
        base_name = f"error_{run_id}"
        conf['output']['base_filename'] = base_name
        conf['output']['output_frequency'] = OUTPUT_FREQ
        
        # Write the modified config back to the main file
        with open(CONFIG_FILE, 'w') as f:
            json.dump(conf, f, indent=4)

        # 2. Cleanup OLD simulation files for THIS specific config
        print("  -> Cleaning up old simulation files for this config...")
        for f in os.listdir(FINAL_OUTPUT_DIR):
            if f.startswith(base_name) and f.endswith(".vtk"):
                try:
                    os.remove(os.path.join(FINAL_OUTPUT_DIR, f))
                except OSError:
                    pass

        # 3. Run Solver (Executed from the current script directory)
        print("  -> Running solver...")
        # The solver path is relative to the current script directory
        ret = subprocess.run([EXECUTABLE], capture_output=True, text=True)
        
        if ret.returncode != 0:
            print("  !!! Solver FAILED. STDERR:")
            print(ret.stderr)
            return None

        # 4. Post-Processing (Analysis)
        h = DOMAIN_LEN_X / (nx + 0.5)
        
        # Filter files for this specific run in the FINAL_OUTPUT_DIR
        all_files = os.listdir(FINAL_OUTPUT_DIR)
        vtk_files = [os.path.join(FINAL_OUTPUT_DIR, f) for f in all_files 
                     if f.startswith(base_name) and f.endswith(".vtk")]
        vtk_files.sort(key=lambda x: get_step_from_filename(os.path.basename(x)))

        print(f"  -> Processing {len(vtk_files)} timesteps...")
        
        l2_history = {'t': [], 'u': [], 'v': [], 'w': [], 'p': [], 'mag': []} 
        
        os.makedirs(ERROR_VTK_DIR, exist_ok=True)
        os.makedirs(ERROR_DIR, exist_ok=True)

        for f_path in vtk_files:
            step = get_step_from_filename(os.path.basename(f_path))
            t = step * dt
            
            grid = pv.read(f_path)
            points = grid.points
            
            vel = grid.point_data['velocity']
            p = grid.point_data['pressure']
            
            # Compute analytical solution
            u_ex = u_ana(points, t, h)
            v_ex = v_ana(points, t, h)
            w_ex = w_ana(points, t, h)
            p_ex = p_ana(points, t)
            
            # Compute point-wise errors
            err_u = vel[:, 0] - u_ex
            err_v = vel[:, 1] - v_ex
            err_w = vel[:, 2] - w_ex
            err_p = p - p_ex
            err_mag = np.sqrt(err_u**2 + err_v**2 + err_w**2)
            
            # Save Error VTK
            grid.point_data.clear()
            grid.point_data['error_vel_mag'] = err_mag
            grid.point_data['error_p'] = err_p
            
            save_name = f"error_{run_id}_step_{step:05d}.vtk"
            grid.save(os.path.join(ERROR_VTK_DIR, save_name))
            
            if ENABLE_L2_ANALYSIS:
                # Calculate RMS L2 norms
                l2_history['t'].append(t)
                l2_history['u'].append(calculate_rms_error(err_u))
                l2_history['v'].append(calculate_rms_error(err_v))
                l2_history['w'].append(calculate_rms_error(err_w))
                l2_history['p'].append(calculate_rms_error(err_p))
                l2_history['mag'].append(calculate_rms_error(err_mag))

        # --- Save Individual Run Data and Plot ---
        if l2_history['t']:
            
            # A. Save Text Data (All components)
            txt_filename = f"l2_error_data_Nx{nx}_dt{dt}.txt"
            txt_path = os.path.join(ERROR_DIR, txt_filename)
            
            with open(txt_path, 'w') as f:
                f.write(f"# Error Analysis Data Log\n")
                f.write(f"# Nx={nx}, dt={dt}, T_end={t_end}\n")
                f.write(f"# {'Time (s)':>12} | {'RMS_u':>14} | {'RMS_v':>14} | {'RMS_w':>14} | {'RMS_p':>14} | {'RMS_Mag':>14}\n")
                f.write("-" * 105 + "\n")
                
                for j in range(len(l2_history['t'])):
                    f.write(f"  {l2_history['t'][j]:12.6f} | "
                            f"{l2_history['u'][j]:14.6e} | "
                            f"{l2_history['v'][j]:14.6e} | "
                            f"{l2_history['w'][j]:14.6e} | "
                            f"{l2_history['p'][j]:14.6e} | "
                            f"{l2_history['mag'][j]:14.6e}\n")
            
            # B. Plotting (All components vs Time)
            plt.figure(figsize=(10, 6))
            
            plt.plot(l2_history['t'], l2_history['u'], label='RMS Error u', linestyle='-', marker='o', markersize=3, alpha=0.8)
            plt.plot(l2_history['t'], l2_history['v'], label='RMS Error v', linestyle='--', marker='s', markersize=3, alpha=0.8)
            plt.plot(l2_history['t'], l2_history['w'], label='RMS Error w', linestyle='-.', marker='^', markersize=3, alpha=0.8)
            plt.plot(l2_history['t'], l2_history['p'], label='RMS Error p', linestyle=':', marker='x', markersize=4, alpha=0.8)
            plt.plot(l2_history['t'], l2_history['mag'], label='RMS Error |Vel|', color='black', linewidth=2, alpha=0.4)
            
            plt.xlabel('Time (s)')
            plt.ylabel('RMS L2 Error (log scale)')
            plt.yscale('log')
            plt.title(f'RMS L2 Error Evolution (Nx={nx}, dt={dt})')
            plt.legend()
            plt.grid(True, which="both", ls="-", alpha=0.5)
            
            plot_filename = f"l2_plot_all_Nx{nx}_dt{dt}.png"
            plt.savefig(os.path.join(ERROR_DIR, plot_filename))
            plt.close()
            
            print(f"  -> Plot saved to: {plot_filename} and data to: {txt_filename}")

        # 5. Return results for the final comparison plot
        return {
            'nx': nx,
            'dt': dt,
            't': l2_history['t'],
            'mag': l2_history['mag'],
            'success': True
        }

    except Exception as e:
        print(f"  !!! An error occurred during run Nx={nx}: {e}")
        return {'success': False, 'nx': nx, 'dt': dt}

    finally:
        # Restore original config file content
        with open(CONFIG_FILE, 'w') as f:
            f.write(original_config_content)

def main():
    print("--- SERIAL BATCH SIMULATION SCRIPT ---")
    
    # 1. Validation
    if not (len(NX_LIST) == len(DT_LIST) == len(T_END_LIST)):
        print("ERROR: NX_LIST, DT_LIST, and T_END_LIST must have the same length (coupled mode).")
        sys.exit(1)

    # 2. Prepare directories
    os.makedirs(FINAL_OUTPUT_DIR, exist_ok=True)
    os.makedirs(ERROR_DIR, exist_ok=True)
    os.makedirs(ERROR_VTK_DIR, exist_ok=True)
    
    # 3. Read and store original config content once
    try:
        with open(CONFIG_FILE, 'r') as f:
            original_config_content = f.read()
    except FileNotFoundError:
        print(f"Error: {CONFIG_FILE} not found.")
        sys.exit(1)

    # 4. Run Serially
    start_time = time.time()
    results = []
    
    tasks = list(zip(NX_LIST, DT_LIST, T_END_LIST))
    
    for i, (nx, dt, tend) in enumerate(tasks):
        res = run_single_case(nx=nx, dt=dt, t_end=tend, 
                              original_config_content=original_config_content, 
                              run_index=i)
        if res:
            results.append(res)
        
    total_time = time.time() - start_time
    print(f"\nAll serial simulations completed in {total_time:.2f} seconds.")

    # 5. Aggregate Results and Plot (Comparison Plot: Mag Error only)
    print("\n--- Generating Final Comparison Plot ---")
    
    plt.figure(figsize=(10, 6))
    
    for res in results:
        if not res or not res['success']:
            # Already printed error message in run_single_case
            continue
            
        nx = res['nx']
        dt = res['dt']
        t = res['t']
        mag = res['mag']
        
        # Add to plot (using only magnitude for comparison plot clarity)
        plt.plot(t, mag, label=f'Nx={nx}, dt={dt}', marker='o', markersize=3, alpha=0.7)

    plt.xlabel('Time (s)')
    plt.ylabel('RMS Velocity Magnitude Error (log scale)')
    plt.yscale('log')
    plt.title('Convergence Analysis (Serial Runs) - Velocity Magnitude Comparison')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.3)
    
    final_plot_path = os.path.join(ERROR_DIR, "convergence_comparison_Vmag.png")
    plt.savefig(final_plot_path)
    print(f"Comparison plot saved to: {final_plot_path}")

if __name__ == "__main__":
    main()