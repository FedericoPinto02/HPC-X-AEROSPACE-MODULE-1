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

# --- USER CONFIGURATION ---
EXECUTABLE = "../build/main"
CONFIG_FILE = "../data/config.json"
OUTPUT_DIR = "../output/"
# Save inside a folder named 'errorPlots' in the current directory (scripts)
ERROR_DIR = "errorPlots"
ERROR_VTK_DIR = os.path.join(OUTPUT_DIR, "ERROR_VTK")

# Simulation Parameters
NX = 50
DT = 0.001
T_END = 0.05
OUTPUT_FREQ = 1
DOMAIN_LEN_X = 6.0

# Enable L2 norm calculation and plotting
ENABLE_L2_ANALYSIS = True

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
    Calculates the Root Mean Square (RMS) error, consistent with plot_convergence_new.py.
    RMS = sqrt(mean(error^2))
    This normalizes the error by the number of points, making it independent of grid size.
    """
    return np.sqrt(np.mean(error_field**2))

def main():
    # 1. Setup directories
    if os.path.exists(ERROR_DIR):
        shutil.rmtree(ERROR_DIR)
    os.makedirs(ERROR_DIR, exist_ok=True)
    if os.path.exists(ERROR_VTK_DIR):
        shutil.rmtree(ERROR_VTK_DIR)
    os.makedirs(ERROR_VTK_DIR, exist_ok=True)

    # 2. Read and backup original config
    try:
        with open(CONFIG_FILE, 'r') as f:
            raw_original_content = f.read()
        orig_conf = json.loads(raw_original_content)
    except FileNotFoundError:
        print(f"Error: {CONFIG_FILE} not found.")
        sys.exit(1)

    print(f"--- Starting Error Analysis Run (Nx={NX}, dt={DT}) ---")

    try:
        # 3. Modify Config
        conf = copy.deepcopy(orig_conf)
        conf['mesh']['nx'] = NX
        conf['mesh']['input_for_manufactured_solution'] = True
        conf['time']['dt'] = DT
        conf['time']['t_end'] = T_END
        
        # Set unique base name to avoid mixing with other runs
        base_name = f"error_calc_n{NX}"
        conf['output']['base_filename'] = base_name
        conf['output']['output_frequency'] = OUTPUT_FREQ

        with open(CONFIG_FILE, 'w') as f:
            json.dump(conf, f, indent=4)

        # 4. Run Solver
        print("Running solver...")
        ret = subprocess.run([EXECUTABLE], capture_output=True, text=True)
        
        if ret.returncode != 0:
            print("Solver failed. STDERR:")
            print(ret.stderr)
            sys.exit(1)

        # 5. Post-process results
        h = DOMAIN_LEN_X / (NX + 0.5)
        
        # Find and sort valid files
        all_files = os.listdir(OUTPUT_DIR)
        vtk_files = [f for f in all_files if f.startswith(base_name) and f.endswith(".vtk")]
        vtk_files.sort(key=get_step_from_filename)

        print(f"Processing {len(vtk_files)} timesteps...")

        # Storage for L2 norms
        l2_history = {'t': [], 'u': [], 'v': [], 'w': [], 'p': [], 'mag': []}

        for i, f in enumerate(vtk_files):
            full_path = os.path.join(OUTPUT_DIR, f)
            
            # Extract time
            step = get_step_from_filename(f)
            t = step * DT

            # Read VTK
            grid = pv.read(full_path)
            points = grid.points
            
            # Get numerical fields
            vel = grid.point_data['velocity'] # Shape (N, 3)
            p = grid.point_data['pressure']   # Shape (N,)

            # Compute errors (Point-wise difference)
            err_u = vel[:, 0] - u_ana(points, t, h)
            err_v = vel[:, 1] - v_ana(points, t, h)
            err_w = vel[:, 2] - w_ana(points, t, h)
            err_p = p - p_ana(points, t)

            # Clear old data and add error fields
            grid.point_data.clear()
            grid.point_data['error_u'] = err_u
            grid.point_data['error_v'] = err_v
            grid.point_data['error_w'] = err_w
            grid.point_data['error_p'] = err_p
            
            # Add scalar magnitude of velocity error for easier visualization
            err_mag = np.sqrt(err_u**2 + err_v**2 + err_w**2)
            grid.point_data['error_vel_mag'] = err_mag

            # Save to error directory
            save_name = f"error_step_{step:05d}.vtk"
            grid.save(os.path.join(ERROR_VTK_DIR, save_name))

            # Calculate and store L2 norms if enabled
            if ENABLE_L2_ANALYSIS:
                # UPDATED: Use RMS (Root Mean Square) to match convergence plots
                # This divides by sqrt(N), making values independent of grid resolution
                l2_u = calculate_rms_error(err_u)
                l2_v = calculate_rms_error(err_v)
                l2_w = calculate_rms_error(err_w)
                l2_p = calculate_rms_error(err_p)
                l2_mag = calculate_rms_error(err_mag)

                l2_history['t'].append(t)
                l2_history['u'].append(l2_u)
                l2_history['v'].append(l2_v)
                l2_history['w'].append(l2_w)
                l2_history['p'].append(l2_p)
                l2_history['mag'].append(l2_mag)

                print(f"  Step {step} (t={t:.4f}) | RMS L2 Errors -> u: {l2_u:.2e}, v: {l2_v:.2e}, w: {l2_w:.2e}, p: {l2_p:.2e}, |V|: {l2_mag:.2e}")
            elif i % 10 == 0:
                print(f"  Processed step {step} (t={t:.4f})")

        # Plotting L2 norms
        if ENABLE_L2_ANALYSIS and len(l2_history['t']) > 0:
            print("\nPlotting RMS L2 error norms...")
            plt.figure(figsize=(10, 6))
            
            # Use different styles to distinguish overlapping lines (e.g. if u and v errors are identical)
            plt.plot(l2_history['t'], l2_history['u'], label='RMS Error u', linestyle='-', marker='o', markersize=3, alpha=0.8)
            plt.plot(l2_history['t'], l2_history['v'], label='RMS Error v', linestyle='--', marker='s', markersize=3, alpha=0.8)
            plt.plot(l2_history['t'], l2_history['w'], label='RMS Error w', linestyle='-.', marker='^', markersize=3, alpha=0.8)
            plt.plot(l2_history['t'], l2_history['p'], label='RMS Error p', linestyle=':', marker='x', markersize=4, alpha=0.8)
            
            # Plot magnitude with a distinct style
            plt.plot(l2_history['t'], l2_history['mag'], label='RMS Error |Vel|', color='black', linewidth=2, alpha=0.4)
            
            plt.xlabel('Time (s)')
            plt.ylabel('RMS L2 Error')
            plt.yscale('log')
            
            # Update title with detailed simulation parameters
            n_steps_plotted = len(l2_history['t'])
            plt.title(f'RMS L2 Error Evolution\nNx={NX}, dt={DT}, T_end={T_END}, Steps={n_steps_plotted}')
            
            plt.legend()
            plt.grid(True, which="both", ls="-", alpha=0.5)
            
            # Save plot with discretization details in filename
            # UPDATED: Included T_END and n_steps_plotted in the filename
            plot_filename = f"l2_error_Nx{NX}_dt{DT}_Tend{T_END}_steps{n_steps_plotted}.png"
            plot_path = os.path.join(ERROR_DIR, plot_filename)
            plt.savefig(plot_path)
            print(f"Plot saved to {plot_path}")
            plt.show()

        print(f"\nDone! Error files saved to: {ERROR_DIR}")
        print("You can now open this folder in ParaView.")

    except Exception as e:
        print(f"\nAn error occurred: {e}")

    finally:
        # 6. Restore original config
        print("Restoring original configuration...")
        with open(CONFIG_FILE, 'w') as f:
            f.write(raw_original_content)

if __name__ == "__main__":
    main()