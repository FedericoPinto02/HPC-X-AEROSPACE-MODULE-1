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

# --- USER CONFIGURATION ---
EXECUTABLE = "../build/main"
CONFIG_FILE = "../data/config.json"
OUTPUT_DIR = "../output/"
ERROR_DIR = os.path.join(OUTPUT_DIR, "error_analysis")

# Simulation Parameters
NX = 50
DT = 0.001
T_END = 0.1
OUTPUT_FREQ = 1
DOMAIN_LEN_X = 6.0

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

def main():
    # 1. Setup directories
    if os.path.exists(ERROR_DIR):
        shutil.rmtree(ERROR_DIR)
    os.makedirs(ERROR_DIR, exist_ok=True)

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

            # Compute errors
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
            grid.point_data['error_vel_mag'] = np.sqrt(err_u**2 + err_v**2 + err_w**2)

            # Save to error directory
            save_name = f"error_step_{step:05d}.vtk"
            grid.save(os.path.join(ERROR_DIR, save_name))

            if i % 10 == 0:
                print(f"  Processed step {step} (t={t:.4f})")

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