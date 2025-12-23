import json
import subprocess
import os
import sys
import copy
import math

# Attempt to import the plotting script
try:
    import plot_convergence_new
except ImportError:
    print("Error: Could not find 'plot_convergence.py' in the same directory.")
    sys.exit(1)

# --- GLOBAL CONFIGURATION (TO BE EDITED BY USER) ---
EXECUTABLE_PATH = "../build/main"      # Path to your C++ executable
CONFIG_PATH = "../data/config.json"    # Path to the JSON config file
OUTPUT_DIR = "../output/"              # Output directory (must match JSON)

# --- STUDY MODE SELECTION (EDIT THIS LINE) ---
# Choose one of the following modes:
# 1. SPATIAL_ONLY: Fix dt, vary Nx. Studies spatial order (Expected O(h^2)).
# 2. TEMPORAL_ONLY: Fix Nx (fine grid), vary dt. Studies temporal order (Expected O(dt^2)).
# 3. SPACE_TIME_COUPLED: Vary both (dt ~ 1/Nx). Maintains constant CFL (Expected O(h^2) or O(dt^2)).
STUDY_MODE = "SPACE_TIME_COUPLED" 
#STUDY_MODE = "SPATIAL_ONLY" 
#STUDY_MODE = "TEMPORAL_ONLY"

# --- BASELINE & RESOLUTION PARAMETERS ---
BASE_NX = 20           # Reference resolution
BASE_DT = 0.05        # Reference timestep
T_END = 1.5           # Fixed physical end time for ALL simulations

# --- RESOLUTION SETS ---
# Set used when STUDY_MODE is SPATIAL_ONLY or SPACE_TIME_COUPLED
NX_VALUES_FOR_SPATIAL_STUDY = [20, 25, 30, 35, 40, 45, 50, 60, 80, 100, 120] 

# Set used when STUDY_MODE is TEMPORAL_ONLY
# Multipliers relative to BASE_DT (e.g., 0.5 means dt = 0.0005)
DT_MULTIPLIERS_FOR_TEMPORAL_STUDY = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
# BASE_DT = 0.2  
# T_END = 1.0
FIXED_NX_FOR_TEMPORAL_STUDY = 70


def calculate_params_for_mode(loop_val, base_dt, study_mode):
    """Calculates the actual Nx, dt, and number of steps based on the study mode."""
    
    nx = BASE_NX 
    dt = base_dt
    
    if study_mode == "SPATIAL_ONLY":
        nx = loop_val
        dt = base_dt # Fixed dt
    
    elif study_mode == "SPACE_TIME_COUPLED":
        nx = loop_val
        # dt scales as 1/Nx (constant CFL for diffusion/viscous terms)
        dt = base_dt * (BASE_NX / nx)

    elif study_mode == "TEMPORAL_ONLY":
        # loop_val is the dt multiplier
        nx = FIXED_NX_FOR_TEMPORAL_STUDY
        dt = base_dt * loop_val
        
    else:
        raise ValueError(f"Unknown STUDY_MODE: {study_mode}")

    # Calculate exact number of steps and actual dt to hit T_END
    num_steps = int(round(T_END / dt))
    
    # Ensure num_steps is at least 1
    if num_steps < 1:
        num_steps = 1
        dt_actual = T_END 
    else:
        dt_actual = T_END / num_steps

    return nx, dt_actual, num_steps

def main():
    # 1. Read the original config
    try:
        with open(CONFIG_PATH, 'r') as f:
            raw_original_content = f.read()
        original_config = json.loads(raw_original_content)
    except FileNotFoundError:
        print(f"Error: Configuration file not found at {CONFIG_PATH}")
        sys.exit(1)

    simulation_data_for_plot = []
    
    # Select the loop variable based on the study mode
    if STUDY_MODE in ["SPATIAL_ONLY", "SPACE_TIME_COUPLED"]:
        loop_values = NX_VALUES_FOR_SPATIAL_STUDY
        loop_name = "Nx"
    elif STUDY_MODE == "TEMPORAL_ONLY":
        loop_values = DT_MULTIPLIERS_FOR_TEMPORAL_STUDY
        loop_name = "dt Multiplier"
    else:
        print(f"Error: Invalid STUDY_MODE '{STUDY_MODE}'.")
        sys.exit(1)


    print("=" * 60)
    print(f"{'CONVERGENCE STUDY PIPELINE':^60}")
    print("=" * 60)
    print(f"STUDY MODE: {STUDY_MODE}")
    print(f"Target T_END: {T_END}")
    print(f"Parameters to test ({loop_name}): {loop_values}")
    if STUDY_MODE == "TEMPORAL_ONLY":
        print(f"Fixed Spatial Grid: Nx = {FIXED_NX_FOR_TEMPORAL_STUDY}")
    print("-" * 60)

    try:
        # --- PHASE 1: RUN SIMULATIONS ---
        for i, loop_val in enumerate(loop_values):
            
            # --- CALCULATE PARAMETERS ---
            nx, dt_actual, num_steps = calculate_params_for_mode(loop_val, BASE_DT, STUDY_MODE)

            if STUDY_MODE == "TEMPORAL_ONLY":
                label_val = f"dt_mult={loop_val} (dt={dt_actual:.6e})"
            else:
                label_val = f"Nx={nx} (dt={dt_actual:.6e})"
            
            print(f"\n[Run {i+1}/{len(loop_values)}] {label_val} | Steps={num_steps}...")

            # A. Create a DEEP COPY of the original config
            current_config = copy.deepcopy(original_config)

            # B. Apply modifications
            current_config['mesh']['nx'] = nx
            current_config['mesh']['input_for_manufactured_solution'] = True
            
            current_config['time']['dt'] = dt_actual
            current_config['time']['t_end'] = T_END

            # 4. Modify Output settings
            # We only save the LAST step to save space and simplify analysis
            if STUDY_MODE == "TEMPORAL_ONLY":
                 # >>> NUOVA LOGICA PER NOME FILE SEMPLICE <<<
                 dt_identifier = int(loop_val * 10000) 
                 # Es. 1.0 -> 10000, 0.5 -> 05000
                 base_name = f"temporal_n{nx}_mult_{dt_identifier:05d}_output"
            else:
                 # Use nx for filename in spatial studies
                 base_name = f"spatial_n{nx}_output"
                 
            current_config['output']['base_filename'] = base_name
            current_config['output']['output_frequency'] = num_steps # Save only at the end

            # C. Overwrite config.json on disk
            with open(CONFIG_PATH, 'w') as f:
                json.dump(current_config, f, indent=4)

            # D. Launch C++ Solver
            # ret = subprocess.run([EXECUTABLE_PATH], capture_output=True, text=True)
            ret = subprocess.run(["mpirun", "-np", "1", EXECUTABLE_PATH], capture_output=True, text=True)

            if ret.returncode != 0:
                print(f"  [ERROR] C++ execution failed for {label_val}:")
                print(f"  {'-'*20} STDERR {'-'*20}")
                print(ret.stderr)
                print(f"  {'-'*48}")
                continue
            else:
                print("  > Simulation completed successfully.")

            # E. Verify output existence
            expected_filename = f"{base_name}_0000_{num_steps:04d}.vtk"
            full_path = os.path.join(OUTPUT_DIR, expected_filename)

            if os.path.exists(full_path):
                print(f"  > Verified output: {expected_filename}")
                simulation_data_for_plot.append({
                    "nx": nx, 
                    "file": full_path,
                    "dt": dt_actual, 
                    "mode": STUDY_MODE
                })
            else:
                print(f"  [WARNING] Expected output file missing: {full_path}")

    except KeyboardInterrupt:
        print("\n\n[!] Process interrupted by user.")

    finally:
        # --- RESTORE PHASE (Always executed) ---
        print("\n" + "-" * 60)
        print("Restoring original configuration...")
        with open(CONFIG_PATH, 'w') as f:
            f.write(raw_original_content)
        print("Configuration restored.")

    # --- PHASE 2: ANALYSIS & PLOTTING ---
    print("\n" + "=" * 60)
    print(f"{'POST-PROCESSING & ANALYSIS':^60}")
    print("=" * 60)

    if len(simulation_data_for_plot) < 2:
        print("[ERROR] Not enough valid data points (min 2) for convergence plot.")
    else:
        print(f"Calling analysis script in {STUDY_MODE} mode...")
        plot_convergence_new.run_analysis(custom_simulations=simulation_data_for_plot)

if __name__ == "__main__":
    main()