import json
import subprocess
import os
import sys
import copy  # Necessary for deep copying nested dictionaries

# Attempt to import the plotting script
try:
    import plot_convergence
except ImportError:
    print("Error: Could not find 'plot_convergence.py' in the same directory.")
    sys.exit(1)

# --- USER CONFIGURATION ---
EXECUTABLE_PATH = "../build/main"      # Path to your C++ executable
CONFIG_PATH = "../data/config.json"    # Path to the JSON config file
OUTPUT_DIR = "../output/"              # Output directory (must match JSON)

# Time parameters to calculate expected filename suffix (e.g., _0030.vtk)
# Ensure these match the JSON logic.
EXPECTED_LAST_STEP = 30  # t_end / dt = 0.03 / 0.001 = 30

# Resolutions to test (Nx)
NX_VALUES = [20, 25, 30] 

def main():
    # 1. Read the original config
    # We read it as raw text to preserve formatting/comments for restoration
    try:
        with open(CONFIG_PATH, 'r') as f:
            raw_original_content = f.read()
        
        # Parse JSON from string to allow modifications
        original_config = json.loads(raw_original_content)
            
    except FileNotFoundError:
        print(f"Error: Configuration file not found at {CONFIG_PATH}")
        sys.exit(1)

    simulation_data_for_plot = []

    print("=" * 60)
    print(f"{'CONVERGENCE STUDY PIPELINE':^60}")
    print("=" * 60)
    print(f"Original configuration backed up.")
    print(f"Target Resolutions (Nx): {NX_VALUES}")
    print("-" * 60)

    try:
        # --- PHASE 1: RUN SIMULATIONS ---
        for i, nx in enumerate(NX_VALUES):
            print(f"\n[Run {i+1}/{len(NX_VALUES)}] Processing Nx = {nx}...")

            # A. Create a DEEP COPY of the original config
            current_config = copy.deepcopy(original_config)

            # B. Apply modifications
            # 1. Set Nx
            current_config['mesh']['nx'] = nx
            # 2. Force manufactured solution flag
            current_config['mesh']['input_for_manufactured_solution'] = True
            # 3. Modify base filename
            # This lets C++ save directly as "n80_simulation_output_0030.vtk"
            base_name = f"n{nx}_simulation_output"
            current_config['output']['base_filename'] = base_name

            # C. Overwrite config.json on disk
            with open(CONFIG_PATH, 'w') as f:
                json.dump(current_config, f, indent=4)

            # D. Launch C++ Solver
            # print(f"  > Executing solver...") # Optional verbosity
            ret = subprocess.run([EXECUTABLE_PATH], capture_output=True, text=True)
            
            if ret.returncode != 0:
                print(f"  [ERROR] C++ execution failed for Nx={nx}:")
                print(f"  {'-'*20} STDERR {'-'*20}")
                print(ret.stderr)
                print(f"  {'-'*48}")
                continue # Skip to next resolution
            else:
                print("  > Simulation completed successfully.")

            # E. Verify output existence
            expected_filename = f"{base_name}_{EXPECTED_LAST_STEP:04d}.vtk"
            full_path = os.path.join(OUTPUT_DIR, expected_filename)

            if os.path.exists(full_path):
                print(f"  > Verified output: {expected_filename}")
                simulation_data_for_plot.append({
                    "nx": nx,
                    "file": full_path
                })
            else:
                print(f"  [WARNING] Expected output file missing: {full_path}")
                print("  > Check 'output_frequency' in JSON vs 't_end'.")

    except KeyboardInterrupt:
        print("\n\n[!] Process interrupted by user.")

    finally:
        # --- RESTORE PHASE (Always executed) ---
        print("\n" + "-" * 60)
        print("Restoring original configuration...")
        with open(CONFIG_PATH, 'w') as f:
            f.write(raw_original_content) # Write back exact original string
        print("Configuration restored.")

    # --- PHASE 2: ANALYSIS & PLOTTING ---
    print("\n" + "=" * 60)
    print(f"{'POST-PROCESSING & ANALYSIS':^60}")
    print("=" * 60)

    if len(simulation_data_for_plot) < 2:
        print("[ERROR] Not enough valid data points (min 2) for convergence plot.")
    else:
        # Call the plotting script
        plot_convergence.run_analysis(custom_simulations=simulation_data_for_plot)

if __name__ == "__main__":
    main()