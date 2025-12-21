import subprocess
import numpy as np
import pyvista as pv
import os
import sys
import json
import glob
from scipy.spatial import cKDTree

# =============================================================================
# CONFIGURATION
# =============================================================================
EXECUTABLE = "../build/main"
CONFIG_PATH = "../data/config.json"
OUTPUT_DIR = "../output"
TOLERANCE  = 1e-5
GEOMETRY_TOLERANCE = 1e-5

# FORCE DIMENSIONS HERE to avoid config.json ambiguity
NX = 32
NY = 32
NZ = 32
# =============================================================================

def clean_directory(base_name):
    for f in glob.glob(os.path.join(OUTPUT_DIR, f"{base_name}*.vtk")):
        try: os.remove(f)
        except: pass

def run_simulation(n_procs, label):
    print(f"\n[{label.upper()}] Running with NP={n_procs} on Grid {NX}x{NY}x{NZ}...")
    base_name = f"verify_{label}"

    with open(CONFIG_PATH, 'r') as f: config = json.load(f)

    # 1. ENFORCE GRID DIMENSIONS
    if 'grid' not in config: config['grid'] = {}
    config['mesh']['nx'] = NX
    config['mesh']['ny'] = NY
    config['mesh']['nz'] = NZ

    # 2. ENFORCE I/O
    config['output']['base_filename'] = base_name
    config['output']['output_frequency'] = 1
    config['time']['dt'] = 0.01
    config['time']['t_end'] = 0.1

    with open(CONFIG_PATH, 'w') as f: json.dump(config, f, indent=4)

    clean_directory(base_name)

    cmd = ["mpirun", "-np", str(n_procs), EXECUTABLE]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)

    file_pattern = f"{base_name}_*_0001.vtk"
    files = sorted(glob.glob(os.path.join(OUTPUT_DIR, file_pattern)))

    if not files: raise FileNotFoundError(f"No output files for {label}!")

    # Combine and Clean
    # Note: clean() tolerance might need adjustment if dx changes significantly
    mesh = pv.read(files).combine().clean(tolerance=GEOMETRY_TOLERANCE)
    return mesh

def main():
    try:
        mesh_serial = run_simulation(1, "serial")
        mesh_parallel = run_simulation(4, "parallel") # Change to 4 to test splitting

        print("\n[GEOMETRY CHECK]")
        print(f"  Serial Points:   {mesh_serial.n_points}")
        print(f"  Parallel Points: {mesh_parallel.n_points}")

        if mesh_serial.n_points != mesh_parallel.n_points:
            print("  [WARNING] Point counts differ. Possible grid decomposition mismatch.")

        print("\n[DATA VALIDATION]")
        tree = cKDTree(mesh_serial.points)
        dist, indices = tree.query(mesh_parallel.points)

        if np.max(dist) > 1e-5:
            print(f"  [FATAL] Geometric mismatch. Max dist: {np.max(dist):.4e}")
            sys.exit(1)

        passed = True
        for field in ["pressure", "velocity"]:
            if field not in mesh_serial.point_data: continue

            data_serial = mesh_serial.point_data[field][indices]
            data_parallel = mesh_parallel.point_data[field]

            diff = np.abs(data_serial - data_parallel)
            if diff.ndim > 1: diff = np.max(diff, axis=1)

            max_err = np.max(diff)
            status = "PASS" if max_err < TOLERANCE else "FAIL"
            print(f"  Field {field:10}: Max Error = {max_err:.4e} [{status}]")
            if status == "FAIL": passed = False

        if passed:
            print("SUCCESS: Parallel solver is consistent.")
            sys.exit(0)
        else:
            print("FAILURE: Inconsistencies detected.")
            sys.exit(1)

    except Exception as e:
        print(f"\n[ERROR] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()