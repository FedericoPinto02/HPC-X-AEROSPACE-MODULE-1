#!/usr/bin/env python3
import os
import re
import json
import sys
import csv
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

# =============================================================================
# === USER CONFIGURATION ===
# =============================================================================

# Default Physical parameters
RE_NUMBER = 1.0

# Paths
DATA_DIR = "../output/"
CONFIG_FILE = "../data/config.json"

# File prefix
FILE_PREFIX = "findTot7"

# Fields
FIELD_P = "pressure"
FIELD_DIV = "divU"

# Probe configuration
NUM_PROBES_TOTAL = 1000  
NUM_PROBES_TO_PLOT = 50
RANDOM_SEED = 42

# Percentuale di probe da forzare sui bordi "destri" (Max X, Max Y, Max Z)
PROBE_BOUNDARY_RATIO = 0.2 

DEFAULT_DT = 0.07

# SAMPLING FREQUENCY
# 1 = Read every file. 2 = Read every 2nd file, etc.
VTK_SAMPLE_STEP = 4

# =============================================================================
# === UTILITIES ===
# =============================================================================

def get_step_from_filename(filename):
    match = re.search(r'_(\d+)\.vtk$', filename)
    return int(match.group(1)) if match else -1

def get_config_json():
    try:
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as f:
                return json.load(f)
    except Exception:
        pass
    return {}

def get_dt_from_config():
    conf = get_config_json()
    try:
        return float(conf.get('time', {}).get('dt', DEFAULT_DT))
    except Exception:
        return DEFAULT_DT

def get_re_from_config(default_re):
    conf = get_config_json()
    # Check simple keys
    for key in ('Re', 're', 'RE', 're_number'):
        if conf.get(key) is not None: return float(conf[key])
    # Check nested physics
    if isinstance(conf.get('physics'), dict):
        for key in ('Re', 're'):
            if conf['physics'].get(key) is not None: return float(conf['physics'][key])
    return default_re

# --- ANALYTICAL FUNCTIONS ---

def p_analytical(t_array, x, y, z, re_number):
    t_arr = np.asarray(t_array)
    amplitude = 3.0 / re_number
    # Formula: p = A * cos(x) * cos(t + y) * cos(z)
    return amplitude * np.cos(x) * np.cos(y + t_arr) * np.cos(z)

def vel_analytical(t, x, y, z):
    u = np.sin(x) * np.cos(t + y) * np.sin(z)
    v = np.cos(x) * np.sin(t + y) * np.sin(z)
    w = 2 * np.cos(x) * np.cos(t + y) * np.cos(z)
    return u, v, w

def vel_spatial_coeffs(x, z):
    u_c = np.sin(x) * np.sin(z)
    v_c = np.cos(x) * np.sin(z)
    w_c = 2 * np.cos(x) * np.cos(z)
    return u_c, v_c, w_c

# =============================================================================
# === MAIN ===
# =============================================================================

def main():
    global RE_NUMBER
    print(f"--- PIT POINT ANALYSIS (Divergence Trends & T0 Exclusion) ---")

    RE_NUMBER = get_re_from_config(RE_NUMBER)
    print(f"Using Re = {RE_NUMBER}")

    if not os.path.exists(DATA_DIR):
        print(f"ERROR: data directory {DATA_DIR} not found.")
        sys.exit(1)

    all_files = os.listdir(DATA_DIR)
    vtk_files = [os.path.join(DATA_DIR, f) for f in all_files if f.startswith(FILE_PREFIX) and f.endswith(".vtk")]
    if not vtk_files:
        print(f"No VTK files found with prefix '{FILE_PREFIX}'")
        sys.exit(1)

    vtk_files.sort(key=lambda p: get_step_from_filename(os.path.basename(p)))
    
    # --- APPLY SAMPLING FREQUENCY ---
    if VTK_SAMPLE_STEP > 1:
        print(f"-> Subsampling VTK files: keeping 1 every {VTK_SAMPLE_STEP} files.")
        vtk_files = vtk_files[::VTK_SAMPLE_STEP]
    
    dt = get_dt_from_config()

    # --- READ FIRST FILE & DETECT DATA & SETUP DIRS ---
    try:
        first_file = vtk_files[0]
        first_grid = pv.read(first_file)
        step0 = get_step_from_filename(os.path.basename(first_file))
        
        dims = first_grid.dimensions
        nx_val = dims[0]
        
    except Exception as e:
        print(f"Error reading first mesh: {e}")
        sys.exit(1)

    # --- DYNAMIC OUTPUT DIRECTORY ---
    output_folder_name = f"probe_analysis_Nx{nx_val}_dt{dt}"
    full_output_dir = os.path.join(DATA_DIR, "../scripts", output_folder_name)
    os.makedirs(full_output_dir, exist_ok=True)
    print(f"-> Output Directory: {full_output_dir}")

    # Auto-detect logic
    data_source_type = "UNKNOWN"
    available_points = None
    
    if FIELD_P in first_grid.point_data:
        data_source_type = "POINT_DATA"
        print(f"-> DETECTED: '{FIELD_P}' is in POINT_DATA. Using grid nodes.")
        available_points = first_grid.points
    elif FIELD_P in first_grid.cell_data:
        data_source_type = "CELL_DATA"
        print(f"-> DETECTED: '{FIELD_P}' is in CELL_DATA. Using cell centers.")
        available_points = first_grid.cell_centers().points
    else:
        print(f"ERROR: Field '{FIELD_P}' not found.")
        sys.exit(1)

    n_locations = len(available_points)
    
    # --- SELEZIONE PROBE (MISTA INTERNA + BORDI) ---
    rng = np.random.default_rng(RANDOM_SEED)
    
    if n_locations > NUM_PROBES_TOTAL:
        # Trova i limiti del dominio
        max_x = np.max(available_points[:, 0])
        max_y = np.max(available_points[:, 1])
        max_z = np.max(available_points[:, 2])
        
        tol = 1e-4 # Tolleranza per identificare i bordi
        
        # Indici sui bordi "destri" (coordinate massime)
        on_max_x = np.where(np.abs(available_points[:, 0] - max_x) < tol)[0]
        on_max_y = np.where(np.abs(available_points[:, 1] - max_y) < tol)[0]
        on_max_z = np.where(np.abs(available_points[:, 2] - max_z) < tol)[0]
        
        # Unione di tutti gli indici di bordo (senza duplicati)
        boundary_candidates = np.unique(np.concatenate((on_max_x, on_max_y, on_max_z)))
        
        # Indici interni (o non sui massimi)
        internal_candidates = np.setdiff1d(np.arange(n_locations), boundary_candidates)
        
        # Calcolo quanti prenderne
        n_boundary_target = int(NUM_PROBES_TOTAL * PROBE_BOUNDARY_RATIO)
        n_internal_target = NUM_PROBES_TOTAL - n_boundary_target
        
        # Selezione casuale dai gruppi
        n_boundary_actual = min(len(boundary_candidates), n_boundary_target)
        selected_boundary = rng.choice(boundary_candidates, size=n_boundary_actual, replace=False)
        
        n_internal_actual = NUM_PROBES_TOTAL - n_boundary_actual
        if n_internal_actual > len(internal_candidates):
             selected_internal = internal_candidates
        else:
            selected_internal = rng.choice(internal_candidates, size=n_internal_actual, replace=False)
            
        print(f"   Selecting Probes: {len(selected_boundary)} on Max Boundaries, {len(selected_internal)} Internal/Other.")
        
        probe_indices = np.concatenate((selected_boundary, selected_internal))
        probe_indices.sort()
    else:
        print(f"   Using all {n_locations} locations.")
        probe_indices = np.arange(n_locations)
    
    probe_points = available_points[probe_indices]
    num_probes = len(probe_indices)
    
    np.savetxt(os.path.join(full_output_dir, "probe_coordinates.txt"), probe_points, header="X,Y,Z", delimiter=",")

    # --- SAMPLING LOOP ---
    probe_histories = [{'p': [], 'div': []} for _ in range(num_probes)]
    times = []

    print(f"-> Sampling data ({len(vtk_files)} steps)...", end='', flush=True)

    for idx, f_path in enumerate(vtk_files):
        step = get_step_from_filename(os.path.basename(f_path))
        t = (step - step0) * dt
        times.append(t)

        try:
            grid = pv.read(f_path)
            
            # Select correct data dictionary
            current_data = grid.point_data if data_source_type == "POINT_DATA" else grid.cell_data

            if FIELD_P in current_data:
                p_vals = current_data[FIELD_P][probe_indices]
            else:
                p_vals = np.zeros(num_probes)

            # DivU logic
            div_key = None
            if FIELD_DIV in current_data: div_key = FIELD_DIV
            elif 'divergence' in current_data: div_key = 'divergence'
            
            div_vals = current_data[div_key][probe_indices] if div_key else np.zeros(num_probes)
            
            for i in range(num_probes):
                probe_histories[i]['p'].append(p_vals[i])
                probe_histories[i]['div'].append(div_vals[i])

        except Exception:
            for i in range(num_probes):
                probe_histories[i]['p'].append(0.0)
                probe_histories[i]['div'].append(0.0)

        if idx % 10 == 0: print(".", end='', flush=True)

    print("\n-> Sampling complete.")
    times = np.array(times)

    # --- ANALYSIS & PLOTTING ---
    print("-> Performing statistical analysis (skipping first step for stats)...")
    summary_data = []

    # Collectors for global summary plots (stats only)
    all_l2_errors_p = []
    all_max_errors_p = []
    all_l2_div = []
    all_max_div = []
    all_mean_div = []
    all_std_div = []

    for i in range(num_probes):
        # Raw sequences (including t=0)
        p_seq = np.array(probe_histories[i]['p'])
        div_seq = np.array(probe_histories[i]['div'])
        x_p, y_p, z_p = probe_points[i]

        # Analytical Reference (full time)
        p_ana_seq = p_analytical(times, x_p, y_p, z_p, RE_NUMBER)
        
        # Calculate Error Signal (full time)
        error_seq = p_seq - p_ana_seq

        # --- STATS CALCULATION (EXCLUDING T=0) ---
        if len(p_seq) > 1:
            valid_slice = slice(1, None) # Skip index 0
        else:
            valid_slice = slice(0, None) # Fallback if only 1 step
            
        p_seq_stat = p_seq[valid_slice]
        p_ana_stat = p_ana_seq[valid_slice]
        error_seq_stat = error_seq[valid_slice]
        div_seq_stat = div_seq[valid_slice]
        
        # Pressure Stats
        l2_error_p = np.sqrt(np.mean(error_seq_stat**2))
        max_abs_error_p = np.max(np.abs(error_seq_stat))
        
        # Divergence Stats (Advanced)
        # Analytical Div should be 0, so Div value is the error itself
        l2_div = np.sqrt(np.mean(div_seq_stat**2))
        max_div = np.max(np.abs(div_seq_stat))
        mean_div = np.mean(div_seq_stat) # Arithmetic Mean (Bias)
        std_div = np.std(div_seq_stat)   # Standard Deviation (Fluctuation)

        # Simple Amplitude (Full sequence usually better for peaks, but consistency matters)
        # Using full sequence for amplitude to capture peaks if simulation is short
        amp_sim = (np.max(p_seq) - np.min(p_seq)) / 2.0
        amp_ana = (np.max(p_ana_seq) - np.min(p_ana_seq)) / 2.0

        p_sim_t0 = p_seq[0] if len(p_seq) > 0 else 0.0

        all_l2_errors_p.append(l2_error_p)
        all_max_errors_p.append(max_abs_error_p)
        all_l2_div.append(l2_div)
        all_max_div.append(max_div)
        all_mean_div.append(mean_div)
        all_std_div.append(std_div)

        summary_data.append({
            'ID': i + 1,
            'Index': probe_indices[i],
            'X': x_p, 'Y': y_p, 'Z': z_p,
            'P_Sim_t0': p_sim_t0,
            'L2_Error_P': l2_error_p,
            'Max_Abs_Error_P': max_abs_error_p,
            'L2_DivU': l2_div,
            'Max_DivU': max_div,
            'Mean_DivU': mean_div,
            'Std_DivU': std_div,
            'Amp_Simple_Sim': amp_sim,
            'Amp_Simple_Ana': amp_ana
        })

        # --- INDIVIDUAL PLOTS ---
        if i < NUM_PROBES_TO_PLOT:
            coord_str = f"({x_p:.2f}, {y_p:.2f}, {z_p:.2f})"
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
            
            fig.suptitle(f'Probe #{i+1} {coord_str}\n'
                         f'ErrP L2: {l2_error_p:.2e} | Div L2: {l2_div:.2e} (Stats ex. t0)', 
                         color='black')

            # Plot 1: Pressures
            ax1.plot(times, p_seq, 'b-', label='Sim', alpha=0.8, linewidth=1.5)
            ax1.plot(times, p_ana_seq, 'g--', label='Analytical', linewidth=1.5)
            ax1.set_ylabel('Pressure')
            ax1.legend(loc='upper right')
            ax1.grid(True, alpha=0.3)

            # Plot 2: Pressure Error
            ax2.plot(times, error_seq, 'purple', label='Error (Sim-Ana)')
            # Mark the excluded region visually?
            if len(times) > 1:
                ax2.axvspan(times[0], times[1], color='red', alpha=0.1, label='Excluded from Stats')
            ax2.set_ylabel('Press. Error')
            ax2.grid(True, alpha=0.3)
            ax2.legend(loc='upper right', fontsize='small')

            # Plot 3: Divergence
            ax3.plot(times, div_seq, 'orange', label='Div U')
            ax3.axhline(mean_div, color='darkorange', linestyle='--', alpha=0.8, linewidth=1, label=f'Mean (ex t0): {mean_div:.2e}')
            ax3.set_ylabel('Divergence')
            ax3.set_xlabel('Time (s)')
            ax3.grid(True, alpha=0.3)
            ax3.legend(loc='upper right', fontsize='small')
            
            max_t = times[-1]
            num_pis = int(max_t / np.pi) + 1
            for ax in (ax1, ax2, ax3):
                for k in range(1, num_pis + 1):
                    val = k * np.pi
                    if val <= max_t:
                        ax.axvline(val, color='gray', linestyle='-', alpha=0.2)

            plt.tight_layout()
            plt.savefig(os.path.join(full_output_dir, f"probe_{i+1}_trace.png"), dpi=100)
            plt.close()

    # --- CSV WRITING ---
    csv_path = os.path.join(full_output_dir, "analysis_summary.csv")
    keys = summary_data[0].keys()
    with open(csv_path, 'w', newline='') as f:
        dict_writer = csv.DictWriter(f, keys)
        dict_writer.writeheader()
        dict_writer.writerows(summary_data)

    # --- GLOBAL SUMMARY PLOT (3 PANELS) ---
    print("-> Generating Global Summary Plot...")
    
    fig_sum, (ax_p_hist, ax_div_hist, ax_scat) = plt.subplots(1, 3, figsize=(18, 5))
    
    # 1. Histogram of Pressure L2 Errors
    ax_p_hist.hist(all_l2_errors_p, bins=20, color='skyblue', edgecolor='black')
    ax_p_hist.set_title(f"Pressure L2 Error Dist.")
    ax_p_hist.set_xlabel("L2 Error (Pa)")
    ax_p_hist.set_ylabel("Count")
    ax_p_hist.grid(True, linestyle=':', alpha=0.5)

    # 2. Histogram of Divergence L2
    ax_div_hist.hist(all_l2_div, bins=20, color='orange', edgecolor='black')
    ax_div_hist.set_title(f"Divergence L2 Dist.")
    ax_div_hist.set_xlabel("L2 Divergence (1/s)")
    ax_div_hist.grid(True, linestyle=':', alpha=0.5)

    # 3. Scatter Plot: Pressure Error vs ID (Color = Mean Divergence Bias)
    ids = np.arange(1, num_probes + 1)
    # Using abs(Mean Div) for color to see bias magnitude
    sc = ax_scat.scatter(ids, all_l2_errors_p, c=np.abs(all_mean_div), cmap='plasma', s=20)
    ax_scat.set_title("P Error vs ID (Color = Abs Mean DivU)")
    ax_scat.set_xlabel("Probe ID")
    ax_scat.set_ylabel("Pressure L2 Error")
    ax_scat.grid(True, linestyle=':', alpha=0.5)
    plt.colorbar(sc, ax=ax_scat, label="Abs Mean DivU")
    
    sum_plot_path = os.path.join(full_output_dir, "global_summary.png")
    plt.tight_layout()
    plt.savefig(sum_plot_path, dpi=150)
    plt.close()

    # --- TEXT REPORT ---
    report_path = os.path.join(full_output_dir, "summary_report.txt")
    
    max_p_err_idx = np.argmax(all_l2_errors_p)
    max_div_idx = np.argmax(all_l2_div)

    with open(report_path, "w") as f:
        f.write("=== PROBE ANALYSIS REPORT ===\n")
        f.write(f"Total Probes Analyzed: {num_probes}\n")
        f.write(f"  - Boundary Probes (Max X/Y/Z): Included\n")
        f.write(f"Files Processed: {len(vtk_files)} (Sample Step: {VTK_SAMPLE_STEP})\n")
        f.write(f"Statistics Exclude First Step (t=0): YES\n\n")
        
        f.write("--- PRESSURE Statistics (L2 Error) ---\n")
        f.write(f"Average: {np.mean(all_l2_errors_p):.6e}\n")
        f.write(f"Median:  {np.median(all_l2_errors_p):.6e}\n")
        f.write(f"Std Dev: {np.std(all_l2_errors_p):.6e}\n")
        f.write(f"Worst Probe ID: {max_p_err_idx + 1} (Val: {all_l2_errors_p[max_p_err_idx]:.6e})\n")
        f.write(f"  -> Location: {probe_points[max_p_err_idx]}\n\n")

        f.write("--- DIVERGENCE Statistics (L2 Norm) ---\n")
        f.write(f"Average: {np.mean(all_l2_div):.6e}\n")
        f.write(f"Median:  {np.median(all_l2_div):.6e}\n")
        f.write(f"Worst Probe ID: {max_div_idx + 1} (Val: {all_l2_div[max_div_idx]:.6e})\n")
        f.write(f"  -> Location: {probe_points[max_div_idx]}\n")
        f.write(f"  -> Mean Div (Bias): {all_mean_div[max_div_idx]:.6e}\n")
        f.write(f"  -> Std Div (Noise): {all_std_div[max_div_idx]:.6e}\n")

    print(f"-> Analysis Done.")
    print(f"   CSV: {csv_path}")
    print(f"   Summary Plot: {sum_plot_path}")
    print(f"   Text Report: {report_path}")

if __name__ == "__main__":
    main()