#!/usr/bin/env python3

"""
This script performs a spatial convergence analysis for a Method of
Manufactured Solutions (MMS) test.

It reads multiple VTK files, each from a simulation with a different
grid resolution. For each file, it computes the L2 RMS error
of velocity components (u, v, w) and pressure (p) against a
known analytical solution.

This version is adapted for staggered grid solvers where the VTK output
is co-located (e.g., at cell centers). It compares the numerical
co-located data against the analytical solution evaluated at the
appropriate staggered locations (e.g., cell faces).

Script modified from a previous version to:
1. Use only L2 RMS (absolute) errors.
2. Add O(h) and O(h^2) reference lines.
3. Translate all comments and outputs to English.
4. Remove all divergence calculations.
5. Save plots to high-definition image files instead of showing them.
"""

# --- 0. IMPORTS ---
#!pip install pyvista numpy matplotlib
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re  # Import regular expressions for parsing

# --- 1. USER CONFIGURATION ---
# IMPORTANT: Update these values for your specific simulation

# List of simulation results to analyze.
# Add or remove dictionaries as needed.
# "nx": Number of cells in the x-direction (used to calculate h)
# "file": Path to the VTK output file.
#       *** The filename MUST contain the step number (e.g., _0030.vtk) ***
SIMULATIONS = [
{"nx": 30, "file": "../output/n30simulation_output_0030.vtk"},
{"nx": 20, "file": "../output/n20simulation_output_0030.vtk"},
# {"nx": 40, "file": "../output/n40simulation_output_0030.vtk"},
# {"nx": 45, "file": "../output/n45simulation_output_0030.vtk"},
# {"nx": 35, "file": "../output/n35simulation_output_0030.vtk"},
# {"nx": 60, "file": "../output/n60simulation_output_0030.vtk"},
# {"nx": 50, "file": "../output/n50simulation_output_0030.vtk"},
# {"nx": 80, "file": "../output/n80simulation_output_0030.vtk"},
]

# Physical domain length in the x-direction (e.g., L_x)
DOMAIN_LENGTH_X = 6

# Add the simulation timestep (dt)
# This is VITAL for calculating the correct time.
DT = 0.001

# Field names in the VTK file
VELOCITY_FIELD_NAME = "velocity"
PRESSURE_FIELD_NAME = "pressure"


# Helper function to get step number from filename
def get_step_from_filename(filepath: str) -> int:
    """
    Parses a filename (e.g., "sim_0319.vtk") and extracts
    the step number (319).
    """
    # This regex looks for one or more digits (\d+) right before ".vtk"
    # It assumes the format is like "_123.vtk" or "abc_123.vtk"
    match = re.search(r'_(\d+)\.vtk$', os.path.basename(filepath))

    if match:
        return int(match.group(1))
    else:
        # If no match, raise an error
        raise ValueError(f"Could not parse step number from filename: {filepath}. "
                         "Expected format like '..._123.vtk'")


# --- 2. ANALYTICAL SOLUTION FUNCTIONS (STAGGERED) ---
# (These functions assume correct staggered evaluation)

def compute_analytical_u_at_x_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    """
    Computes analytical u-velocity.
    Assumes 'points' are cell centers (i, j, k).
    Evaluates 'u' at the x-face center (i+0.5, j, k)
    by shifting x coordinates by +h/2.
    """
    x = points[:, 0] + h / 2.0  # Staggered shift
    y = points[:, 1]
    z = points[:, 2]

    # u = sin(t) * sin(x) * sin(y) * sin(z)
    u_ana = np.sin(t) * np.sin(x) * np.sin(y) * np.sin(z)
    return u_ana


def compute_analytical_v_at_y_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    """
    Computes analytical v-velocity.
    Assumes 'points' are cell centers (i, j, k).
    Evaluates 'v' at the y-face center (i, j+0.5, k)
    by shifting y coordinates by +h/2.
    """
    x = points[:, 0]
    y = points[:, 1] + h / 2.0  # Staggered shift
    z = points[:, 2]

    # v = sin(t) * cos(x) * cos(y) * cos(z)
    v_ana = np.sin(t) * np.cos(x) * np.cos(y) * np.cos(z)
    return v_ana


def compute_analytical_w_at_z_faces(points: np.ndarray, t: float, h: float) -> np.ndarray:
    """
    Computes analytical w-velocity.
    Assumes 'points' are cell centers (i, j, k).
    Evaluates 'w' at the z-face center (i, j, k+0.5)
    by shifting z coordinates by +h/2.
    """
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2] + h / 2.0  # Staggered shift

    # w = sin(t) * cos(x) * sin(y) * (sin(z) + cos(z))
    w_ana = np.sin(t) * np.cos(x) * np.sin(y) * (np.sin(z) + np.cos(z))
    return w_ana


def compute_analytical_p_at_centers(points: np.ndarray, t: float) -> np.ndarray:
    """
    Computes analytical pressure.
    Assumes 'points' are cell centers (i, j, k) and evaluates p there.
    """
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    # p = - sin(t) * 3.0 * cos(x) * sin(y) * (sin(z) - cos(z))
    p_ana = -np.sin(t) * 3.0 * np.cos(x) * np.sin(y) * (np.sin(z) - np.cos(z))
    return p_ana

# --- 3. ERROR CALCULATION HELPERS ---

def calculate_l2_rms(data: np.ndarray) -> float:
    """
    Calculates the L2 Root Mean Square (RMS) value of a data field.
    This works for both scalar (N,) and vector (N, 3) fields.
    """
    data_sq = data**2

    if data.ndim == 2:
        # For vectors, sum component squares
        data_sq = np.sum(data_sq, axis=1)

    mean_squared_value = np.mean(data_sq)
    return np.sqrt(mean_squared_value + 1e-16)


def calculate_l2_rms_error(field_numerical: np.ndarray, field_analytical: np.ndarray) -> float:
    """
    Calculates the L2 RMS of the error (field_numerical - field_analytical).
    This is now an absolute error (the RMS norm of the error field).
    """
    error_field = field_numerical - field_analytical

    rms_error = calculate_l2_rms(error_field)

    # We no longer divide by the analytical norm
    return rms_error


def plot_convergence(h_values: list, error_data: dict, title: str, save_filename: str):
    h_arr = np.array(h_values)

    if len(h_arr) == 0:
        print("No data to plot.", file=sys.stderr)
        return

    # --- FIX 1: ORDINARE I DATI ---
    # È fondamentale ordinare h dal più piccolo al più grande per evitare linee a zig-zag
    sort_idx = np.argsort(h_arr)
    h_arr = h_arr[sort_idx]

    plottable_error_data = {}
    for label, errors in error_data.items():
        if len(errors) == len(h_arr):
            # Riordiniamo anche gli errori usando gli stessi indici di h
            plottable_error_data[label] = np.array(errors)[sort_idx]
        else:
            print(f"WARNING: Skipping '{label}' (dim mismatch)", file=sys.stderr)

    if not plottable_error_data:
        return

    fig, ax = plt.subplots(figsize=(10, 7))

    markers = ['o', 's', '^', 'D', 'v', '*']

    # Plotting Data
    for i, (label, errors) in enumerate(plottable_error_data.items()):
        ax.plot(h_arr, errors, marker=markers[i % len(markers)], markersize=8, linestyle='-', label=label)

    # --- Reference Lines ---
    if len(h_arr) > 1:
        # Prendiamo il punto più a destra (h più grande) come ancora,
        # perché solitamente è il punto più stabile numericamente per calcolare la pendenza teorica
        h_ref = h_arr[-1]
        # Usiamo l'errore della prima serie associata a quel h
        err_ref = list(plottable_error_data.values())[0][-1]

        h_line = np.array([h_arr.min(), h_arr.max()])

        # O(h^2)
        # C2 = err / h^2
        C_2 = err_ref / (h_ref**2)
        ax.plot(h_line, C_2 * (h_line**2), 'k--', label='O($h^2$) Reference')

        # O(h)
        C_1 = err_ref / (h_ref**1)
        ax.plot(h_line, C_1 * (h_line**1), 'k:', label='O($h$) Reference')

    # --- Formatting ---
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Grid Size $h$ (m)', fontsize=14)
    ax.set_ylabel("L2 RMS Error", fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True, which="both", linestyle=':', linewidth=0.5)
    ax.legend(fontsize=12)

    # --- FIX 2: FORMATTAZIONE ASSI ---
    # Rimuoviamo il FormatStrFormatter rigido.
    # Usiamo un LogFormatter o lasciamo automatico.
    # Se vuoi numeri leggibili ma non scientifici dove possibile:
    from matplotlib.ticker import ScalarFormatter
    formatter = ScalarFormatter()
    formatter.set_scientific(False) # Prova a non usare 1e-2 se c'è spazio

    # Oppure, molto meglio per log-log: NON forzare nulla, Matplotlib gestisce bene i log da solo.
    # Se proprio vuoi vedere le potenze:
    # ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())

    plt.tight_layout()

    print(f"Saving plot to {save_filename}...")
    plt.savefig(save_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
# --- 5. MAIN ANALYSIS LOOP ---

def run_analysis():
    """
    Main function to loop through simulations, compute errors, and plot.
    """

    # Lists to store results
    h_values = []
    errors_u = []
    errors_v = []
    errors_w = []
    errors_p = []
    # errors_div list removed

    print("--- Starting Spatial Convergence Analysis (MMS Staggered) ---")
    print(f"Analyzing {len(SIMULATIONS)} files. Using DT={DT} to find time from filenames.")
    # Print table header (Divergence column removed)
    print(f"{'Nx':>4} | {'h':>10} | {'Time':>10} | {'L2_Err_u':>10} | {'L2_Err_v':>10} | {'L2_Err_w':>10} | {'L2_Err_p':>10}")
    print("-" * 74)

    # Ensure simulations are sorted by grid size (nx)
    sorted_simulations = sorted(SIMULATIONS, key=lambda s: s['nx'])

    for sim in sorted_simulations:
        nx = sim['nx']
        filepath = sim['file']

        # Calculate CURRENT_TIME for this specific file
        try:
            step_number = get_step_from_filename(filepath)
            current_time = step_number * DT

        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            print("Skipping this data point.", file=sys.stderr)
            continue

        # Calculate grid size
        # h = L / (N + 0.5)
        h = DOMAIN_LENGTH_X / (nx + 0.5)
        h_values.append(h)

        # --- Load Data ---
        try:
            grid = pv.read(filepath)
        except FileNotFoundError:
            print(f"ERROR: File not found: {filepath}", file=sys.stderr)
            print("Skipping this data point.", file=sys.stderr)
            continue
        except Exception as e:
            print(f"ERROR: Could not read {filepath}: {e}", file=sys.stderr)
            continue

        # These are the co-located points (e.g., cell centers)
        points = grid.points

        # --- 1. Velocity Error (u, v, w) ---
        try:
            # Numerical velocity [u,v,w] read from the co-located points
            vel_numerical = grid.point_data[VELOCITY_FIELD_NAME]

            # Analytical velocity components evaluated at their
            # proper staggered face locations.
            u_analytical = compute_analytical_u_at_x_faces(points, current_time, h)
            v_analytical = compute_analytical_v_at_y_faces(points, current_time, h)
            w_analytical = compute_analytical_w_at_z_faces(points, current_time, h)

            # Calculate absolute L2 RMS error for each component
            err_u = calculate_l2_rms_error(vel_numerical[:, 0], u_analytical)
            err_v = calculate_l2_rms_error(vel_numerical[:, 1], v_analytical)
            err_w = calculate_l2_rms_error(vel_numerical[:, 2], w_analytical)

            errors_u.append(err_u)
            errors_v.append(err_v)
            errors_w.append(err_w)

        except KeyError:
            print(f"ERROR: Field '{VELOCITY_FIELD_NAME}' not found in {filepath}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"ERROR: Velocity analysis failed for {filepath}: {e}", file=sys.stderr)
            continue

        # --- 2. Pressure Error (p) ---
        try:
            # Numerical pressure read from the co-located points
            p_numerical = grid.point_data[PRESSURE_FIELD_NAME]

            # Analytical pressure evaluated at the cell centers
            p_analytical = compute_analytical_p_at_centers(points, current_time)

            # Calculate absolute L2 RMS error for pressure
            err_p = calculate_l2_rms_error(p_numerical, p_analytical)
            errors_p.append(err_p)

        except KeyError:
            print(f"ERROR: Field '{PRESSURE_FIELD_NAME}' not found in {filepath}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"ERROR: Pressure analysis failed for {filepath}: {e}", file=sys.stderr)
            continue

        # --- 3. Divergence Error (div(u)) ---
        # (This entire section has been removed)


        # --- Print Results to Terminal ---
        # (Divergence value removed from the print statement)
        print(f"{nx:>4} | {h:>10.4e} | {current_time:>10.4f} | {err_u:>10.4e} | {err_v:>10.4e} | {err_w:>10.4e} | {err_p:>10.4e}")

    print("-" * 74)
    print("Analysis complete. Saving plots...\n")

    # --- 6. PLOT RESULTS ---

    if not h_values:
        print("No data was successfully processed. Exiting.", file=sys.stderr)
        return

    # Plot 1: Velocity
    vel_error_dict = {
        'u-velocity': errors_u,
        'v-velocity': errors_v,
        'w-velocity': errors_w
    }
    # ### MODIFICATION: Added save_filename ###
    plot_convergence(h_values,
                     vel_error_dict,
                     'Convergence of Velocity Components (Staggered)',
                     save_filename="../results/convergence_velocity.png")

    # Plot 2: Pressure
    p_error_dict = {
        'Pressure': errors_p
    }
    # ### MODIFICATION: Added save_filename ###
    plot_convergence(h_values,
                     p_error_dict,
                     'Convergence of Pressure (Staggered)',
                     save_filename="../results/convergence_pressure.png")

    # Plot 3: Divergence
    # (This plot has been removed)

    print("\nAll plots saved successfully.")


if __name__ == "__main__":

    if len(sys.argv) > 1 and ".vtk" in sys.argv[1]:
        print("---")
        print("WARNING: This script is not run by dragging files onto it.")
        print("Please edit the 'SIMULATIONS' list inside the script")
        print("and then run it directly.")
        print("---")

    run_analysis()