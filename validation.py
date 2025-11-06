

import pyvista as pv
import numpy as np
import sys


VTK_INPUT_FILE = "solution.vtk" 
VELOCITY_FIELD_NAME = "U"
PRESSURE_FIELD_NAME = "p"

def compute_analytical_velocity(points: np.ndarray) -> np.ndarray:
    
    x = points[:, 0]
    y = points[:, 1]
    # z = points[:, 2] 
    
    pi = np.pi  
    u_analytical = np.sin(pi * x) * np.cos(pi * y)
    v_analytical = -np.cos(pi * x) * np.sin(pi * y)
    w_analytical = np.zeros_like(x)
    
    # Stack them into a (N, 3) vector array
    return np.stack((u_analytical, v_analytical, w_analytical), axis=1)


def compute_analytical_pressure(points: np.ndarray) -> np.ndarray:
   
    x = points[:, 0]
    y = points[:, 1]
    # z = points[:, 2] # Not used
    
    pi = np.pi
    p_analytical = np.sin(pi * x) * np.sin(pi * y)
    
    return p_analytical


def calculate_rms_error(numerical_data: np.ndarray, analytical_data: np.ndarray) -> float:
    
    # Calculate the difference
    error_field = numerical_data - analytical_data
    
    # Square the error. For vectors, this squares each component.
    error_field_sq = error_field**2
    
    if error_field.ndim == 2:
        error_field_sq = np.sum(error_field_sq, axis=1)
        

    mean_squared_error = np.mean(error_field_sq)
    return np.sqrt(mean_squared_error)


# --- 3. MAIN VALIDATION FUNCTION ---

def main():
   
    print("--- MMS Validation Script Starting ---")
    try:
        grid = pv.read(VTK_INPUT_FILE)
        print(f"Successfully loaded '{VTK_INPUT_FILE}'")
        print(f"Grid type: {grid.is_all_structured}")
        print(f"Total points: {grid.n_points}")
    except FileNotFoundError:
        print(f"ERROR: Input file not found: '{VTK_INPUT_FILE}'")
        print("Please check the VTK_INPUT_FILE variable.")
        sys.exit(1)
        
    points = grid.points
    
    # --- Velocity Validation ---
    print("\n--- Validating Velocity ---")
    try:
        # 1. Get numerical velocity from VTK
        vel_numerical = grid.point_data[VELOCITY_FIELD_NAME]
        print(f"Found numerical velocity field: '{VELOCITY_FIELD_NAME}'")
        
        # 2. Compute analytical velocity
        vel_analytical = compute_analytical_velocity(points)
        
        # 3. Calculate L2 RMS Error
        velocity_error = calculate_rms_error(vel_numerical, vel_analytical)
        
        # 4. Store error fields for visualization
        grid.point_data['velocity_analytical'] = vel_analytical
        grid.point_data['velocity_error_vector'] = vel_numerical - vel_analytical
        grid.point_data['velocity_error_magnitude'] = np.linalg.norm(
            grid.point_data['velocity_error_vector'], axis=1
        )
        
        print(f"Velocity L2 RMS Error: {velocity_error: .6e}")

    except KeyError:
        print(f"ERROR: Velocity field '{VELOCITY_FIELD_NAME}' not found in VTK file.")
        print(f"Available point data fields: {list(grid.point_data.keys())}")
        print("Please update the VELOCITY_FIELD_NAME variable.")
    except Exception as e:
        print(f"An unexpected error occurred during velocity validation: {e}")


    # --- Pressure Validation ---
    print("\n--- Validating Pressure ---")
    try:
        
        if PRESSURE_FIELD_NAME in grid.point_data:
            print(f"Found numerical pressure field: '{PRESSURE_FIELD_NAME}' (Point Data)")
            pressure_numerical = grid.point_data[PRESSURE_FIELD_NAME]
            pressure_points = points
        elif PRESSURE_FIELD_NAME in grid.cell_data:
            print(f"Found numerical pressure field: '{PRESSURE_FIELD_NAME}' (Cell Data)")
            pressure_numerical = grid.cell_data[PRESSURE_FIELD_NAME]
            pressure_points = grid.cell_centers().points
            print(f"Using {pressure_points.shape[0]} cell centers for pressure analysis.")
        else:
            raise KeyError(f"Field '{PRESSURE_FIELD_NAME}' not found in point or cell data.")

        
        # 2. Compute analytical pressure
        pressure_analytical = compute_analytical_pressure(pressure_points)
        
        # 3. Calculate L2 RMS Error
        pressure_error = calculate_rms_error(pressure_numerical, pressure_analytical)
        
        # 4. Store error fields for visualization
        # Note: We save the error back to the original grid data location
        if PRESSURE_FIELD_NAME in grid.point_data:
            grid.point_data['pressure_analytical'] = pressure_analytical
            grid.point_data['pressure_error'] = pressure_numerical - pressure_analytical
        elif PRESSURE_FIELD_NAME in grid.cell_data:
            grid.cell_data['pressure_analytical'] = pressure_analytical
            grid.cell_data['pressure_error'] = pressure_numerical - pressure_analytical

        print(f"Pressure L2 RMS Error: {pressure_error: .6e}")

    except KeyError:
        print(f"ERROR: Pressure field '{PRESSURE_FIELD_NAME}' not found in VTK file.")
        print(f"Available point data: {list(grid.point_data.keys())}")
        print(f"Available cell data:  {list(grid.cell_data.keys())}")
        print("Please update the PRESSURE_FIELD_NAME variable.")
    except Exception as e:
        print(f"An unexpected error occurred during pressure validation: {e}")


    print("--- Validation Script Finished ---")


if __name__ == "__main__":
    main()