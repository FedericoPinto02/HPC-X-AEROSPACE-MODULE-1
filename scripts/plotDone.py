import numpy as np
import matplotlib.pyplot as plt

# --- 1. Simulation Data ---

# List of spatial nodes (NX) from the error analysis (up to Nx=50 in this dataset)
NX_LIST = np.array([20, 30, 40, 50])

# RMS Errors for all fields at T_END = 0.05
# Note: The pressure error (RMS_p) uses a separate scale and convergence behavior
# compared to velocity fields (RMS_u, RMS_v, RMS_w, RMS_Mag).
RMS_ERRORS = {
    # Velocity components (RMS u, v, w) and Magnitude (RMS Mag)
    'RMS_u':   np.array([6.015058e-02, 4.718574e-02, 4.122232e-02, 3.787648e-02]),
    'RMS_v':   np.array([1.255041e-01, 1.271701e-01, 1.268702e-01, 1.263124e-01]),
    'RMS_w':   np.array([1.126325e-01, 1.132682e-01, 1.116860e-01, 1.101379e-01]),
    'RMS_Mag': np.array([1.790404e-01, 1.767156e-01, 1.739801e-01, 1.718133e-01]),

    # Pressure error (RMS p)
    'RMS_p':   np.array([1.673516e+00, 1.253976e+00, 1.025882e+00, 8.947853e-01]),
}



# Grid spacing h is proportional to 1/NX (assuming a fixed domain length L=1)
H_LIST = 1.0 / NX_LIST

# Define plot styles for each series
plot_styles = {
    'RMS_u': {'color': 'blue', 'marker': 'o', 'label': 'RMS u'},
    'RMS_v': {'color': 'orange', 'marker': 's', 'label': 'RMS v'},
    'RMS_w': {'color': 'purple', 'marker': '^', 'label': 'RMS w'},
    'RMS_p': {'color': 'black', 'marker': 'D', 'label': 'RMS p'},
    'RMS_Mag': {'color': 'red', 'marker': 'x', 'label': 'RMS Mag'}
}

# Define velocity and pressure fields for separate plots
VELOCITY_FIELDS = ['RMS_u', 'RMS_v', 'RMS_w', 'RMS_Mag']
PRESSURE_FIELD = 'RMS_p'

# --- 2. Plotting Function ---

def plot_convergence(fields_to_plot, title, anchor_field):
    """Generates a log-log convergence plot for the specified fields."""
    
    # Anchor the reference lines to the coarsest data point of the specified anchor field
    h_ref = H_LIST[0]
    error_ref = RMS_ERRORS[anchor_field][0]

    # First-Order Convergence (O(h^1)): Error = C * h^1
    ref_order_1 = error_ref * (H_LIST / h_ref)

    # Second-Order Convergence (O(h^2)): Error = C * h^2
    ref_order_2 = error_ref * (H_LIST / h_ref)**2
    
    # Set up the figure and axes
    plt.figure(figsize=(10, 7))
    plt.title(title, fontsize=14)
    plt.xlabel('Grid Spacing $h = 1/N_x$', fontsize=12)
    plt.ylabel('RMS Error (Log Scale)', fontsize=12)

    # Plot the simulation data (Convergence Curves)
    for label in fields_to_plot:
        errors = RMS_ERRORS[label]
        style = plot_styles[label]
        # Use '-' for velocity/magnitude fields, use '--' for pressure (if it was grouped)
        line_style = '-'
        if label == 'RMS_p':
            line_style = '-' 
        
        plt.loglog(H_LIST, errors,
                   marker=style['marker'],
                   linestyle=line_style,
                   color=style['color'],
                   linewidth=2,
                   markersize=8,
                   label=style['label'])

    # Plot the First-Order reference line (Slope = 1) - Emphasized
    plt.loglog(H_LIST, ref_order_1,
               '--',          # Dashed line
               color='gray', # Changed color to be neutral
               linewidth=3,   # Increased linewidth for emphasis
               label='First Order ($O(h^1)$)')

    # Plot the Second-Order reference line (Slope = 2) - Emphasized
    plt.loglog(H_LIST, ref_order_2,
               '-.',          # Dash-dot line
               color='red',   # Changed color for emphasis
               linewidth=3,   # Increased linewidth for emphasis
               label='Second Order ($O(h^2)$)')

    # Add grid and legend
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=10, loc='lower left')
    plt.gca().invert_xaxis() # Invert x-axis so that convergence (increasing Nx) moves left to right
    
    plt.tight_layout()
    plt.show()

# --- 3. Generate Plots ---

# Plot 1: Velocity Errors
# Anchor reference lines to the velocity magnitude error
plot_convergence(
    fields_to_plot=VELOCITY_FIELDS,
    title='Convergence Study: Velocity Errors vs. Grid Spacing $h$',
    anchor_field='RMS_Mag'
)

# Plot 2: Pressure Error
# Anchor reference lines to the pressure error
plot_convergence(
    fields_to_plot=[PRESSURE_FIELD],
    title='Convergence Study: Pressure Error (RMS P) vs. Grid Spacing $h$',
    anchor_field='RMS_p'
)