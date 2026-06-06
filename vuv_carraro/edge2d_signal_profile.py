"""
author: R. Cicioni
last updated: 06-06-2026

This code plots the 2D signal map from a JT-60SA edge2d simulation.
This code can be run in both JDC and EFGW servers, but the path to the TRAN file must be updated accordingly (see RUN_DIR variable).
To change the quantity of interest, simply change the value of the 'PROFILE' variable to any other signal available in the TRAN file.
"""
import sys
import os
import matplotlib
matplotlib.use("TkAgg") 
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sys.path.append('/home/lgarzot/Work/Python/JET') #jintrac python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #jintrac python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import numpy as np
import jetto_binary_tools
import eproc as ep
import matplotlib.pyplot as plt
import datetime

#RUN_DIR = "/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/may2926/seq#2" #example run on JDC server
RUN_DIR = "/pfs/work/g2rcicio/edge2d/runs/run_sa_fk_c_WAr_p20_s19_a3_50ms_std" #example run on EFGW server
TRAN_FILE = os.path.join(RUN_DIR, "tran")
PROFILE = 'QPARTOT' #signal name for the profile (can be changed to any other signal available in the TRAN file)

rmesh=ep.data(TRAN_FILE,'RMESH') #centrid of the mesh in R direction
rmesh=rmesh.data[0:rmesh.nPts]
zmesh=ep.data(TRAN_FILE,'ZMESH') #centrid of the mesh in Z direction
zmesh=-zmesh.data[0:zmesh.nPts]

rvertp=ep.data(TRAN_FILE,'RVERTP') #vertices of the mesh in R direction (quintuple structure: 4 vertices + 1 centroid)
rvertp=rvertp.data[0:rvertp.nPts]
zvertp=ep.data(TRAN_FILE,'ZVERTP') #vertices of the mesh in Z direction (quintuple structure: 4 vertices + 1 centroid)
zvertp=-zvertp.data[0:zvertp.nPts]

profile_obj=ep.data(TRAN_FILE, PROFILE)
profile_name=profile_obj.desc if hasattr(profile_obj, 'desc') else PROFILE
profile_units=profile_obj.units if hasattr(profile_obj, 'units') else 'N/A'
profile=profile_obj.data[0:profile_obj.nPts]  # Extract the actual array data
print(f"Profile units: {profile_units}")

R0_obj = ep.data(TRAN_FILE, 'R0')
R0 = R0_obj.data[0]
print(f"R0 (major radius): {R0:.4f} m")

# Load separatrix position
R_sep_obj = ep.data(TRAN_FILE, 'RSEPX')
R_sep = R_sep_obj.data[0]
print(f"R_sep (separatrix): {R_sep:.4f} m")

# ============================================================================
# Extract nodes and quads from rvertp/zvertp quintuple structure
# ============================================================================
# Data structure: quintuple organization
# - Every 5 consecutive elements form one quad: 4 vertices + 1 centroid
# - Vertices are at indices [0,1,2,3], centroid at index [4]

n_total = len(rvertp)
n_quads = n_total // 5  # Number of quads

print(f"Extraxting mesh structure from rvertp/zvertp...")
print(f"Total number of points (rvertp,zvertp): {n_total}")
print(f"Number of quads: {n_quads}")

# Initialize nodes array (4 vertices per quad)
nodes = np.zeros((4 * n_quads, 2))
quads = np.zeros((n_quads, 4), dtype=int)

# Extract vertices (skip centroids which are at indices 5k+4)
node_idx = 0
for quad_idx in range(n_quads):
    # Extract the 4 vertices from this quintuple
    for vertex_in_quad in range(4):
        global_idx = 5 * quad_idx + vertex_in_quad
        nodes[node_idx, 0] = rvertp[global_idx]  # R coordinate
        nodes[node_idx, 1] = zvertp[global_idx]  # Z coordinate
        quads[quad_idx, vertex_in_quad] = node_idx
        node_idx += 1

# Printing info about the extracted mesh structure
print(f"\nExtracted mesh structure:")
print(f"  nodes shape: {nodes.shape}")
print(f"  quads shape: {quads.shape}")

print(f"\n Successfully extracted nodes and quads arrays.")
print(f"{'='*70}\n")

# ============================================================================
# Map centroids to mesh points (rmesh, zmesh) and extract PROFILE values
# ============================================================================
from scipy.spatial import cKDTree

# Extract centroids coordinates from rvertp/zvertp
centroid_indices = np.arange(4, len(rvertp), 5)
centroids = np.column_stack([rvertp[centroid_indices], zvertp[centroid_indices]])

# Build KDTree from mesh points for fast nearest-neighbor search
mesh_points = np.column_stack([rmesh, zmesh])
tree = cKDTree(mesh_points)

# Find nearest mesh point for each centroid and extract PROFILE values
distances, indices = tree.query(centroids)
values = profile[indices]

print(f"{profile_name} values mapping:")
print(f"  Number of centroids: {len(centroids)}")
print(f"  Number of mesh points (rmesh): {len(rmesh)}")
print(f"  {profile_name} values extracted: {len(values)}")
print(f"  {profile_name} range: {values.min():.2f} - {values.max():.2f}")
print(f"  Max distance to nearest mesh point: {distances.max():.2e}")

print(f"\nFirst 10 centroids and their {profile_name} values:")
print("-" * 72)
for i in range(min(10, len(centroids))):
    r_cent, z_cent = centroids[i]
    r_mesh, z_mesh = mesh_points[indices[i]]
    print(f"{i:>6}  {r_cent:>10.4f}  {z_cent:>10.4f}  {r_mesh:>10.4f}  {z_mesh:>10.4f}  {values[i]:>12.4f}  {distances[i]:>10.2e}")

print(f"\n Successfully created values array.")
print(f"{'='*70}\n")

# ============================================================================
# Line of Sight intersection algorithm (based on Sutherland-Hodgman algorithm)
# ============================================================================

def cross2d(a, b):
    """Scalar cross product of 2-vectors."""
    return a[0] * b[1] - a[1] * b[0]


def line_segment_intersect_t(origin, direction, p1, p2):
    """
    Intersect the infinite line X = origin + t*direction
    with the finite segment [p1, p2].
    Returns t and s parameters, or (None, None) if parallel.
    """
    d = direction
    v = p2 - p1
    denom = cross2d(d, v)
    if abs(denom) < 1e-12:
        return None, None
    w = p1 - origin
    t = cross2d(w, v) / denom
    s = cross2d(w, d) / denom
    return t, s


def ray_quad_interval(origin, direction, quad_nodes):
    """
    Find interval [t_enter, t_exit] where line is inside quadrilateral.
    Returns (t_enter, t_exit) or None if line misses quad.
    """
    n = len(quad_nodes)
    ts = []
    for k in range(n):
        p1 = quad_nodes[k]
        p2 = quad_nodes[(k + 1) % n]
        t, s = line_segment_intersect_t(origin, direction, p1, p2)
        if t is not None and -1e-9 <= s <= 1 + 1e-9:
            ts.append(t)

    if len(ts) < 2:
        return None

    ts = sorted(ts)
    t_enter, t_exit = ts[0], ts[-1]
    if t_exit - t_enter < 1e-10:
        return None
    return t_enter, t_exit


def extract_profile_los(nodes, quads, values, origin, direction):
    """
    For each quad, compute the segment of the line inside it.
    Returns sorted list of (t_enter, t_exit, value).
    """
    segments = []
    for qi, quad in enumerate(quads):
        quad_nodes = nodes[quad]
        result = ray_quad_interval(origin, direction, quad_nodes)
        if result is not None:
            t_enter, t_exit = result
            segments.append((t_enter, t_exit, values[qi]))

    segments.sort(key=lambda s: s[0])
    return segments


# ============================================================================
# Extract profile along outer mid-plane (Z=0, from R0 to R_max)
# ============================================================================
R_max = rmesh.max()
point1 = np.array([R0, 0.0])
point2 = np.array([R_max, 0.0])

origin = point1
direction = point2 - point1

print(f"\n{'='*70}")
print(f"LINE OF SIGHT {profile_name} PROFILE (OUTER MID-PLANE)")
print(f"{'='*70}")
print(f"Origin: {origin}")
print(f"Direction: {direction}")
print(f"Point 1 (R0): {point1}")
print(f"Point 2 (R_max): {point2}")

# Extract profile
segments = extract_profile_los(nodes, quads, values, origin, direction)

print(f"\nNumber of intersected cells: {len(segments)}")
print("-" * 42)
for t_enter, t_exit, val in segments:
    print(f"{t_enter:10.4f}  {t_exit:10.4f}  {t_exit-t_enter:10.4f}  {val:10.4f}")

print(f"\n Successfully extracted LoS profile.")
print(f"{'='*70}\n")

# ============================================================================
# Visualization: 2D mesh with signal + 1D mid-plane profile
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# --- Left: 2D mesh with line of sight ---
ax = axes[0]
ax.set_title(f"2D {profile_name} Map + Outer Mid-plane LoS", fontsize=12, fontweight='bold')
ax.set_aspect("equal")

norm = plt.Normalize(vmin=values.min(), vmax=values.max())
cmap = plt.cm.viridis

patches = []
colors = []
for qi, quad in enumerate(quads):
    poly = Polygon(nodes[quad], closed=True)
    patches.append(poly)
    colors.append(values[qi])

pc = PatchCollection(patches, cmap=cmap, norm=norm, edgecolors="black", linewidths=0.3, alpha=0.9)
pc.set_array(np.array(colors))
ax.add_collection(pc)
plt.colorbar(pc, ax=ax, label=f"{profile_name} ({profile_units})")

# Draw line of sight
t_all = [s[0] for s in segments] + [s[1] for s in segments]
if t_all:
    t_lo = min(t_all) - 0.1
    t_hi = max(t_all) + 0.1
else:
    t_lo, t_hi = -0.5, 1.5

los_start = origin + t_lo * direction
los_end = origin + t_hi * direction
ax.plot([los_start[0], los_end[0]], [los_start[1], los_end[1]],
        "r-", linewidth=2.5, label="Mid-plane LoS", zorder=10)

# Mark line endpoints
ax.plot(*point1, "go", markersize=12, label=f"R0 = {R0:.3f} m", zorder=11, markeredgecolor='darkgreen', markeredgewidth=2)
ax.plot(*point2, "mo", markersize=12, label=f"R_max = {R_max:.3f} m", zorder=11, markeredgecolor='purple', markeredgewidth=2)

# Mark separatrix position
z_min = nodes[:, 1].min()
z_max = nodes[:, 1].max()
ax.axvline(R_sep, color='cyan', linestyle='--', linewidth=2.5, label=f"Separatrix (R_sep = {R_sep:.3f} m)", zorder=9)

# Mark intersected quad centers
for t_enter, t_exit, val in segments:
    t_mid = 0.5 * (t_enter + t_exit)
    pt = origin + t_mid * direction
    ax.plot(*pt, "r.", markersize=8, zorder=12)

# Set axis limits to show all patches
ax.autoscale()
ax.autoscale_view()

ax.set_xlabel('R (m)')
ax.set_ylabel('Z (m)')
ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3)

# --- Right: 1D profile along mid-plane ---
ax2 = axes[1]
ax2.set_title(f"1D {profile_name} Along Outer Mid-plane", fontsize=12, fontweight='bold')

# Convert t parameters to R coordinates and plot as (R, value) 
R_positions = []
profile_values = []

for t_enter, t_exit, val in segments:
    R_mid = origin[0] + 0.5 * (t_enter + t_exit) * direction[0]
    
    # Only include if R is within valid range
    if R_mid >= R0 - 0.01 and R_mid <= R_max + 0.01:
        R_positions.append(R_mid)
        profile_values.append(-val)  # Flip sign to make positive

if len(R_positions) > 0:
    R_positions = np.array(R_positions)
    profile_values = np.array(profile_values)
    
    # Sort by R position
    sort_idx = np.argsort(R_positions)
    R_positions = R_positions[sort_idx]
    profile_values = profile_values[sort_idx]
    
    # Plot clean (no colors, just data)
    ax2.plot(R_positions, profile_values, 'o-', linewidth=2, markersize=6, color='black')
    
    # Mark separatrix position
    ax2.axvline(R_sep, color='cyan', linestyle='--', linewidth=2.5, label=f"Separatrix (R_sep = {R_sep:.3f} m)")
    
    ax2.set_xlabel('R (m)', fontsize=12)
    ax2.set_ylabel(f"{profile_name} ({profile_units})", fontsize=12)
    # Restrict xlim to actual data range
    R_min_data = R_positions.min()
    R_max_data = R_positions.max()
    margin = 0.02 * (R_max_data - R_min_data)
    ax2.set_xlim(R_min_data - margin, R_max_data + margin)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
else:
    ax2.text(0.5, 0.5, 'No valid profile data', ha='center', va='center', 
             transform=ax2.transAxes, fontsize=12)

plt.tight_layout()

# Try to show (may not work in all environments)
try:
    plt.show()
except Exception as e:
    print(f"Could not display plot interactively: {e}")

print(f"\n{profile_name} 2D map plot complete")
print(f"{'='*70}\n")

# ============================================================================
# Fit exponential decay profile in the SOL
# ============================================================================

def fit_exponential_decay(x, y, x_min=None, x_max=None, R_sep=None):
    """
    Fit exponential decay profile: q_parallel(x) = q0 * exp(-x/lambda_q)
    
    Parameters:
    -----------
    x : array-like
        Distance from separatrix (meters)
    y : array-like
        Parallel heat flux (W/m^2)
    x_min : float, optional
        Minimum x for fit window (default: min(x[x>0]))
    x_max : float, optional
        Maximum x for fit window (default: max(x))
    R_sep : float, optional
        Separatrix position (for reference)
    
    Returns:
    --------
    dict with keys:
        - 'lambda_q_m': lambda_q in meters
        - 'lambda_q_mm': lambda_q in millimeters
        - 'q0': q0 in W/m^2
        - 'slope': slope m in units of 1/m
        - 'intercept': intercept ln(q0)
        - 'R2': R-squared of the fit
        - 'n_points': number of points used in fit
        - 'x_fit': x values used in fit
        - 'y_fit': y values used in fit
    """
    
    # Filter valid data
    mask_valid = (x > 0) & (y != 0) & np.isfinite(y) & np.isfinite(x)
    x_valid = x[mask_valid]
    y_valid = y[mask_valid]
    
    print(f"\n{'='*70}")
    print(f"EXPONENTIAL DECAY FIT (q_parallel)")
    print(f"{'='*70}")
    print(f"Total points: {len(x)}")
    print(f"Valid points (x>0, y≠0, finite): {len(x_valid)}")
    
    if len(x_valid) < 3:
        print("ERROR: Not enough valid points for fit (need at least 3)")
        return None
    
    # Set fit window
    if x_min is None:
        x_min = x_valid.min()
    if x_max is None:
        x_max = x_valid.max()
    
    mask_fit = (x_valid >= x_min) & (x_valid <= x_max)
    x_fit = x_valid[mask_fit]
    y_fit = y_valid[mask_fit]
    
    print(f"\nFit window: x ∈ [{x_min:.4f}, {x_max:.4f}] m")
    print(f"Points in fit window: {len(x_fit)}")
    
    if len(x_fit) < 3:
        print("ERROR: Not enough points in fit window (need at least 3)")
        return None
    
    # Linear fit of ln|y| vs x
    ln_y_fit = np.log(np.abs(y_fit))
    
    # Coefficients: [slope, intercept]
    coeffs = np.polyfit(x_fit, ln_y_fit, 1)
    slope = coeffs[0]  # m = -1/lambda_q
    intercept = coeffs[1]  # ln(q0)
    
    q0 = np.exp(intercept)
    
    # Check if slope is negative (physical requirement)
    if slope >= 0:
        print("WARNING: Slope is non-negative. This is not physical for heat flux decay.")
        print(f"Slope = {slope:.6f}")
        return None
    
    lambda_q = -1.0 / slope  # lambda_q in meters
    lambda_q_mm = lambda_q * 1000  # in millimeters
    
    # Calculate R^2
    y_pred = np.exp(intercept + slope * x_fit)
    ss_res = np.sum((np.abs(y_fit) - y_pred) ** 2)
    ss_tot = np.sum((np.abs(y_fit) - np.mean(np.abs(y_fit))) ** 2)
    R2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else -999
    
    print(f"\nFit Results:")
    print(f"  q0 = {q0:.4e} W/m²")
    print(f"  Slope m = {slope:.6e} m⁻¹")
    print(f"  λ_q = {lambda_q:.4f} m = {lambda_q_mm:.2f} mm")
    print(f"  R² = {R2:.6f}")
    
    result = {
        'lambda_q_m': lambda_q,
        'lambda_q_mm': lambda_q_mm,
        'q0': q0,
        'slope': slope,
        'intercept': intercept,
        'R2': R2,
        'n_points': len(x_fit),
        'x_fit': x_fit,
        'y_fit': y_fit,
        'x_valid': x_valid,
        'y_valid': y_valid
    }
    
    return result


def plot_exponential_fit(fit_result, profile_name, profile_units, R_sep=None, x_range=None):
    """
    Plot the original data and exponential fit.
    
    Parameters:
    -----------
    fit_result : dict
        Output from fit_exponential_decay
    profile_name : str
        Name of the profile
    profile_units : str
        Units of the profile
    R_sep : float, optional
        Separatrix position (for reference)
    x_range : tuple, optional
        (x_min, x_max) for plot range
    """
    
    if fit_result is None:
        return
    
    x_valid = fit_result['x_valid']
    y_valid = fit_result['y_valid']
    x_fit = fit_result['x_fit']
    y_fit = fit_result['y_fit']
    
    q0 = fit_result['q0']
    slope = fit_result['slope']
    lambda_q = fit_result['lambda_q_m']
    R2 = fit_result['R2']
    
    # Generate smooth curve for fit
    x_smooth = np.linspace(x_fit.min(), x_fit.max(), 200)
    y_fit_smooth = q0 * np.exp(-x_smooth / lambda_q)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot data
    ax.semilogy(x_valid, np.abs(y_valid), 'o', markersize=8, color='black', 
                label='Data', alpha=0.7)
    
    # Plot fit
    ax.semilogy(x_smooth, y_fit_smooth, '-', linewidth=2.5, color='red',
                label=f'Fit: $q_0$ exp$(-x/\\lambda_q)$')
    
    # Formatting
    ax.set_xlabel(f'Distance from separatrix (m)', fontsize=12)
    ax.set_ylabel(f'{profile_name} ({profile_units})', fontsize=12)
    ax.set_title(f'Exponential decay: $\\lambda_q$ = {lambda_q*1000:.1f} mm (R² = {R2:.4f})', 
                 fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=11, loc='upper right')
    
    if x_range is not None:
        ax.set_xlim(x_range[0], x_range[1])
    
    plt.tight_layout()
    return fig


# ============================================================================
# Apply exponential fit to mid-plane profile data
# ============================================================================

if len(R_positions) > 0:
    # Filter data for R > R_sep (SOL only)
    sol_mask = R_positions > R_sep
    R_positions_sol = R_positions[sol_mask]
    profile_values_sol = profile_values[sol_mask]
    
    if len(R_positions_sol) > 0:
        # Prepare data: x = R - R_sep (distance from separatrix)
        x_data = R_positions_sol - R_sep
        y_data = profile_values_sol
        
        print(f"\n{'='*70}")
        print(f"MID-PLANE PROFILE DATA FOR SOL FIT")
        print(f"{'='*70}")
        print(f"R_sep = {R_sep:.4f} m")
        print(f"R range (SOL only, R > R_sep): {R_positions_sol.min():.4f} - {R_positions_sol.max():.4f} m")
        print(f"x = R - R_sep range: {x_data.min():.4f} - {x_data.max():.4f} m")
        print(f"y ({profile_name}) range: {y_data.min():.4e} - {y_data.max():.4e} {profile_units}")
        
        # Perform fit (user can adjust x_min, x_max for SOL window)
        # For now, use a reasonable default: take points in SOL
        x_min_fit = 0.0
        x_max_fit = x_data.max()
        
        fit_result = fit_exponential_decay(x_data, y_data, x_min=x_min_fit, x_max=x_max_fit)
        
        if fit_result is not None:
            # Plot the fit
            fig_fit = plot_exponential_fit(fit_result, profile_name, profile_units)
            
            try:
                plt.show()
            except Exception as e:
                print(f"Could not display fit plot: {e}")
        else:
            print("Exponential fit could not be performed.")
    else:
        print(f"No data in SOL (R > R_sep = {R_sep:.4f} m)")
