"""
author: R. Cicioni
last updated: 29-05-2025

Module for extracting and visualizing Line of Sight (LoS) profiles from JETTO/JINTRAC data.
Can be used as a standalone script or imported as a module.
"""
import sys
import os
import matplotlib
matplotlib.use("TkAgg") 
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sys.path.append('/home/lgarzot/Work/Python/JET') #python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import numpy as np
import jetto_binary_tools
import eproc as ep
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree


# ============================================================================
# HELPER FUNCTIONS FOR LINE OF SIGHT ANALYSIS
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


def extract_profile(nodes, quads, values, origin, direction):
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


def extract_mesh_structure(rvertp, zvertp):
    """
    Extract nodes and quads from rvertp/zvertp quintuple structure.
    
    Data structure: quintuple organization
    - Every 5 consecutive elements form one quad: 4 vertices + 1 centroid
    - Vertices are at indices [0,1,2,3], centroid at index [4]
    
    Returns:
        nodes (ndarray): Array of shape (4*n_quads, 2) with [R, Z] coordinates
        quads (ndarray): Array of shape (n_quads, 4) with vertex indices
    """
    n_total = len(rvertp)
    n_quads = n_total // 5

    print(f"Extracting mesh structure from rvertp/zvertp...")
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

    print(f"✓ Successfully extracted {len(nodes)} nodes and {len(quads)} quads.\n")
    return nodes, quads


def map_values_to_mesh(rmesh, zmesh, rvertp, zvertp, values):
    """
    Map values from mesh points to quad centroids using nearest-neighbor search.
    
    Args:
        rmesh, zmesh: Mesh point coordinates
        rvertp, zvertp: Mesh vertex coordinates (quintuple structure)
        values: Values at mesh points
        
    Returns:
        mapped_values: Values at each quad's centroid
    """
    # Extract centroids coordinates from rvertp/zvertp
    centroid_indices = np.arange(4, len(rvertp), 5)
    centroids = np.column_stack([rvertp[centroid_indices], zvertp[centroid_indices]])

    # Build KDTree from mesh points for fast nearest-neighbor search
    mesh_points = np.column_stack([rmesh, zmesh])
    tree = cKDTree(mesh_points)

    # Find nearest mesh point for each centroid and extract values
    distances, indices = tree.query(centroids)
    mapped_values = values[indices]

    print(f"Value mapping statistics:")
    print(f"  Number of centroids: {len(centroids)}")
    print(f"  Number of mesh points: {len(rmesh)}")
    print(f"  Values range: {mapped_values.min():.2f} - {mapped_values.max():.2f}")
    print(f"  Max distance to nearest mesh point: {distances.max():.2e}\n")

    return mapped_values


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def extract_los_profile(tran_file, quantity_name, point1, point2, 
                       plot=True, verbose=True):
    """
    Extract and analyze Line of Sight profile from JETTO/JINTRAC data.
    
    Args:
        tran_file (str): Path to TRAN file
        quantity_name (str): Name of the quantity to extract (e.g., 'TEVE', 'NEDE')
        point1 (array-like): Starting point [R, Z] of the line of sight
        point2 (array-like): Ending point [R, Z] of the line of sight
        plot (bool): Whether to create visualization plots
        verbose (bool): Whether to print detailed information
        
    Returns:
        dict: Dictionary containing:
            - 'segments': List of (t_enter, t_exit, value) tuples
            - 'nodes': Mesh nodes array
            - 'quads': Mesh quads array
            - 'values': Values at each quad
            - 'quantity_name': Name of the extracted quantity
            - 'origin': Line of sight origin
            - 'direction': Line of sight direction
    """
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"JINTRAC Line of Sight Profile Extraction")
        print(f"{'='*70}")
        print(f"TRAN file: {tran_file}")
        print(f"Quantity: {quantity_name}")
        print(f"Point 1 (origin): {point1}")
        print(f"Point 2 (direction): {point2}\n")
    
    # Load mesh and quantity data
    rmesh = ep.data(tran_file, 'RMESH')
    rmesh = rmesh.data[0:rmesh.nPts]
    
    zmesh = ep.data(tran_file, 'ZMESH')
    zmesh = -zmesh.data[0:zmesh.nPts]
    
    rvertp = ep.data(tran_file, 'RVERTP')
    rvertp = rvertp.data[0:rvertp.nPts]
    
    zvertp = ep.data(tran_file, 'ZVERTP')
    zvertp = -zvertp.data[0:zvertp.nPts]
    
    # Load the specified quantity
    quantity_data = ep.data(tran_file, quantity_name)
    quantity_values = quantity_data.data[0:quantity_data.nPts]
    
    if verbose:
        print(f"Loaded {quantity_name}: {len(quantity_values)} values")
        print(f"Range: {quantity_values.min():.2e} - {quantity_values.max():.2e}\n")
    
    # Extract mesh structure
    nodes, quads = extract_mesh_structure(rvertp, zvertp)
    
    # Map quantity values to quad centroids
    values = map_values_to_mesh(rmesh, zmesh, rvertp, zvertp, quantity_values)
    
    # Convert to parametric form: origin + t*direction
    point1 = np.array(point1)
    point2 = np.array(point2)
    origin = point1
    direction = point2 - point1
    
    if verbose:
        print(f"{'='*70}")
        print(f"LINE OF SIGHT TEMPERATURE PROFILE")
        print(f"{'='*70}")
        print(f"Origin: {origin}")
        print(f"Direction: {direction}")
        print(f"Direction magnitude: {np.linalg.norm(direction):.4f}\n")
    
    # Extract profile
    segments = extract_profile(nodes, quads, values, origin, direction)
    
    if verbose:
        print(f"Number of intersected cells: {len(segments)}")
        if len(segments) > 0:
            print(f"\n{'t_enter':>10}  {'t_exit':>10}  {'Δt':>10}  {quantity_name:>12}")
            print("-" * 44)
            for t_enter, t_exit, val in segments:
                print(f"{t_enter:10.4f}  {t_exit:10.4f}  {t_exit-t_enter:10.4f}  {val:12.4e}")
        print()
    
    # Visualization
    if plot:
        visualize_los_profile(nodes, quads, values, segments, origin, direction,
                             point1, point2, quantity_name)
    
    return {
        'segments': segments,
        'nodes': nodes,
        'quads': quads,
        'values': values,
        'quantity_name': quantity_name,
        'origin': origin,
        'direction': direction,
        'rmesh': rmesh,
        'zmesh': zmesh,
    }


def visualize_los_profile(nodes, quads, values, segments, origin, direction,
                         point1, point2, quantity_name):
    """
    Create visualization of 2D mesh with LoS and 1D profile plot.
    
    Args:
        nodes, quads: Mesh structure
        values: Values at each quad
        segments: List of (t_enter, t_exit, value) tuples
        origin, direction: Line of sight parametric equation
        point1, point2: Line of sight endpoints
        quantity_name: Name of the plotted quantity
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # --- Left: 2D mesh with line of sight ---
    ax = axes[0]
    ax.set_title("Mesh + Line of Sight", fontsize=12, fontweight='bold')
    ax.set_aspect("equal")

    norm = plt.Normalize(vmin=values.min(), vmax=values.max())
    cmap = plt.cm.viridis

    patches = []
    colors = []
    for qi, quad in enumerate(quads):
        poly = Polygon(nodes[quad], closed=True)
        patches.append(poly)
        colors.append(values[qi])

    pc = PatchCollection(patches, cmap=cmap, norm=norm, edgecolors="black", 
                        linewidths=0.3, alpha=0.9)
    pc.set_array(np.array(colors))
    ax.add_collection(pc)
    plt.colorbar(pc, ax=ax, label=f"{quantity_name}")

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
            "r-", linewidth=2.5, label="Line of Sight", zorder=10)

    # Mark line endpoints
    ax.plot(*point1, "go", markersize=12, label=f"Point 1", zorder=11, 
           markeredgecolor='darkgreen', markeredgewidth=2)
    ax.plot(*point2, "mo", markersize=12, label=f"Point 2", zorder=11, 
           markeredgecolor='purple', markeredgewidth=2)

    # Mark intersected quad centers
    for t_enter, t_exit, val in segments:
        t_mid = 0.5 * (t_enter + t_exit)
        pt = origin + t_mid * direction
        ax.plot(*pt, "r.", markersize=8, zorder=12)

    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # --- Right: 1D profile ---
    ax2 = axes[1]
    ax2.set_title(f"1D Profile Along Line of Sight", fontsize=12, fontweight='bold')

    for t_enter, t_exit, val in segments:
        color = cmap(norm(val))
        ax2.fill_between([t_enter, t_exit], [val, val], step="post", alpha=0.7, color=color)
        ax2.hlines(val, t_enter, t_exit, colors=color, linewidths=2.5)
        ax2.vlines([t_enter, t_exit], 0, val, colors="gray", linewidths=0.6, linestyles="--")

    ax2.set_xlabel("Parameter t along ray")
    ax2.set_ylabel(f"{quantity_name}")
    ax2.set_xlim(t_lo, t_hi)
    if len(values) > 0:
        ax2.set_ylim(values.min() - 0.05 * (values.max() - values.min()), 
                     values.max() + 0.1 * (values.max() - values.min()))
    ax2.grid(True, alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax2, label=f"{quantity_name}")

    plt.tight_layout()
    plt.show()


# ============================================================================
# MAIN SCRIPT EXECUTION
# ============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Extract Line of Sight profile from JETTO/JINTRAC TRAN file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python jintrac_los_profile.py --tran /path/to/tran --quantity TEVE --point1 4.5 0.0 --point2 3.0 2.0
        """)
    
    parser.add_argument('--tran', type=str, required=True,
                       help='Path to TRAN file')
    parser.add_argument('--quantity', type=str, default='TEVE',
                       help='Quantity name to extract (default: TEVE)')
    parser.add_argument('--point1', type=float, nargs=2, default=[4.5, 0.0],
                       help='Starting point coordinates [R Z] (default: 4.5 0.0)')
    parser.add_argument('--point2', type=float, nargs=2, default=[3.0, 2.0],
                       help='Ending point coordinates [R Z] (default: 3.0 2.0)')
    parser.add_argument('--no-plot', action='store_true',
                       help='Disable plotting')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress verbose output')
    
    args = parser.parse_args()
    
    result = extract_los_profile(
        tran_file=args.tran,
        quantity_name=args.quantity,
        point1=args.point1,
        point2=args.point2,
        plot=not args.no_plot,
        verbose=not args.quiet
    )




