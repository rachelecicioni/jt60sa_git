"""
author: R. Cicioni
last updated: 29-05-2025

This code reconstruct the LoS signal from the 2D temperature map of a JT-60SA edge2d simulation.
This code can be run in both JDC and EFGW servers, but the path to the TRAN file must be updated accordingly (see RUN_DIR variable).
1) To change the quantity of interest, simply change the value of the 'PROFILE' variable to any other signal available in the TRAN file.
2) To change the line of sight, simply change the coordinates of 'point1' and 'point2' variables.
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

RUN_DIR = "/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/may2926/seq#2" #example run on JDC server
#RUN_DIR = "/pfs/work/g2rcicio/edge2d/runs/run_sa_fk_c_WAr_p20_s19_a3_50ms_std" #example run on EFGW server
TRAN_FILE = os.path.join(RUN_DIR, "tran")
PROFILE = 'CHII' #signal name for the temperature profile (can be changed to any other signal available in the TRAN file)
# Define line of sight endpoints
point1 = np.array([4.5, 0.0]) #line of sight endpoint 1 (can be changed to any other coordinates)
point2 = np.array([3.0, 2.0]) #line of sight endpoint 2 (can be changed to any other coordinates)

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

'''
print(f"\n NODES (first 20 rows) ")
print(f"{'Index':>6}  {'R (m)':>14}  {'Z (m)':>14}")
print("-" * 38)
for i in range(min(20, len(nodes))):
    print(f"{i:>6}  {nodes[i, 0]:>14.6f}  {nodes[i, 1]:>14.6f}")

print(f"\n QUADS (first 10 rows) ")
print(f"{'Quad':>6}  {'Vertex 0':>10}  {'Vertex 1':>10}  {'Vertex 2':>10}  {'Vertex 3':>10}")
print("-" * 54)
for i in range(min(10, len(quads))):
    print(f"{i:>6}  {quads[i, 0]:>10}  {quads[i, 1]:>10}  {quads[i, 2]:>10}  {quads[i, 3]:>10}")
'''
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
#print(f"{'Quad':>6}  {'R cent':>10}  {'Z cent':>10}  {'R mesh':>10}  {'Z mesh':>10}  {{{profile_name}}:>12}  {'Dist':>10}")
print("-" * 72)
for i in range(min(10, len(centroids))):
    r_cent, z_cent = centroids[i]
    r_mesh, z_mesh = mesh_points[indices[i]]
    print(f"{i:>6}  {r_cent:>10.4f}  {z_cent:>10.4f}  {r_mesh:>10.4f}  {z_mesh:>10.4f}  {values[i]:>12.4f}  {distances[i]:>10.2e}")

print(f"\n Successfully created values array.")
print(f"{'='*70}\n")

# ============================================================================
# LoS intersection algorithm (bsed on Sutherland-Hodgman algorithm)
# ============================================================================

def cross2d(a, b):
    #Scalar cross product of 2-vectors.
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


# ============================================================================
# Apply LoS to extract profile along the ray defined by point1 and point2
# ============================================================================
# Convert to parametric form: origin + t*direction
origin = point1
direction = point2 - point1

print(f"\n{'='*70}")
print(f"LINE OF SIGHT {profile_name} PROFILE")
print(f"{'='*70}")
print(f"Origin: {origin}")
print(f"Direction: {direction}")
print(f"Point 1: {point1}")
print(f"Point 2: {point2}")

# Extract profile
segments = extract_profile(nodes, quads, values, origin, direction)

print(f"\nNumber of intersected cells: {len(segments)}")
#print(f"\n{'t_enter':>10}  {'t_exit':>10}  {'Δt':>10}  {{{profile_name}}:>10}")
print("-" * 42)
for t_enter, t_exit, val in segments:
    print(f"{t_enter:10.4f}  {t_exit:10.4f}  {t_exit-t_enter:10.4f}  {val:10.4f}")

# ============================================================================
# Visualization: mesh with line of sight + 1D temperature profile
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# --- Left: 2D mesh with line of sight ---
ax = axes[0]
ax.set_title(f"{profile_name} + Line of Sight", fontsize=12, fontweight='bold')
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
        "r-", linewidth=2.5, label="Line of Sight", zorder=10)

# Mark line endpoints
ax.plot(*point1, "go", markersize=12, label="Point 1", zorder=11, markeredgecolor='darkgreen', markeredgewidth=2)
ax.plot(*point2, "mo", markersize=12, label="Point 2", zorder=11, markeredgecolor='purple', markeredgewidth=2)

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
ax2.set_title(f"1D {profile_name} Along Line of Sight", fontsize=12, fontweight='bold')

for t_enter, t_exit, val in segments:
    color = cmap(norm(val))
    ax2.fill_between([t_enter, t_exit], [val, val], step="post", alpha=0.7, color=color)
    ax2.hlines(val, t_enter, t_exit, colors=color, linewidths=2.5)
    ax2.vlines([t_enter, t_exit], 0, val, colors="gray", linewidths=0.6, linestyles="--")

ax2.set_xlabel("Parameter t along ray")
ax2.set_ylabel(f"{profile_name} ({profile_units})")
ax2.set_xlim(t_lo, t_hi)
ax2.set_ylim(values.min() - 0.05 * (values.max() - values.min()), 
             values.max() + 0.1 * (values.max() - values.min()))
ax2.grid(True, alpha=0.3)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ax=ax2, label=f"{profile_name} ({profile_units})")

plt.tight_layout()
plt.show()

print(f"\n{profile_name} Line of Sight analysis complete")
print(f"{'='*70}\n")


