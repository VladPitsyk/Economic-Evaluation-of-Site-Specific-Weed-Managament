"""
Rotated grid generation and per-cell weed count aggregation.



Description:
    Generates a field-aligned rectangular grid over a field boundary polygon,
    aligning grid cell orientation to the crop row direction defined by a tramline
    reference layer. Individual weed detection points (dicots and monocots) from
    UAV-derived weed maps are aggregated by spatial join to produce per-cell weed
    counts. The resulting grid layer is exported as a GeoPackage file ready for
    import into QGIS, where it serves as input to the bio-economic simulation
    (02_bioeconomic_simulation.py).

    The grid is constructed by rotating the field boundary to axis-aligned
    orientation, generating a regular rectangular grid over the bounding box,
    then rotating the grid back to match the field orientation. This ensures
    grid cells align with the direction of crop rows, which is the relevant
    spatial unit for site-specific herbicide application.

Inputs:
    - Field boundary polygon (shapefile or GeoPackage)
    - Tramline reference line defining crop row direction (shapefile or GeoPackage)
    - Dicot weed detection point layer (shapefile or GeoPackage)
    - Monocot weed detection point layer (shapefile or GeoPackage)

Outputs:
    - GeoPackage containing the rotated grid with per-cell weed counts:
        DicotCount  : number of detected dicot plants per cell
        MonoCount   : number of detected monocot plants per cell
        TotalCount  : sum of DicotCount and MonoCount
        Inflevel    : categorical infestation level (None / Low / Medium /
                      High / Very High), based on TotalCount thresholds

Requirements:
    Python >= 3.8
    geopandas >= 0.12, numpy, pandas, shapely

Usage:
    Set all paths in the USER CONFIGURATION section, then run:
        python 01_grid_creation_weed_counts.py
"""

from pathlib import Path
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Folder containing all input files
DATA_FOLDER = Path("C:/path/to/field/data/")

# Input filenames (shapefile or GeoPackage)
FIELD_BOUNDARY_FILE = "field_boundary.shp"   # Field boundary polygon
TRAMLINE_FILE = "tramline.gpkg"         # Reference line for crop row direction
# Dicot weed detection points (Amaranthus spp.)
DICOT_FILE = "dicot_detections.shp"
# Monocot weed detection points (E. crus-galli)
MONOCOT_FILE = "monocot_detections.shp"

# Grid cell size in metres. Tested at 1.0, 0.50, and 0.25 m in the paper.
CELL_SIZE = 0.25

# Target projected CRS for metric calculations (UTM Zone 33N, EPSG:32633)
# Adjust to the appropriate UTM zone for your study area.
# All input layers are reprojected to this CRS before processing.
TARGET_CRS = "EPSG:32633"


# =============================================================================
# FUNCTIONS
# =============================================================================

def get_rotation_angle(gdf_line):
    """
    Calculate the orientation angle of the crop row direction from a tramline.

    For MultiLineString inputs, the longest component is used to determine
    the primary row direction.

    Parameters
    ----------
    gdf_line : GeoDataFrame
        Single-feature layer containing the tramline reference geometry.

    Returns
    -------
    float : Angle in degrees (clockwise from east).
    """
    line = gdf_line.geometry.iloc[0]
    if line.geom_type == "MultiLineString":
        line = max(line.geoms, key=lambda g: g.length)
    x0, y0 = line.coords[0]
    x1, y1 = line.coords[-1]
    return np.degrees(np.arctan2(y1 - y0, x1 - x0))


def create_rotated_grid(gdf_boundary, angle_deg, cell_size):
    """
    Generate a rectangular grid aligned to the crop row direction.

    The field boundary is temporarily rotated to axis-aligned orientation,
    a regular grid is generated over the bounding box, and the grid is then
    rotated back to match the original field orientation.

    Parameters
    ----------
    gdf_boundary : GeoDataFrame
        Field boundary polygon layer.
    angle_deg    : float
        Rotation angle in degrees (output of get_rotation_angle).
    cell_size    : float
        Grid cell side length in metres.

    Returns
    -------
    GeoDataFrame : Grid cell polygons with the same CRS as gdf_boundary.
    """
    centroid = gdf_boundary.geometry.union_all().centroid

    # Rotate boundary to axis-aligned orientation for grid generation
    boundary_rot = gdf_boundary.rotate(-angle_deg, origin=centroid)
    minx, miny, maxx, maxy = boundary_rot.total_bounds

    # Generate regular grid over bounding box
    xs = np.arange(np.floor(minx), np.ceil(maxx), cell_size)
    ys = np.arange(np.floor(miny), np.ceil(maxy), cell_size)

    polygons = [
        Polygon([(x, y), (x + cell_size, y),
                 (x + cell_size, y + cell_size), (x, y + cell_size)])
        for x in xs for y in ys
    ]

    # Rotate grid back to original field orientation
    grid = gpd.GeoDataFrame({'geometry': polygons}, crs=gdf_boundary.crs)
    grid = grid.rotate(angle_deg, origin=centroid)
    return gpd.GeoDataFrame({'geometry': grid}, crs=gdf_boundary.crs)


def classify_infestation(row):
    """
    Assign a categorical infestation level based on total weed count per cell.

    Thresholds represent the total number of detected weed plants (both classes
    combined) per grid cell, regardless of cell size.

    Parameters
    ----------
    row : pd.Series
        DataFrame row containing a 'TotalCount' column.

    Returns
    -------
    str : Infestation level label.
    """
    total = row['TotalCount']
    if total >= 15:
        return 'Very High'
    elif total >= 10:
        return 'High'
    elif total >= 5:
        return 'Medium'
    elif total >= 1:
        return 'Low'
    return 'None'


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":

    # --- Load input layers ---
    print("Loading input layers...")
    gdf_boundary = gpd.read_file(DATA_FOLDER / FIELD_BOUNDARY_FILE)
    gdf_line = gpd.read_file(DATA_FOLDER / TRAMLINE_FILE)

    # Weed detection layers are optional; empty layers are substituted if absent
    try:
        gdf_dicots = gpd.read_file(DATA_FOLDER / DICOT_FILE)
    except Exception:
        print("  No dicot detection file found — assuming zero dicots.")
        gdf_dicots = gpd.GeoDataFrame(geometry=[], crs=gdf_boundary.crs)

    try:
        gdf_monocots = gpd.read_file(DATA_FOLDER / MONOCOT_FILE)
    except Exception:
        print("  No monocot detection file found — assuming zero monocots.")
        gdf_monocots = gpd.GeoDataFrame(geometry=[], crs=gdf_boundary.crs)

    # --- Reproject all layers to metric CRS ---
    # Required so that CELL_SIZE is interpreted in metres rather than degrees.
    print(f"Reprojecting to {TARGET_CRS}...")
    for name, gdf in [("boundary", gdf_boundary), ("tramline", gdf_line)]:
        if gdf.crs.to_epsg() != int(TARGET_CRS.split(":")[1]):
            print(f"  Converting {name} from {gdf.crs} to {TARGET_CRS}")

    gdf_boundary = gdf_boundary.to_crs(TARGET_CRS)
    gdf_line = gdf_line.to_crs(TARGET_CRS)
    if not gdf_dicots.empty:
        gdf_dicots = gdf_dicots.to_crs(TARGET_CRS)
    if not gdf_monocots.empty:
        gdf_monocots = gdf_monocots.to_crs(TARGET_CRS)

    # --- Generate rotated grid ---
    print(f"Generating {CELL_SIZE} m grid aligned to crop row direction...")
    angle = get_rotation_angle(gdf_line)
    print(f"  Detected crop row angle: {angle:.2f} degrees")

    grid = create_rotated_grid(gdf_boundary, angle, CELL_SIZE)

    # Clip grid to field boundary
    print("  Clipping grid to field boundary...")
    grid = gpd.clip(grid, gdf_boundary)

    # Assign unique cell identifier
    grid['grid_id'] = range(1, len(grid) + 1)
    print(f"  Grid contains {len(grid)} cells after clipping.")

    # --- Count weed detections per cell ---
    # Spatial join assigns each point detection to the grid cell it falls in.
    # Points on cell boundaries are assigned to both adjacent cells (intersects
    # predicate), which may marginally inflate counts at boundaries.

    print("Counting dicot detections per cell...")
    if not gdf_dicots.empty:
        joined_d = gpd.sjoin(gdf_dicots, grid, how="inner",
                             predicate="intersects")
        grid['DicotCount'] = (grid['grid_id']
                              .map(joined_d.groupby("grid_id").size())
                              .fillna(0).astype(int))
    else:
        grid['DicotCount'] = 0

    print("Counting monocot detections per cell...")
    if not gdf_monocots.empty:
        joined_m = gpd.sjoin(gdf_monocots, grid,
                             how="inner", predicate="intersects")
        grid['MonoCount'] = (grid['grid_id']
                             .map(joined_m.groupby("grid_id").size())
                             .fillna(0).astype(int))
    else:
        grid['MonoCount'] = 0

    # --- Derive summary columns ---
    grid['TotalCount'] = grid['DicotCount'] + grid['MonoCount']
    grid['Inflevel'] = grid.apply(classify_infestation, axis=1)

    # --- Summary statistics ---
    n_weed_cells = (grid['TotalCount'] > 0).sum()
    pct_infested = n_weed_cells / len(grid) * 100
    print(f"\nGrid summary:")
    print(f"  Total cells:    {len(grid)}")
    print(f"  Infested cells: {n_weed_cells} ({pct_infested:.1f}%)")
    print(f"  Mean DicotCount:  {grid['DicotCount'].mean():.3f}")
    print(f"  Mean MonoCount:   {grid['MonoCount'].mean():.3f}")

    # --- Export ---
    output_file = DATA_FOLDER / f"grid_{CELL_SIZE}m_weed_counts.gpkg"
    print(f"\nSaving to {output_file}...")
    grid.to_file(output_file, driver="GPKG")
    print("Done. Load the output file into QGIS as input to 02_bioeconomic_simulation.py")
