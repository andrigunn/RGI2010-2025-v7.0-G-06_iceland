"""
Create RGI7-like glacier outlines for ALL glaciers in Iceland.

For glaciers with local ice divide data (6 major ice caps), uses those divides.
For all other glaciers, uses RGI7 ice divides adjusted to the local 2025 outlines.

Input:
  - 2025 glacier outlines (385 polygons): the reference extent
  - Local ice divide folders (6 folders with LineString rings)
  - Full Iceland RGI7 data (568 polygons): schema template + fallback divides

Output:
  - Single shapefile with all basins, matching RGI7 attribute schema
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import os
import warnings
from pathlib import Path
import shapely
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.validation import make_valid, explain_validity


# ── Paths ──────────────────────────────────────────────────────────────────────
OUTLINES_PATH = (
    r"\\lv.is\lv\joklar\gögn\útlínur\utlinur-polygons-1890-2025"
    r"\iceland-glaciers-outlines-2025-glacier-boundary-polygon.shp"
)
ICE_DIVIDES_BASE = (
    r"\\lv.is\lv\joklar\gögn\vatna-ísaskil\20250701-vatna-isaskil-JHI"
)
RGI7_PATH = (
    r"\\lv.is\lv\joklar\gögn\vatna-ísaskil"
    r"\RGI2000-v7.0-G-06_iceland\RGI2000-v7.0-G-06_iceland.shp"
)
OUTPUT_DIR = (
    r"\\lv.is\lv\joklar\verkefni\2025 - AMOC-Jöklar-Orka"
    r"\RGI2010-2025-v7.0-G-06_iceland"
)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "RGI2025-v7.0-G-06_iceland.shp")

RGI_ID_BASE = "RGI2025-v7.0-G-06"


# ── Geometry helpers ───────────────────────────────────────────────────────────

def ensure_valid(geom):
    """Make geometry valid if needed."""
    if geom is None or geom.is_empty:
        return geom
    if not geom.is_valid:
        geom = make_valid(geom)
    return geom


def ensure_polygon(geom):
    """Convert MultiPolygon/GeometryCollection to single Polygon (largest)."""
    if geom is None or geom.is_empty:
        return geom
    if geom.geom_type == 'MultiPolygon':
        return max(geom.geoms, key=lambda g: g.area)
    elif geom.geom_type == 'GeometryCollection':
        polys = [g for g in geom.geoms if g.geom_type in ('Polygon', 'MultiPolygon')]
        if polys:
            largest = max(polys, key=lambda g: g.area)
            if largest.geom_type == 'MultiPolygon':
                return max(largest.geoms, key=lambda g: g.area)
            return largest
        return None
    return geom


def line_to_polygon(line_geom):
    """Convert a closed/near-closed LineString to a valid Polygon."""
    poly = Polygon(line_geom.coords)
    if not poly.is_valid:
        poly = make_valid(poly)
    return ensure_polygon(poly)


def clip_and_clean(geom, outline_poly):
    """Clip geometry to outline and ensure it's a valid single Polygon."""
    clipped = ensure_valid(geom).intersection(ensure_valid(outline_poly))
    if clipped.is_empty:
        return None
    clipped = ensure_valid(clipped)
    clipped = ensure_polygon(clipped)
    return clipped


def assign_uncovered_to_nearest(outline_poly, basins):
    """Assign uncovered areas in the outline to the nearest basin."""
    if not basins:
        return basins

    union_basins = unary_union(basins)
    uncovered = outline_poly.difference(union_basins)

    if uncovered.is_empty:
        return basins

    centroids = [b.centroid for b in basins]

    if uncovered.geom_type == 'MultiPolygon':
        parts = list(uncovered.geoms)
    elif uncovered.geom_type == 'Polygon':
        parts = [uncovered]
    elif hasattr(uncovered, 'geoms'):
        parts = [g for g in uncovered.geoms
                 if g.geom_type in ('Polygon', 'MultiPolygon') and not g.is_empty]
    else:
        return basins

    updated = list(basins)
    for part in parts:
        if part.is_empty or part.area < 1:
            continue
        distances = [part.centroid.distance(c) for c in centroids]
        nearest_idx = np.argmin(distances)
        updated[nearest_idx] = unary_union([updated[nearest_idx], part])

    for i, geom in enumerate(updated):
        updated[i] = ensure_polygon(geom)

    return updated


def clean_overlaps(geometries, names):
    """
    Fix self-intersections and remove pairwise overlaps.
    Returns cleaned list of geometries.
    """
    # Fix validity
    for i, geom in enumerate(geometries):
        if geom is not None and not geom.is_valid:
            geometries[i] = ensure_polygon(make_valid(geom))

    # Remove overlaps in multiple passes
    for _pass in range(3):
        found = False
        for i in range(len(geometries)):
            if geometries[i] is None:
                continue
            for j in range(i + 1, len(geometries)):
                if geometries[j] is None:
                    continue
                overlap = geometries[i].intersection(geometries[j])
                if overlap.is_empty or overlap.area <= 0:
                    continue
                found = True
                # Subtract from smaller basin
                if geometries[i].area >= geometries[j].area:
                    fix = geometries[j].difference(geometries[i])
                    geometries[j] = ensure_polygon(ensure_valid(fix))
                else:
                    fix = geometries[i].difference(geometries[j])
                    geometries[i] = ensure_polygon(ensure_valid(fix))
        if not found:
            break

    return geometries


# ── RGI7-style global overlap resolution ───────────────────────────────────────

def polygonize(geom):
    """Convert geometry to (Multi)Polygon, removing all zero-area components."""
    try:
        geoms = [x for x in geom.geoms if x.area]
    except AttributeError:
        if geom.area:
            return geom
        raise ValueError("Geometry has zero area")
    if not geoms:
        raise ValueError("Geometry has zero area")
    if len(geoms) == 1:
        return geoms[0]
    return MultiPolygon(geoms)


def compute_self_overlaps(gs):
    """Compute pairwise overlaps using spatial index (RGI7 method)."""
    sindex = gs.sindex
    matches = sindex.query(gs, "intersects")
    is_unique_pair = (matches[0] != matches[1]) & (matches[0] < matches[1])
    pairs = matches[:, is_unique_pair].transpose()

    overlaps = []
    for i, j in pairs:
        overlap = gs.iloc[i].intersection(gs.iloc[j])
        try:
            overlap = polygonize(overlap)
        except ValueError:
            continue
        overlaps.append({
            "i": i, "j": j,
            "i_area_fraction": overlap.area / gs.iloc[i].area,
            "j_area_fraction": overlap.area / gs.iloc[j].area,
            "geometry": overlap,
        })

    columns = ["i", "j", "i_area_fraction", "j_area_fraction", "geometry"]
    return gpd.GeoDataFrame(overlaps, columns=columns, crs=gs.crs)


def resolve_self_overlaps(overlaps, geoms, min_area=10000):
    """
    Resolve overlaps using the RGI7 method: difference overlap from the
    polygon that minimizes total boundary length.
    Returns dict of {positional_index: fixed_geometry} and set of indices to remove.
    """
    fixed = {}
    to_remove = set()
    for _, row in overlaps.iterrows():
        gi = [row["i"], row["j"]]
        # Skip if already marked for removal
        if gi[0] in to_remove or gi[1] in to_remove:
            continue
        g = [fixed[i] if i in fixed else geoms.iloc[i] for i in gi]

        # Try differencing
        d = [None, None]
        for k in range(2):
            diff = g[k].difference(row["geometry"])
            try:
                d[k] = polygonize(diff)
            except ValueError:
                d[k] = None  # zero-area result

        # If one difference is zero-area, the polygon is fully contained
        # Remove the smaller (contained) polygon
        if d[0] is None and d[1] is None:
            # Both zero - remove the smaller one
            smaller = 0 if g[0].area <= g[1].area else 1
            to_remove.add(gi[smaller])
            continue
        if d[0] is None:
            # g[0] fully inside overlap -> remove g[0]
            to_remove.add(gi[0])
            continue
        if d[1] is None:
            # g[1] fully inside overlap -> remove g[1]
            to_remove.add(gi[1])
            continue

        # Choose difference that minimizes total perimeter
        first_choice = int(
            (d[0].boundary.length + g[1].boundary.length)
            > (g[0].boundary.length + d[1].boundary.length)
        )

        fixed_geom = None
        for choice in (first_choice, int(not first_choice)):
            if isinstance(d[choice], shapely.MultiPolygon):
                parts = [p for p in d[choice].geoms if p.area > min_area]
                if len(parts) == 1:
                    fixed_geom = parts[0]
                    break
            else:
                fixed_geom = d[choice]
                break

        if fixed_geom:
            fixed[gi[choice]] = fixed_geom

    return fixed, to_remove


# ── Processing functions ───────────────────────────────────────────────────────

def process_local_ice_divides(outline_poly, folder_path, outline_crs):
    """
    Process a folder of local ice divide shapefiles (LineString rings).
    Returns (geometries, names, src_files) clipped to the outline.
    """
    shp_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.shp')])
    if not shp_files:
        return [], [], []

    geometries = []
    names = []
    src_files = []

    for shp_file in shp_files:
        gdf = gpd.read_file(os.path.join(folder_path, shp_file))
        if gdf.crs != outline_crs:
            gdf = gdf.to_crs(outline_crs)

        line = gdf.geometry.iloc[0]
        basin_poly = line_to_polygon(line)
        clipped = clip_and_clean(basin_poly, outline_poly)

        if clipped is None or clipped.is_empty:
            continue

        name_part = Path(shp_file).stem.split('-')[-1]
        geometries.append(clipped)
        names.append(name_part)
        src_files.append(shp_file)

    # Assign uncovered areas
    if geometries:
        geometries = assign_uncovered_to_nearest(outline_poly, geometries)

    return geometries, names, src_files


def process_rgi7_ice_divides(outline_poly, rgi7_basins, outline_crs):
    """
    Use RGI7 polygons as ice divides, clipped to the local outline.
    Returns (geometries, names, src_files).
    """
    geometries = []
    names = []
    src_files = []

    for _, rgi_row in rgi7_basins.iterrows():
        rgi_geom = ensure_valid(rgi_row.geometry)
        clipped = clip_and_clean(rgi_geom, outline_poly)

        if clipped is None or clipped.is_empty or clipped.area < 1:
            continue

        name = rgi_row.get('glac_name', None)
        if name is None or (isinstance(name, float) and np.isnan(name)):
            name = rgi_row['rgi_id']

        geometries.append(clipped)
        names.append(name)
        src_files.append(f"RGI7:{rgi_row['rgi_id']}")

    # Assign uncovered areas
    if geometries:
        geometries = assign_uncovered_to_nearest(outline_poly, geometries)

    return geometries, names, src_files


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Creating RGI7-like outlines for all Iceland glaciers")
    print("=" * 70)

    # --- Load data ---
    print("\nLoading 2025 glacier outlines...")
    outlines = gpd.read_file(OUTLINES_PATH)
    outlines['geometry'] = outlines.geometry.apply(ensure_valid)
    outline_crs = outlines.crs
    print(f"  {len(outlines)} glacier outlines loaded")

    print("Loading RGI7 Iceland data...")
    rgi7 = gpd.read_file(RGI7_PATH)
    rgi7_proj = rgi7.to_crs(outline_crs)
    rgi7_proj['geometry'] = rgi7_proj.geometry.apply(ensure_valid)
    rgi7_cols = [c for c in rgi7.columns if c != 'geometry']
    print(f"  {len(rgi7)} RGI7 polygons loaded")

    # --- Identify local ice divide folders ---
    print("\nScanning for local ice divide folders...")
    ice_divide_folders = {}
    for item in os.listdir(ICE_DIVIDES_BASE):
        item_path = os.path.join(ICE_DIVIDES_BASE, item)
        if os.path.isdir(item_path) and '-ísaskil-' in item:
            # Extract glacier name (first part before -ísaskil-)
            glacier_prefix = item.split('-ísaskil-')[0]
            ice_divide_folders[glacier_prefix] = item_path
            n_shps = len([f for f in os.listdir(item_path) if f.endswith('.shp')])
            print(f"  {item}: {n_shps} ice divide files")

    # --- Match ice divide folders to outlines by spatial intersection ---
    print("\nMatching ice divide folders to outlines...")
    outline_to_folder = {}  # outline_idx -> folder_path
    for glacier_name, folder_path in ice_divide_folders.items():
        # Read first shapefile to get centroid location
        shps = [f for f in os.listdir(folder_path) if f.endswith('.shp')]
        if not shps:
            continue
        gdf = gpd.read_file(os.path.join(folder_path, shps[0]))
        if gdf.crs != outline_crs:
            gdf = gdf.to_crs(outline_crs)
        centroid = gdf.geometry.iloc[0].centroid

        # Find which outline contains this centroid
        matched = False
        for idx, row in outlines.iterrows():
            if row.geometry.contains(centroid):
                outline_to_folder[idx] = (glacier_name, folder_path)
                area = row.geometry.area / 1e6
                print(f"  {glacier_name} -> outline idx={idx} ({area:.1f} km²)")
                matched = True
                break
        if not matched:
            # Fallback: nearest outline
            dists = outlines.geometry.distance(centroid)
            nearest = dists.idxmin()
            outline_to_folder[nearest] = (glacier_name, folder_path)
            area = outlines.loc[nearest].geometry.area / 1e6
            print(f"  {glacier_name} -> nearest outline idx={nearest} ({area:.1f} km²)")

    # --- Process each outline ---
    print("\n" + "=" * 70)
    print("Processing all glacier outlines...")
    print("=" * 70)

    all_rows = []
    counter = 0

    for outline_idx, outline_row in outlines.iterrows():
        outline_poly = outline_row.geometry
        outline_area = outline_poly.area / 1e6

        if outline_idx in outline_to_folder:
            # ── Local ice divides available ──
            glacier_name, folder_path = outline_to_folder[outline_idx]
            print(f"\n[{outline_idx}] {glacier_name} ({outline_area:.1f} km²) - LOCAL ice divides")

            geometries, names, src_files = process_local_ice_divides(
                outline_poly, folder_path, outline_crs
            )
        else:
            # ── Use RGI7 as fallback ──
            # Find RGI7 polygons that intersect this outline
            rgi7_intersecting = rgi7_proj[rgi7_proj.intersects(outline_poly)]

            if len(rgi7_intersecting) == 0:
                # No RGI7 match - use the outline as a single basin
                print(f"\n[{outline_idx}] Outline ({outline_area:.2f} km²) - NO RGI7 match, single basin")
                geometries = [outline_poly]
                names = [f"glacier_{outline_idx}"]
                src_files = ["outline_only"]
            elif len(rgi7_intersecting) == 1:
                # Single RGI7 polygon - use outline directly
                rgi_row = rgi7_intersecting.iloc[0]
                rgi_name = rgi_row.get('glac_name', None)
                if rgi_name is None or (isinstance(rgi_name, float) and np.isnan(rgi_name)):
                    rgi_name = rgi_row['rgi_id']
                print(f"\n[{outline_idx}] {rgi_name} ({outline_area:.2f} km²) - single RGI7 basin")
                geometries = [outline_poly]
                names = [rgi_name]
                src_files = [f"RGI7:{rgi_row['rgi_id']}"]
            else:
                # Multiple RGI7 polygons - use them as ice divides
                first_name = rgi7_intersecting.iloc[0].get('glac_name', None)
                if first_name is None or (isinstance(first_name, float) and np.isnan(first_name)):
                    first_name = f"outline_{outline_idx}"
                print(f"\n[{outline_idx}] {first_name} area ({outline_area:.2f} km²) - "
                      f"{len(rgi7_intersecting)} RGI7 basins")

                geometries, names, src_files = process_rgi7_ice_divides(
                    outline_poly, rgi7_intersecting, outline_crs
                )

        if not geometries:
            print(f"  WARNING: No basins created for outline {outline_idx}, skipping")
            continue

        # Clean overlaps within this glacier
        geometries = clean_overlaps(geometries, names)

        # Build rows for this glacier
        for geom, name, src_file in zip(geometries, names, src_files):
            if geom is None or geom.is_empty:
                continue
            counter += 1

            row = {}
            for col in rgi7_cols:
                dtype = rgi7.dtypes[col]
                if dtype == 'object':
                    row[col] = None
                elif dtype == 'int64':
                    row[col] = 0
                elif dtype == 'float64':
                    row[col] = 0.0
                else:
                    row[col] = None

            row['rgi_id'] = f"{RGI_ID_BASE}-{counter:05d}"
            row['o1region'] = '06'
            row['o2region'] = '06-01'
            row['glac_name'] = name
            row['src_date'] = '2025-01-01T00:00:00'
            row['area_km2'] = geom.area / 1e6

            all_rows.append({**row, 'geometry': geom, 'src_file': src_file})

        print(f"  -> {len(geometries)} basins created")

    # --- Build output GeoDataFrame ---
    print(f"\n{'=' * 70}")
    print(f"Building output ({len(all_rows)} total basins)...")

    output_gdf = gpd.GeoDataFrame(all_rows, crs=outline_crs)

    # Reproject to EPSG:4326
    output_gdf = output_gdf.to_crs("EPSG:4326")

    # Final validity fix after reprojection
    output_gdf['geometry'] = output_gdf.geometry.apply(
        lambda g: make_valid(g) if not g.is_valid else g
    )
    output_gdf['geometry'] = output_gdf.geometry.apply(ensure_polygon)

    # --- RGI7-style overlap resolution ---
    print("Resolving global overlaps (RGI7 method)...")
    gs = output_gdf.geometry
    overlaps = compute_self_overlaps(gs)
    if len(overlaps) > 0:
        print(f"  Found {len(overlaps)} overlapping pairs, resolving...")
        fixed, to_remove = resolve_self_overlaps(overlaps, gs)
        for idx, geom in fixed.items():
            output_gdf.at[output_gdf.index[idx], 'geometry'] = geom
        if to_remove:
            print(f"  Removing {len(to_remove)} fully-contained duplicate basins")
            drop_indices = [output_gdf.index[i] for i in to_remove]
            output_gdf = output_gdf.drop(drop_indices).reset_index(drop=True)
        # Verify
        remaining = compute_self_overlaps(output_gdf.geometry)
        print(f"  After resolution: {len(remaining)} overlaps remaining")
    else:
        print("  No overlaps found")

    # --- Filter tiny basins ---
    # Project to compute area in m², then filter
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        areas_m2 = output_gdf.to_crs("EPSG:3057").geometry.area
    tiny_mask = areas_m2 < 10000  # < 0.01 km²
    n_tiny = tiny_mask.sum()
    if n_tiny > 0:
        print(f"Removing {n_tiny} tiny basins (< 0.01 km²)...")
        output_gdf = output_gdf[~tiny_mask].reset_index(drop=True)

    # --- Renumber rgi_id after filtering ---
    for i, idx in enumerate(output_gdf.index):
        output_gdf.at[idx, 'rgi_id'] = f"{RGI_ID_BASE}-{i+1:05d}"

    # Recompute area_km2 after overlap fixes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        output_gdf['area_km2'] = output_gdf.to_crs("EPSG:3057").geometry.area / 1e6

    # Compute centroids in WGS84 (suppress geographic CRS warning - we want lon/lat)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        centroids = output_gdf.geometry.centroid
    output_gdf['cenlon'] = centroids.x
    output_gdf['cenlat'] = centroids.y

    # Column order
    col_order = [c for c in rgi7_cols if c in output_gdf.columns]
    col_order += ['src_file', 'geometry']
    output_gdf = output_gdf[col_order]

    # --- Write output ---
    print(f"Writing to {OUTPUT_FILE}...")
    output_gdf.to_file(OUTPUT_FILE, promote_to_multi=False, geometry_type="Polygon")

    print(f"\nDone! Created {len(output_gdf)} basin outlines.")
    print(f"Output: {OUTPUT_FILE}")

    # Summary
    local_count = output_gdf[~output_gdf.src_file.str.startswith('RGI7:') &
                              (output_gdf.src_file != 'outline_only')].shape[0]
    rgi7_count = output_gdf[output_gdf.src_file.str.startswith('RGI7:')].shape[0]
    single_count = output_gdf[output_gdf.src_file == 'outline_only'].shape[0]
    print(f"\n  Local ice divides: {local_count} basins")
    print(f"  RGI7 ice divides:  {rgi7_count} basins")
    print(f"  Single-basin (no match): {single_count} basins")
    print(f"  Total area: {output_gdf.area_km2.sum():.1f} km²")

    return output_gdf


if __name__ == "__main__":
    main()
