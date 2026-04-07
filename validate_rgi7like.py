"""
Validate RGI7-like shapefiles using the same overlap detection approach
as the official RGI7 scripts (overlaps_helpers.py).

Usage:
    python validate_rgi7like.py <shapefile>
"""

import sys
import geopandas as gpd
import pandas as pd
import numpy as np
import shapely
from shapely.validation import make_valid


# ── From RGI7 overlaps_helpers.py ──────────────────────────────────────────────

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
    return shapely.MultiPolygon(geoms)


def compute_self_overlaps(gs):
    """
    Compute overlaps between pairs of input geometries.
    Returns GeoDataFrame with i, j, i_area_fraction, j_area_fraction, geometry.
    """
    sindex = gs.sindex
    print("  Finding intersecting geometries...")
    matches = sindex.query(gs, "intersects")
    is_unique_pair = (matches[0] != matches[1]) & (matches[0] < matches[1])
    pairs = matches[:, is_unique_pair].transpose()
    print(f"  Found {len(pairs)} intersecting pairs to check")

    overlaps = []
    for counter, (i, j) in enumerate(pairs):
        if counter % 100 == 0 or counter == len(pairs) - 1:
            print(f"  [{len(pairs)}] {counter + 1}", end="\r", flush=True)
        overlap = gs.iloc[i].intersection(gs.iloc[j])
        try:
            overlap = polygonize(overlap)
        except ValueError:
            continue
        overlaps.append({
            "i": i,
            "j": j,
            "i_area_fraction": overlap.area / gs.iloc[i].area,
            "j_area_fraction": overlap.area / gs.iloc[j].area,
            "geometry": overlap,
        })
    print()

    columns = ["i", "j", "i_area_fraction", "j_area_fraction", "geometry"]
    return gpd.GeoDataFrame(overlaps, columns=columns, crs=gs.crs)


# ── Validation ─────────────────────────────────────────────────────────────────

def validate(path):
    print("=" * 70)
    print(f"Validating: {path}")
    print("=" * 70)

    gdf = gpd.read_file(path)
    print(f"\nBasic info:")
    print(f"  Total features: {len(gdf)}")
    print(f"  CRS: {gdf.crs}")

    # 1. Geometry types
    print(f"\n--- Geometry Types ---")
    type_counts = gdf.geometry.geom_type.value_counts()
    for gtype, count in type_counts.items():
        status = "PASS" if gtype == "Polygon" else "WARN"
        print(f"  [{status}] {gtype}: {count}")
    multi = (gdf.geometry.geom_type != 'Polygon').sum()
    if multi == 0:
        print(f"  [PASS] All geometries are Polygon")
    else:
        print(f"  [FAIL] {multi} non-Polygon geometries found")

    # 2. Validity
    print(f"\n--- Geometry Validity ---")
    invalid_mask = ~gdf.geometry.is_valid
    n_invalid = invalid_mask.sum()
    if n_invalid == 0:
        print(f"  [PASS] All {len(gdf)} geometries are valid")
    else:
        print(f"  [FAIL] {n_invalid} invalid geometries")
        for idx in gdf[invalid_mask].index[:5]:
            from shapely.validation import explain_validity
            print(f"    rgi_id={gdf.loc[idx, 'rgi_id']}: {explain_validity(gdf.loc[idx, 'geometry'])}")

    # 3. Empty geometries
    print(f"\n--- Empty Geometries ---")
    empty = gdf.geometry.is_empty.sum()
    if empty == 0:
        print(f"  [PASS] No empty geometries")
    else:
        print(f"  [FAIL] {empty} empty geometries")

    # 4. Duplicate rgi_id
    print(f"\n--- Duplicate rgi_id ---")
    dupes = gdf['rgi_id'].duplicated().sum()
    if dupes == 0:
        print(f"  [PASS] All rgi_id values are unique")
    else:
        print(f"  [FAIL] {dupes} duplicate rgi_id values")

    # 5. Self-overlaps (the main RGI7 check)
    print(f"\n--- Self-Overlaps (RGI7 method) ---")
    gs = gdf.geometry
    overlaps = compute_self_overlaps(gs)
    if len(overlaps) == 0:
        print(f"  [PASS] No overlaps found")
    else:
        print(f"  [WARN] {len(overlaps)} overlapping pairs found")
        # Classify by severity
        tiny = overlaps[(overlaps['i_area_fraction'] < 0.001) & (overlaps['j_area_fraction'] < 0.001)]
        small = overlaps[((overlaps['i_area_fraction'] >= 0.001) | (overlaps['j_area_fraction'] >= 0.001)) &
                         (overlaps['i_area_fraction'] < 0.01) & (overlaps['j_area_fraction'] < 0.01)]
        significant = overlaps[(overlaps['i_area_fraction'] >= 0.01) | (overlaps['j_area_fraction'] >= 0.01)]

        print(f"    Tiny (<0.1% of both):     {len(tiny)}")
        print(f"    Small (0.1%-1%):           {len(small)}")
        print(f"    Significant (>1%):         {len(significant)}")

        if len(significant) > 0:
            print(f"\n  Significant overlaps:")
            for _, row in significant.head(10).iterrows():
                i, j = int(row['i']), int(row['j'])
                print(f"    {gdf.iloc[i]['rgi_id']} & {gdf.iloc[j]['rgi_id']}: "
                      f"{row['i_area_fraction']*100:.2f}% / {row['j_area_fraction']*100:.2f}%")

    # 6. Area check
    print(f"\n--- Area Summary ---")
    if gdf.crs.is_geographic:
        gdf_proj = gdf.to_crs("EPSG:3057")
        areas = gdf_proj.geometry.area / 1e6
    else:
        areas = gdf.geometry.area / 1e6
    print(f"  Total area: {areas.sum():.1f} km²")
    print(f"  Min basin:  {areas.min():.4f} km²")
    print(f"  Max basin:  {areas.max():.1f} km²")
    print(f"  Mean basin: {areas.mean():.2f} km²")
    tiny_basins = (areas < 0.01).sum()
    if tiny_basins > 0:
        print(f"  [WARN] {tiny_basins} basins smaller than 0.01 km²")

    # 7. Required RGI7 columns
    print(f"\n--- RGI7 Schema ---")
    required = ['rgi_id', 'o1region', 'o2region', 'glims_id', 'anlys_id',
                'subm_id', 'src_date', 'cenlon', 'cenlat', 'utm_zone',
                'area_km2', 'primeclass', 'conn_lvl', 'surge_type',
                'term_type', 'glac_name', 'is_rgi6', 'termlon', 'termlat',
                'zmin_m', 'zmax_m', 'zmed_m', 'zmean_m', 'slope_deg',
                'aspect_deg', 'aspect_sec', 'dem_source', 'lmax_m']
    missing = [c for c in required if c not in gdf.columns]
    if missing:
        print(f"  [FAIL] Missing columns: {missing}")
    else:
        print(f"  [PASS] All {len(required)} RGI7 columns present")

    print(f"\n{'=' * 70}")
    all_pass = (multi == 0 and n_invalid == 0 and empty == 0 and
                dupes == 0 and len(overlaps) == 0 and not missing)
    if all_pass:
        print("RESULT: ALL CHECKS PASSED")
    else:
        issues = []
        if multi > 0: issues.append(f"{multi} non-Polygon")
        if n_invalid > 0: issues.append(f"{n_invalid} invalid")
        if empty > 0: issues.append(f"{empty} empty")
        if dupes > 0: issues.append(f"{dupes} duplicate IDs")
        if len(overlaps) > 0: issues.append(f"{len(overlaps)} overlaps")
        if missing: issues.append(f"missing columns: {missing}")
        print(f"RESULT: ISSUES FOUND — {', '.join(issues)}")
    print("=" * 70)
    return all_pass


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python validate_rgi7like.py <shapefile>")
        sys.exit(1)
    validate(sys.argv[1])
