# RGI7-like glacier outlines for Iceland

Creates individual ice-flow basin outlines for all Icelandic glaciers, following the RGI7 attribute schema and naming conventions.

## Output files

| File | Basins | Outline year | rgi_id prefix |
|------|--------|-------------|---------------|
| `RGI2010-v7.0-G-06_iceland.shp` | 499 | 2014 | `RGI2010-v7.0-G-06-XXXXX` |
| `RGI2025-v7.0-G-06_iceland.shp` | 572 | 2025 | `RGI2025-v7.0-G-06-XXXXX` |

CRS: EPSG:4326 (WGS84). Geometry type: Polygon.

## Scripts

- `create_rgi7_like_iceland.py` — generates **2014** outlines
- `create_rgi7_like_iceland_2025.py` — generates **2025** outlines
- `validate_rgi7like.py` — RGI7-style validation (geometry, overlaps, schema)

## Processing flow

1. **Load data** — glacier outlines (EPSG:3057), RGI7 Iceland polygons, local ice divide folders
2. **Match ice divides to outlines** — 6 major ice caps matched by centroid containment:
   - Drangajökull, Eyjafjallajökull, Hofsjökull, Langjökull, Mýrdalsjökull, Vatnajökull
3. **Process each outline** — three cases:
   - **Local ice divides available**: convert LineString rings → Polygons, clip to outline, assign uncovered areas to nearest basin
   - **Multiple RGI7 polygons overlap**: use RGI7 as fallback ice divides, clip to outline
   - **Single/no RGI7 match**: use outline directly as a single basin
4. **Per-glacier overlap cleaning** — multi-pass pairwise overlap removal (subtract from smaller basin)
5. **Reproject to EPSG:4326** — with post-reprojection validity fix (`make_valid`)
6. **Global overlap resolution** — RGI7 method: spatial index to find all intersecting pairs, resolve by differencing from the polygon that minimizes total boundary length, remove fully-contained duplicates
7. **Filter tiny basins** — remove basins < 0.01 km² (10,000 m²)
8. **Renumber rgi_id** — sequential numbering after filtering
9. **Write output** — Shapefile with 28 RGI7 attribute columns + `src_file`

## What is populated

- `rgi_id` — sequential counter
- `o1region` — `06` (Iceland)
- `o2region` — `06-01`
- `glac_name` — from ice divide filename or RGI7
- `src_date` — outline year (`2014-01-01` or `2025-01-01`)
- `area_km2` — computed from geometry in EPSG:3057
- `cenlon`, `cenlat` — centroid in WGS84
- `src_file` — source tracking (ice divide filename, `RGI7:<id>`, or `outline_only`)

## What is NOT populated

- `glims_id`, `anlys_id`, `subm_id` — left as null (no GLIMS submission)
- `utm_zone` — not computed
- `termlon`, `termlat` — terminus coordinates not computed
- `zmin_m`, `zmax_m`, `zmed_m`, `zmean_m` — requires DEM intersection
- `slope_deg`, `aspect_deg`, `aspect_sec` — requires DEM intersection
- `dem_source`, `lmax_m` — not computed
- `primeclass`, `conn_lvl`, `surge_type`, `term_type`, `is_rgi6` — left as default (0)

## Input data

- **Glacier outlines**: `\\lv.is\lv\joklar\gögn\útlínur\utlinur-polygons-1890-2025\`
- **Local ice divides**: `\\lv.is\lv\joklar\gögn\vatna-ísaskil\20250701-vatna-isaskil-JHI\`
- **RGI7 reference**: `\\lv.is\lv\joklar\gögn\vatna-ísaskil\RGI2000-v7.0-G-06_iceland\`

## Validation

Run `python validate_rgi7like.py <shapefile>` to check:
- Geometry types (all Polygon), validity, no empty geometries
- Unique rgi_id values
- Self-overlaps using RGI7 spatial index method (tiny boundary artifacts may remain)
- Area statistics
- RGI7 schema completeness (28 required columns)
