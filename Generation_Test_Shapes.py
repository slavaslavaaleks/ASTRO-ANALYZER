import os
import sys
import csv
import copy
import json
import math
import random
import shutil
import argparse
import threading
from collections import deque
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Dict, Optional

import numpy as np
import tifffile

try:
    from scipy.ndimage import gaussian_filter
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
    _HAS_TK = True
except Exception:
    _HAS_TK = False


@dataclass
class CellSpec:
    cell_id: int = 0
    n_soma_vox: int = 0
    n_process_vox: int = 0
    soma_volume_um3: float = 0.0
    process_volume_um3: float = 0.0
    total_volume_um3: float = 0.0
    n_branches: int = 0
    centroid_z: float = 0.0
    centroid_y: float = 0.0
    centroid_x: float = 0.0
    bbox_z: int = 0
    bbox_y: int = 0
    bbox_x: int = 0
    nuc_volume_um3: float = 0.0
    nuc_centroid_z: float = 0.0
    nuc_centroid_y: float = 0.0
    nuc_centroid_x: float = 0.0
    seed_z: int = 0
    seed_y: int = 0
    seed_x: int = 0


@dataclass
class GeneratorConfig:
    z_dim: int = 20
    y_dim: int = 2048
    x_dim: int = 2048
    xy_um: float = 0.144
    z_um: float = 1.0

    n_images: int = 5
    n_cells_per_image: int = 20
    seed: int = 42

    cell_volume_min_um3: float = 100.0
    cell_volume_max_um3: float = 1500.0

    soma_fraction_min: float = 0.35
    soma_fraction_max: float = 0.60

    n_branches_min: int = 3
    n_branches_max: int = 6
    branch_length_um_min: float = 5.0
    branch_length_um_max: float = 12.0
    branch_radius_um_min: float = 0.6
    branch_radius_um_max: float = 1.2

    z_step_penalty: float = 1.2

    nucleus_fraction: float = 0.30

    cell_kind: str = "astrocyte"

    background_level: float = 0.0
    background_gradient_xy: float = 0.0
    background_gradient_z: float = 0.0

    noise_sigma: float = 0.0
    intensity_gradient: float = 0.0
    blur_sigma_um: float = 0.0
    vignette_strength: float = 0.0
    vignette_center_scale: float = 1.0

    normalize_output: bool = False

    _downsample_factor: int = 1


STRESS_LEVELS = {"none": 0.00, "weak": 0.25, "medium": 0.60, "strong": 1.00}

STRESS_FACTORS = {
    "noise":      ["none", "weak", "medium", "strong"],
    "psf":        ["none", "weak", "medium", "strong"],
    "gradient":   ["none", "weak", "medium", "strong"],
    "vignette":   ["none", "weak", "medium", "strong"],
    "downsample": ["none", "weak", "medium", "strong"],
}


def _interp(level_name, min_val, max_val):
    t = STRESS_LEVELS[level_name]
    return min_val + t * (max_val - min_val)


def _vox_to_um3(n_vox, xy_um, z_um):
    return n_vox * (xy_um ** 2) * z_um


def _stamp_ellipsoid(full_shape, center, radii, mask, box):
    z0, z1, y0, y1, x0, x1 = box
    cz, cy, cx = center
    rz, ry, rx = radii

    local_z0 = max(z0, int(math.floor(cz - rz)))
    local_z1 = min(z1, int(math.ceil(cz + rz)) + 1)
    local_y0 = max(y0, int(math.floor(cy - ry)))
    local_y1 = min(y1, int(math.ceil(cy + ry)) + 1)
    local_x0 = max(x0, int(math.floor(cx - rx)))
    local_x1 = min(x1, int(math.ceil(cx + rx)) + 1)

    if local_z1 <= local_z0 or local_y1 <= local_y0 or local_x1 <= local_x0:
        return 0

    zz = (np.arange(local_z0, local_z1) - cz) / max(rz, 1e-6)
    yy = (np.arange(local_y0, local_y1) - cy) / max(ry, 1e-6)
    xx = (np.arange(local_x0, local_x1) - cx) / max(rx, 1e-6)
    m = (zz[:, None, None] ** 2 +
         yy[None, :, None] ** 2 +
         xx[None, None, :] ** 2) <= 1.0

    mask[local_z0 - z0:local_z1 - z0,
    local_y0 - y0:local_y1 - y0,
    local_x0 - x0:local_x1 - x0] |= m
    return int(m.sum())


def _stamp_cylinder(mask, box, start, direction, length, radius):
    z0, z1, y0, y1, x0, x1 = box
    sz, sy, sx = start
    dz, dy, dx = direction
    r = max(1.0, radius)
    r2 = r * r

    n_steps = max(2, int(math.ceil(length * 2)))
    added = 0
    for t in np.linspace(0.0, length, n_steps):
        cz = sz + dz * t
        cy = sy + dy * t
        cx = sx + dx * t
        lz0 = max(z0, int(math.floor(cz - r)))
        lz1 = min(z1, int(math.ceil(cz + r)) + 1)
        ly0 = max(y0, int(math.floor(cy - r)))
        ly1 = min(y1, int(math.ceil(cy + r)) + 1)
        lx0 = max(x0, int(math.floor(cx - r)))
        lx1 = min(x1, int(math.ceil(cx + r)) + 1)
        if lz1 <= lz0 or ly1 <= ly0 or lx1 <= lx0:
            continue
        zz = (np.arange(lz0, lz1) - cz) ** 2
        yy = (np.arange(ly0, ly1) - cy) ** 2
        xx = (np.arange(lx0, lx1) - cx) ** 2
        m = (zz[:, None, None] + yy[None, :, None] + xx[None, None, :]) <= r2
        mask[lz0 - z0:lz1 - z0,
        ly0 - y0:ly1 - y0,
        lx0 - x0:lx1 - x0] |= m
        added += int(m.sum())
    return added


def _pick_radii(target_vox, z_penalty):
    factor = (4.0 / 3.0) * math.pi / max(z_penalty, 1.0)
    r = max(2.0, (target_vox / max(factor, 1e-6)) ** (1.0 / 3.0))
    rz = max(1.5, r / max(z_penalty, 1.0))
    return rz, r, r


def _grid_centers(shape, n_cells, rng, min_sep_yx, min_sep_z):
    z_dim, y_dim, x_dim = shape
    nz = max(1, int(math.ceil(z_dim / max(min_sep_z, 1))))
    ny = max(1, int(math.ceil(y_dim / max(min_sep_yx, 1))))
    nx = max(1, int(math.ceil(x_dim / max(min_sep_yx, 1))))
    while nz * ny * nx < n_cells:
        if nx <= ny and nx <= nz:
            nx += 1
        elif ny <= nz:
            ny += 1
        else:
            nz += 1

    cell_d = z_dim / nz
    cell_h = y_dim / ny
    cell_w = x_dim / nx

    centers = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                cz = (iz + 0.5) * cell_d + rng.uniform(-cell_d * 0.2, cell_d * 0.2)
                cy = (iy + 0.5) * cell_h + rng.uniform(-cell_h * 0.2, cell_h * 0.2)
                cx = (ix + 0.5) * cell_w + rng.uniform(-cell_w * 0.2, cell_w * 0.2)
                cz = int(max(2, min(z_dim - 3, cz)))
                cy = int(max(2, min(y_dim - 3, cy)))
                cx = int(max(2, min(x_dim - 3, cx)))
                centers.append((cz, cy, cx))

    rng.shuffle(centers)
    return centers[:n_cells]


def make_cell_local(shape, cfg, cell_id, rng, center):
    z_dim, y_dim, x_dim = shape
    vox_vol = (cfg.xy_um ** 2) * cfg.z_um

    cell_vol_um3 = rng.uniform(cfg.cell_volume_min_um3,
                               cfg.cell_volume_max_um3)
    total_vox = max(50, int(round(cell_vol_um3 / vox_vol)))

    soma_frac = rng.uniform(cfg.soma_fraction_min, cfg.soma_fraction_max)
    soma_vox = max(30, int(total_vox * soma_frac))
    process_vox = max(0, total_vox - soma_vox)

    if cfg.cell_kind == "astrocyte":
        n_branches = rng.randint(cfg.n_branches_min, cfg.n_branches_max)
    else:
        n_branches = rng.randint(max(2, cfg.n_branches_min - 1),
                                 max(3, cfg.n_branches_max - 2))

    rz_soma, ry_soma, rx_soma = _pick_radii(soma_vox, cfg.z_step_penalty)

    max_len_vox = int(cfg.branch_length_um_max / cfg.xy_um)
    min_len_vox = int(cfg.branch_length_um_min / cfg.xy_um)

    reach = max_len_vox + int(max(rz_soma, ry_soma, rx_soma)) + 10
    cz, cy, cx = center
    z0 = 0
    z1 = z_dim
    y0 = max(0, cy - reach)
    y1 = min(y_dim, cy + reach + 1)
    x0 = max(0, cx - reach)
    x1 = min(x_dim, cx + reach + 1)

    box = (z0, z1, y0, y1, x0, x1)
    lz = z1 - z0
    ly = y1 - y0
    lx = x1 - x0

    soma_mask = np.zeros((lz, ly, lx), dtype=bool)
    process_mask = np.zeros((lz, ly, lx), dtype=bool)

    n_soma = _stamp_ellipsoid(
        shape,
        center=(cz - z0, cy - y0, cx - x0),
        radii=(rz_soma, ry_soma, rx_soma),
        mask=soma_mask,
        box=(0, lz, 0, ly, 0, lx),
    )
    if n_soma == 0:
        return None

    placed_process_vox = 0
    n_branches_done = 0

    for br in range(n_branches):
        remaining = process_vox - placed_process_vox
        if remaining < 30:
            break

        len_um = rng.uniform(cfg.branch_length_um_min,
                             cfg.branch_length_um_max)
        len_vox = max(min_len_vox, int(len_um / cfg.xy_um))

        r_um = rng.uniform(cfg.branch_radius_um_min,
                           cfg.branch_radius_um_max)
        r_vox = max(2.0, r_um / cfg.xy_um)

        branch_vol_estimate = math.pi * r_vox ** 2 * len_vox
        budget = remaining * 1.2
        if branch_vol_estimate > budget:
            r_vox = max(2.0, math.sqrt(budget / max(math.pi * len_vox, 1.0)))

        if r_vox < 1.5:
            continue

        angle = rng.uniform(0, 2 * math.pi)
        elev = rng.uniform(-0.35, 0.35)
        direction = (
            math.sin(elev),
            math.cos(angle) * math.cos(elev),
            math.sin(angle) * math.cos(elev),
        )

        start = (cz - z0 + direction[0] * rz_soma * 0.8,
                 cy - y0 + direction[1] * ry_soma * 0.8,
                 cx - x0 + direction[2] * rx_soma * 0.8)

        added = _stamp_cylinder(
            process_mask, (0, lz, 0, ly, 0, lx),
            start, direction, length=len_vox, radius=r_vox)

        process_mask &= (~soma_mask)
        if added < math.pi * r_vox ** 2 * 3:
            continue

        placed_process_vox = int(process_mask.sum())
        n_branches_done += 1

    all_mask = soma_mask | process_mask
    real_vox = int(all_mask.sum())
    if real_vox < 0.30 * total_vox:
        return None

    coords = np.argwhere(all_mask)
    if len(coords) == 0:
        return None

    cz_l, cy_l, cx_l = coords.mean(axis=0)
    centroid = (cz_l + z0, cy_l + y0, cx_l + x0)

    zs = coords[:, 0]; ys = coords[:, 1]; xs = coords[:, 2]
    bbox_z = int(zs.max() - zs.min() + 1)
    bbox_y = int(ys.max() - ys.min() + 1)
    bbox_x = int(xs.max() - xs.min() + 1)

    n_soma = int(soma_mask.sum())
    n_proc = int(process_mask.sum())

    nuc_radii = (rz_soma * 0.55, ry_soma * 0.55, rx_soma * 0.55)
    nuc_mask = np.zeros_like(soma_mask)
    _stamp_ellipsoid(
        shape,
        center=(cz - z0, cy - y0, cx - x0),
        radii=nuc_radii,
        mask=nuc_mask,
        box=(0, lz, 0, ly, 0, lx),
    )
    nuc_mask &= soma_mask
    nuc_vox = int(nuc_mask.sum())
    if nuc_vox > 0:
        nc = np.argwhere(nuc_mask).mean(axis=0)
        nuc_centroid = (nc[0] + z0, nc[1] + y0, nc[2] + x0)
    else:
        nuc_centroid = (0.0, 0.0, 0.0)

    spec = CellSpec(
        cell_id=cell_id,
        n_soma_vox=n_soma,
        n_process_vox=n_proc,
        soma_volume_um3=_vox_to_um3(n_soma, cfg.xy_um, cfg.z_um),
        process_volume_um3=_vox_to_um3(n_proc, cfg.xy_um, cfg.z_um),
        total_volume_um3=_vox_to_um3(real_vox, cfg.xy_um, cfg.z_um),
        n_branches=n_branches_done,
        centroid_z=float(centroid[0]),
        centroid_y=float(centroid[1]),
        centroid_x=float(centroid[2]),
        bbox_z=bbox_z, bbox_y=bbox_y, bbox_x=bbox_x,
        nuc_volume_um3=_vox_to_um3(nuc_vox, cfg.xy_um, cfg.z_um),
        nuc_centroid_z=float(nuc_centroid[0]),
        nuc_centroid_y=float(nuc_centroid[1]),
        nuc_centroid_x=float(nuc_centroid[2]),
        seed_z=int(cz), seed_y=int(cy), seed_x=int(cx),
    )
    spec._box = box
    spec._soma_mask = soma_mask
    spec._process_mask = process_mask
    spec._nuc_mask = nuc_mask
    return spec

def _paste_into_stack(spec, cells_stack, dapi_stack, intensity_soma, intensity_proc):
    z0, z1, y0, y1, x0, x1 = spec._box
    sub = cells_stack[z0:z1, y0:y1, x0:x1]
    sub_soma = spec._soma_mask
    sub_proc = spec._process_mask
    sub[sub_soma] = np.maximum(sub[sub_soma], intensity_soma)
    sub[sub_proc] = np.maximum(sub[sub_proc], intensity_proc)

    sub_d = dapi_stack[z0:z1, y0:y1, x0:x1]
    sub_nuc = spec._nuc_mask
    sub_d[sub_nuc] = np.maximum(sub_d[sub_nuc], 1.0)


def _add_background(cells, dapi, cfg, rng):
    z_dim, y_dim, x_dim = cells.shape

    base = float(cfg.background_level)
    gxy = float(cfg.background_gradient_xy)
    gz = float(cfg.background_gradient_z)

    bg = np.full((z_dim, y_dim, x_dim), base, dtype=np.float32)

    if gxy > 0:
        yy = np.linspace(-1.0, 1.0, y_dim, dtype=np.float32)
        xx = np.linspace(-1.0, 1.0, x_dim, dtype=np.float32)
        radial = np.sqrt(yy[:, None] ** 2 + xx[None, :] ** 2)
        radial = np.clip(radial / math.sqrt(2.0), 0.0, 1.0)
        radial = 1.0 - gxy * radial
        bg = bg * radial[None, :, :].astype(np.float32)

    if gz > 0:
        zz = np.linspace(1.0 - gz, 1.0, z_dim, dtype=np.float32)
        bg = bg * zz[:, None, None]

    cells = cells + bg
    dapi = dapi + bg * 0.5

    return cells, dapi

def _generate_one_image(cfg, img_idx):
    rng = random.Random(cfg.seed * 1000003 + img_idx)
    rng_np = np.random.default_rng(cfg.seed * 1000003 + img_idx)

    shape = (cfg.z_dim, cfg.y_dim, cfg.x_dim)
    cells_stack = np.zeros(shape, dtype=np.float32)
    dapi_stack = np.zeros(shape, dtype=np.float32)

    max_diam_xy_vox = max(30, int(
        (cfg.cell_volume_max_um3 / ((cfg.xy_um ** 2) * cfg.z_um)) ** (1 / 3)
        * (4 / 3 * math.pi) ** (1 / 3) * 2.0
    ))
    min_sep_yx = max_diam_xy_vox + 4
    min_sep_z = max(3, cfg.z_dim // 3)

    centers = _grid_centers(shape, cfg.n_cells_per_image, rng,
                            min_sep_yx, min_sep_z)
    specs = []
    cell_id = 0
    for center in centers:
        spec = make_cell_local(shape, cfg, cell_id + 1, rng, center)
        if spec is None:
            continue
        cell_id += 1
        spec.cell_id = cell_id
        specs.append(spec)
        soma_int = rng.uniform(0.85, 1.0)
        proc_int = rng.uniform(0.60, 0.80)
        _paste_into_stack(spec, cells_stack, dapi_stack, soma_int, proc_int)

    cells_stack, dapi_stack = _add_background(cells_stack, dapi_stack,
                                              cfg, rng_np)

    return cells_stack, dapi_stack, specs


def save_tiff(path, cells, dapi, cfg):
    def to_u16(a):
        if cfg.normalize_output:
            mx = float(np.max(a))
            if mx <= 0:
                return a.astype(np.uint16)
            return (a / mx * 65535.0).astype(np.uint16)
        return np.clip(a * 65535.0, 0, 65535).astype(np.uint16)

    stack = np.stack([to_u16(cells), to_u16(dapi)], axis=0)
    stack = np.transpose(stack, (1, 0, 2, 3))

    with tifffile.TiffWriter(path, imagej=True) as tw:
        tw.write(
            stack,
            metadata={
                "axes": "ZCYX",
                "unit": "um",
                "spacing": float(cfg.z_um),
                "channels": 2,
            },
            resolution=(1.0 / cfg.xy_um, 1.0 / cfg.xy_um),
            resolutionunit="MICROMETER",
        )


def _write_ground_truth(path, rows):
    if not rows:
        return
    cols = [
        "file", "cell_id",
        "n_soma_vox", "n_process_vox",
        "soma_volume_um3", "process_volume_um3", "total_volume_um3",
        "nuc_volume_um3", "n_branches",
        "centroid_z", "centroid_y", "centroid_x",
        "nuc_centroid_z", "nuc_centroid_y", "nuc_centroid_x",
        "bbox_z", "bbox_y", "bbox_x",
        "seed_z", "seed_y", "seed_x",
        "seed_visible",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def generate_base_bundle(out_dir, cfg, verbose=True, progress_cb=None):
    os.makedirs(out_dir, exist_ok=True)
    all_rows = []

    for img_idx in range(cfg.n_images):
        cells_stack, dapi_stack, specs = _generate_one_image(cfg, img_idx)
        fname = f"stack_{img_idx:03d}.tif"
        save_tiff(os.path.join(out_dir, fname), cells_stack, dapi_stack, cfg)

        for spec in specs:
            d = {k: v for k, v in spec.__dict__.items() if not k.startswith("_")}
            d["file"] = fname
            d["seed_visible"] = True
            all_rows.append(d)

        msg = f"[base {img_idx + 1}/{cfg.n_images}] {fname}: {len(specs)} клеток"
        if verbose:
            print(msg)
        if progress_cb:
            progress_cb(msg)

    _write_ground_truth(os.path.join(out_dir, "ground_truth.csv"), all_rows)
    with open(os.path.join(out_dir, "base_meta.json"), "w") as f:
        json.dump({"cfg": asdict(cfg), "seed": cfg.seed,
                   "shape": [cfg.z_dim, cfg.y_dim, cfg.x_dim],
                   "n_images": cfg.n_images}, f, indent=2)
    return os.path.join(out_dir, "ground_truth.csv"), len(all_rows)


def _apply_noise(cells, dapi, sigma, seed):
    if sigma <= 0:
        return cells, dapi
    rng = np.random.default_rng(seed)
    cells = cells + rng.normal(0, sigma, size=cells.shape).astype(np.float32)
    dapi = dapi + rng.normal(0, sigma, size=dapi.shape).astype(np.float32)
    cells = np.clip(cells, 0, None)
    dapi = np.clip(dapi, 0, None)
    return cells, dapi


def _apply_gradient(cells, dapi, strength):
    if strength <= 0:
        return cells, dapi
    z_dim = cells.shape[0]
    g = np.linspace(1.0 - strength, 1.0, z_dim, dtype=np.float32)
    g = g.reshape(-1, 1, 1)
    return cells * g, dapi * g


def _apply_vignette(cells, dapi, cfg_run, level_t):
    z_dim, y_dim, x_dim = cells.shape

    cy = (y_dim - 1) / 2.0
    cx = (x_dim - 1) / 2.0
    yy = np.arange(y_dim, dtype=np.float32)[:, None] - cy
    xx = np.arange(x_dim, dtype=np.float32)[None, :] - cx
    r2 = yy ** 2 + xx ** 2
    r_max2 = (max(cy, cx) ** 2) or 1.0
    r_norm = np.clip(r2 / r_max2, 0.0, 1.0)

    edge_scale = 1.0 - 0.9 * float(level_t)
    v_mult = 1.0 - (1.0 - edge_scale) * r_norm
    v_mult = np.clip(v_mult, 0.05, 1.0).astype(np.float32)
    v_mult3 = v_mult[None, :, :]

    add_center = float(level_t) * 0.35
    v_add = add_center * (1.0 - r_norm)
    v_add = np.clip(v_add, 0.0, 1.0).astype(np.float32)
    v_add3 = v_add[None, :, :]

    cells = cells * v_mult3 + v_add3
    dapi = dapi * v_mult3 + v_add3 * 0.5

    return cells, dapi


def _apply_psf(cells, dapi, cfg_run, level_t, base_na=0.95, base_lambda_um=0.52):
    if level_t <= 0:
        return cells, dapi
    if not _HAS_SCIPY:
        return cells, dapi

    na_eff = base_na * (1.0 - 0.7 * level_t)
    na_eff = max(na_eff, 0.15)
    lambda_um = base_lambda_um

    airy_radius_um = 0.61 * lambda_um / na_eff
    sigma_xy_um = airy_radius_um / 2.0
    sigma_z_um = airy_radius_um / 1.0

    sigma_xy_vox = sigma_xy_um / cfg_run.xy_um
    sigma_z_vox = sigma_z_um / cfg_run.z_um

    sigma_xy_vox = min(sigma_xy_vox, 12.0)
    sigma_z_vox = min(sigma_z_vox, 6.0)

    cells = gaussian_filter(cells, sigma=(sigma_z_vox, sigma_xy_vox, sigma_xy_vox))
    dapi = gaussian_filter(dapi, sigma=(sigma_z_vox, sigma_xy_vox, sigma_xy_vox))

    return cells, dapi


def _apply_downsample(cells, dapi, factor):
    f = max(1, int(factor))
    if f <= 1:
        return cells, dapi
    cells = cells[::f, ::f, ::f].copy()
    dapi = dapi[::f, ::f, ::f].copy()
    return cells, dapi


def _cfg_for_level(cfg, kind, level):
    c = copy.deepcopy(cfg)
    t = STRESS_LEVELS[level]

    if kind == "noise":
        c.noise_sigma = _interp(level, 0.0, 0.35)
    elif kind == "psf":
        c._psf_level_t = t
    elif kind == "gradient":
        c.intensity_gradient = _interp(level, 0.0, 0.90)
    elif kind == "vignette":
        c.vignette_strength = t
    elif kind == "downsample":
        c._downsample_factor = 1 + int(round(t * 3))
    return c


def _apply_one_stress(cells, dapi, cfg_run, kind, level, seed):
    if kind in ("baseline", "none"):
        return cells, dapi

    if kind == "noise":
        return _apply_noise(cells, dapi, cfg_run.noise_sigma, seed)

    if kind == "psf":
        t = getattr(cfg_run, "_psf_level_t", 0.0)
        return _apply_psf(cells, dapi, cfg_run, t)

    if kind == "gradient":
        return _apply_gradient(cells, dapi, cfg_run.intensity_gradient)

    if kind == "vignette":
        return _apply_vignette(cells, dapi, cfg_run, cfg_run.vignette_strength)

    if kind == "downsample":
        return _apply_downsample(cells, dapi, cfg_run._downsample_factor)

    return cells, dapi


def apply_stress_bundle(base_dir, out_dir, cfg, stress_kind, level_name,
                        verbose=True, progress_cb=None):
    os.makedirs(out_dir, exist_ok=True)
    seed = abs(hash((stress_kind, level_name, cfg.seed))) % (1 << 30)

    for img_idx in range(cfg.n_images):
        cells, dapi, _ = _generate_one_image(cfg, img_idx)

        cfg_run = _cfg_for_level(cfg, stress_kind, level_name)
        if stress_kind == "downsample":
            f = max(1, cfg_run._downsample_factor)
            cfg_run.xy_um = cfg.xy_um * f
            cfg_run.z_um = cfg.z_um * f
            cfg_run.z_dim = max(1, cfg.z_dim // f)

        cells, dapi = _apply_one_stress(cells, dapi, cfg_run,
                                        stress_kind, level_name,
                                        seed + img_idx)

        fname = f"stack_{img_idx:03d}.tif"
        save_tiff(os.path.join(out_dir, fname), cells, dapi, cfg_run)

        msg = f"  [{stress_kind}/{level_name}] {fname}"
        if verbose:
            print(msg)
        if progress_cb:
            progress_cb(msg)

    shutil.copy(os.path.join(base_dir, "ground_truth.csv"),
                os.path.join(out_dir, "ground_truth.csv"))
    with open(os.path.join(out_dir, "stress_meta.json"), "w") as f:
        json.dump({
            "stress_kind": stress_kind,
            "level": level_name,
            "level_value": STRESS_LEVELS[level_name],
        }, f, indent=2)


def build_all_stress_tests(root_dir, base_cfg, verbose=True, progress_cb=None):
    os.makedirs(root_dir, exist_ok=True)
    base_dir = os.path.join(root_dir, "baseline")
    if progress_cb:
        progress_cb("=== Генерация baseline ===")
    generate_base_bundle(base_dir, base_cfg, verbose=verbose,
                         progress_cb=progress_cb)

    for factor, levels in STRESS_FACTORS.items():
        for level in levels:
            out_dir = os.path.join(root_dir, factor, level)
            if level == "none":
                if os.path.exists(out_dir):
                    shutil.rmtree(out_dir)
                os.makedirs(out_dir, exist_ok=True)
                for f in os.listdir(base_dir):
                    if f.endswith(".tif") or f.endswith(".csv"):
                        shutil.copy(os.path.join(base_dir, f),
                                    os.path.join(out_dir, f))
                with open(os.path.join(out_dir, "stress_meta.json"), "w") as f:
                    json.dump({"stress_kind": factor, "level": "none",
                               "note": "copied from baseline"}, f, indent=2)
                if progress_cb:
                    progress_cb(f"=== {factor}/none → копия baseline ===")
                continue

            if progress_cb:
                progress_cb(f"=== {factor}/{level} ===")
            apply_stress_bundle(base_dir, out_dir, base_cfg,
                                stress_kind=factor, level_name=level,
                                verbose=verbose, progress_cb=progress_cb)



def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("out_dir")
    p.add_argument("--build-stress-suite", action="store_true")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--z-dim", type=int, default=20)
    p.add_argument("--y-dim", type=int, default=2048)
    p.add_argument("--x-dim", type=int, default=2048)
    p.add_argument("--xy-um", type=float, default=0.144)
    p.add_argument("--z-um", type=float, default=1.0)
    p.add_argument("--n-images", type=int, default=5)
    p.add_argument("--n-cells", type=int, default=20)
    p.add_argument("--vol-min-um3", type=float, default=100.0)
    p.add_argument("--vol-max-um3", type=float, default=1500.0)
    p.add_argument("--bg", type=float, default=0.0,
                   help="Уровень фона (0..1)")
    p.add_argument("--bg-xy", type=float, default=0.0,
                   help="XY-градиент фона (0..1)")
    p.add_argument("--bg-z", type=float, default=0.0,
                   help="Z-градиент фона (0..1)")
    p.add_argument("--kind", default="astrocyte",
                   choices=["astrocyte", "microglia"])
    p.add_argument("--normalize", action="store_true")
    return p.parse_args()


def cli_main():
    args = parse_args()
    cfg = GeneratorConfig(
        seed=args.seed,
        z_dim=args.z_dim, y_dim=args.y_dim, x_dim=args.x_dim,
        xy_um=args.xy_um, z_um=args.z_um,
        n_images=args.n_images,
        n_cells_per_image=args.n_cells,
        cell_volume_min_um3=args.vol_min_um3,
        cell_volume_max_um3=args.vol_max_um3,
        cell_kind=args.kind,
        background_level=args.bg,
        background_gradient_xy=args.bg_xy,
        background_gradient_z=args.bg_z,
        normalize_output=args.normalize,
    )
    if args.build_stress_suite:
        build_all_stress_tests(args.out_dir, cfg, verbose=True)
        return
    gt, n = generate_base_bundle(args.out_dir, cfg, verbose=True)
    print(f"\nГотово. Строк: {n}")
    print(f"ground_truth: {gt}")



if _HAS_TK:
    class MakerApp:
        BG = "#23252b"
        BG_SOFT = "#2b2e36"
        BG_CARD = "#33363f"
        FG = "#d4d4dc"
        FG_DIM = "#8e93a0"
        ACCENT = "#7aa2f7"
        ACCENT_D = "#5d84d6"
        BORDER = "#3d4048"

        def __init__(self):
            self.root = tk.Tk()
            self.root.title("VoxelMaker v0.8 — фон + PSF-стресс")
            self.root.configure(bg=self.BG)
            self.root.geometry("840x800")
            self._busy = False
            self._build_ui()

        def _row(self, parent, label, default, width=10):
            fr = tk.Frame(parent, bg=self.BG_CARD)
            fr.pack(fill="x", pady=3)
            tk.Label(fr, text=label, width=40, anchor="w",
                     bg=self.BG_CARD, fg=self.FG,
                     font=("Segoe UI", 10)).pack(side="left")
            var = tk.StringVar(value=str(default))
            tk.Entry(fr, textvariable=var, width=width,
                     bg=self.BG_SOFT, fg=self.FG,
                     insertbackground=self.ACCENT,
                     relief="flat", bd=0, highlightthickness=1,
                     highlightbackground=self.BORDER,
                     highlightcolor=self.ACCENT,
                     font=("Consolas", 10)).pack(side="left", fill="x",
                                                 expand=True)
            return var

        def _build_ui(self):
            tk.Label(self.root, text="VoxelMaker v0.8",
                     bg=self.BG, fg=self.FG,
                     font=("Segoe UI", 13, "bold")).pack(
                anchor="w", padx=16, pady=(14, 4))
            tk.Label(self.root,
                     text="Фон + PSF-размытие + градиент + виньетка + даунсемплинг. "
                          "Overlap/dropout удалены.",
                     bg=self.BG, fg=self.FG_DIM,
                     font=("Segoe UI", 9), wraplength=800,
                     justify="left").pack(anchor="w", padx=16, pady=(0, 10))

            card = tk.Frame(self.root, bg=self.BG_CARD, highlightthickness=1,
                            highlightbackground=self.BORDER,
                            highlightcolor=self.BORDER)
            card.pack(fill="both", expand=True, padx=16, pady=(0, 10))
            inner = tk.Frame(card, bg=self.BG_CARD)
            inner.pack(fill="both", expand=True, padx=14, pady=12)

            fr = tk.Frame(inner, bg=self.BG_CARD)
            fr.pack(fill="x", pady=4)
            tk.Label(fr, text="Папка выхода:", width=40, anchor="w",
                     bg=self.BG_CARD, fg=self.FG,
                     font=("Segoe UI", 10)).pack(side="left")
            self.var_dir = tk.StringVar(
                value=os.path.join(os.path.expanduser("~"), "voxel_validation"))
            tk.Entry(fr, textvariable=self.var_dir,
                     bg=self.BG_SOFT, fg=self.FG,
                     insertbackground=self.ACCENT,
                     relief="flat", bd=0, highlightthickness=1,
                     highlightbackground=self.BORDER,
                     highlightcolor=self.ACCENT,
                     font=("Consolas", 10)).pack(side="left", fill="x",
                                                 expand=True)
            tk.Button(fr, text="Обзор…", command=self._choose_dir,
                      bg=self.BG_SOFT, fg=self.FG,
                      activebackground=self.BORDER, activeforeground=self.FG,
                      relief="flat", bd=0, padx=10, pady=3,
                      font=("Segoe UI", 9)).pack(side="left", padx=(6, 0))

            self.var_z = self._row(inner, "Z (срезы):", 20)
            self.var_y = self._row(inner, "Y (px):", 2048)
            self.var_x = self._row(inner, "X (px):", 2048)
            self.var_xy = self._row(inner, "XY (µm/px):", 0.144)
            self.var_zs = self._row(inner, "Z (µm/slice):", 1.0)
            self.var_ni = self._row(inner, "Изображений:", 5)
            self.var_nc = self._row(inner, "Клеток на изображение:", 20)
            self.var_vmin = self._row(inner, "Мин. объём клетки (µm³):", 100.0)
            self.var_vmax = self._row(inner, "Макс. объём клетки (µm³):", 1500.0)
            self.var_bg = self._row(inner, "Фон (0..1):", 0.0)
            self.var_bgxy = self._row(inner, "XY-градиент фона:", 0.0)
            self.var_bgz = self._row(inner, "Z-градиент фона:", 0.0)
            self.var_kind = self._row(inner, "Тип клеток:", "astrocyte")
            self.var_seed = self._row(inner, "Random seed:", 42)

            self.var_norm = tk.BooleanVar(value=False)
            tk.Checkbutton(inner,
                           text="Нормализовать TIFF по максимуму",
                           variable=self.var_norm,
                           bg=self.BG_CARD, fg=self.FG,
                           selectcolor=self.BG_SOFT,
                           activebackground=self.BG_CARD,
                           font=("Segoe UI", 9)).pack(anchor="w", pady=(4, 0))

            btns = tk.Frame(self.root, bg=self.BG)
            btns.pack(fill="x", padx=16, pady=(0, 8))
            self.btn_baseline = tk.Button(
                btns, text="  Baseline  ", command=self._generate_baseline,
                bg=self.BG_CARD, fg=self.FG,
                activebackground=self.BORDER, activeforeground=self.FG,
                relief="flat", bd=0, padx=14, pady=8,
                font=("Segoe UI", 10))
            self.btn_baseline.pack(side="left")

            self.btn_stress = tk.Button(
                btns, text="  Все стресс-тесты  ",
                command=self._generate_stress_suite,
                bg=self.ACCENT, fg="#1a1c22",
                activebackground=self.ACCENT_D, activeforeground="#1a1c22",
                relief="flat", bd=0, padx=16, pady=8,
                font=("Segoe UI", 10, "bold"))
            self.btn_stress.pack(side="left", padx=(8, 0))

            self.btn_exit = tk.Button(
                btns, text="  Выход  ", command=self.root.quit,
                bg=self.BG_CARD, fg=self.FG,
                activebackground=self.BORDER, activeforeground=self.FG,
                relief="flat", bd=0, padx=14, pady=8,
                font=("Segoe UI", 10))
            self.btn_exit.pack(side="right")

            self.lbl_status = tk.Label(self.root, text="",
                                       bg=self.BG, fg=self.FG_DIM,
                                       font=("Consolas", 9), anchor="w",
                                       justify="left", wraplength=800)
            self.lbl_status.pack(fill="x", padx=16, pady=(0, 14))

        def _choose_dir(self):
            d = filedialog.askdirectory(initialdir=self.var_dir.get() or ".")
            if d:
                self.var_dir.set(d)

        def _collect_cfg(self):
            return GeneratorConfig(
                z_dim=int(self.var_z.get()),
                y_dim=int(self.var_y.get()),
                x_dim=int(self.var_x.get()),
                xy_um=float(self.var_xy.get().replace(",", ".")),
                z_um=float(self.var_zs.get().replace(",", ".")),
                n_images=int(self.var_ni.get()),
                n_cells_per_image=int(self.var_nc.get()),
                cell_volume_min_um3=float(self.var_vmin.get().replace(",", ".")),
                cell_volume_max_um3=float(self.var_vmax.get().replace(",", ".")),
                background_level=float(self.var_bg.get().replace(",", ".")),
                background_gradient_xy=float(self.var_bgxy.get().replace(",", ".")),
                background_gradient_z=float(self.var_bgz.get().replace(",", ".")),
                cell_kind=self.var_kind.get().strip(),
                seed=int(self.var_seed.get()),
                normalize_output=bool(self.var_norm.get()),
            )

        def _set_busy(self, busy, text=""):
            self._busy = busy
            state = "disabled" if busy else "normal"
            for b in (self.btn_baseline, self.btn_stress, self.btn_exit):
                b.configure(state=state)
            self.lbl_status.config(text=text)
            self.root.update_idletasks()

        def _progress(self, msg):
            self.root.after(0, lambda: self.lbl_status.config(text=msg))

        def _on_done(self, text):
            self._set_busy(False, "")
            messagebox.showinfo("Готово", text)

        def _on_error(self, text):
            self._set_busy(False, "")
            messagebox.showerror("Ошибка", text)

        def _generate_baseline(self):
            if self._busy:
                return
            out_dir = self.var_dir.get().strip()
            if not out_dir:
                messagebox.showerror("Ошибка", "Укажите папку выхода")
                return
            cfg = self._collect_cfg()
            self._set_busy(True, "Генерация baseline…")

            def worker():
                try:
                    gt, n = generate_base_bundle(out_dir, cfg, verbose=False,
                                                 progress_cb=self._progress)
                    self.root.after(0, lambda: self._on_done(
                        f"Готово.\nКлеток: {n}\n{gt}"))
                except Exception as e:
                    err = str(e)
                    self.root.after(0, lambda: self._on_error(err))

            threading.Thread(target=worker, daemon=True).start()

        def _generate_stress_suite(self):
            if self._busy:
                return
            out_dir = self.var_dir.get().strip()
            if not out_dir:
                messagebox.showerror("Ошибка", "Укажите папку выхода")
                return
            cfg = self._collect_cfg()
            self._set_busy(True, "Генерация набора стресс-тестов…")

            def worker():
                try:
                    build_all_stress_tests(out_dir, cfg, verbose=False,
                                           progress_cb=self._progress)
                    self.root.after(0, lambda: self._on_done(
                        f"Набор стресс-тестов построен:\n{out_dir}"))
                except Exception as e:
                    err = str(e)
                    self.root.after(0, lambda: self._on_error(err))

            threading.Thread(target=worker, daemon=True).start()

        def run(self):
            self.root.mainloop()


def main():
    if len(sys.argv) == 1 and _HAS_TK:
        MakerApp().run()
        return
    cli_main()


if __name__ == "__main__":
    main()