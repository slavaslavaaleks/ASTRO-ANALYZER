# ASTRO-ANALYZER v1.0.0

**ASTRO-ANALYZER** is a program for automatic segmentation and morphometry of astrocytes, neurons, and microglial cells on 3D stacks and 2D images obtained by confocal and epifluorescence microscopy (tissue sections or cell culture).

The output is a table of morphological parameters for each cell: soma volume, nucleus volume, whole-cell volume with processes, projection areas, mean intensities, elongation, nuclear-cytoplasmic ratio, process branching complexity, and many other features.

---

## I. PIPELINE

For each input file (.tif / .nd2 / .czi / .lif / .ome.tif / .oir / .oif stack), the program performs the following steps:

1. Reads the image and automatically extracts metadata (XY and Z scale, number and parameters of channels, optical parameters). Offers the user a channel selection dialog based on automatic metadata analysis.
2. Performs 3D deconvolution (if enabled by the user) to suppress out-of-focus haze.
3. Asks the user to place soma centers on the maximum intensity projection (MIP) or finds them automatically (Beta).
4. Automatically outlines the soma of each labeled cell in 3D.
5. Automatically finds the nucleus of each cell in the DAPI channel (if enabled by the user).
6. Grows the whole cell (soma + processes) in the cell channel.
7. Segments granules within each cell (if enabled by the user).
8. Allows manual correction or removal of problematic cells (if manual filtering is enabled).
9. Computes morphological parameters for each cell and saves them in CSV format.
10. Generates visualizations of each cell's masks.
11. Aggregates results from all files in the given folder into a cumulative CSV (morphometric parameters and cell counts).

---

## II. INSTALLATION

**Python 3.8+** and the following packages are required for correct operation:

pip install numpy pandas tifffile pillow matplotlib scipy scikit-image nd2 czifile oirfile oiffile readlif bioio bioio-openreadout

Launching the program:

python Astro-Analyzer-v1.0.0.py

---

## III. INTERFACE

The main window is divided into two parts:

- **Left — workspace.** In settings mode, the parameters panel is displayed here. After the analysis starts, the work log appears here. You can freely switch between the settings screen and the log using the "← Settings" button (in the log) and "→ Return to log" button (in settings).
- **Right — dialog panel.** All pop-up panels appear here (channel selection, soma placement, manual filtering, etc.). Instructions for each panel are displayed at its bottom.

### The settings panel includes the following fields:

- **File/folder.** Interactive selection of a file or folder for analysis ("Browse…").
- **Loading mode.** "Single file" (to analyze one stack / image) or "Folder" (to analyze a folder with multiple files).
- **Metadata.** The "Read from file" button scans the selected file and fills in the fields below (XY/Z scale, optical parameters for cells and DAPI) if the metadata was not extracted automatically for some reason. If no metadata is available, the fields keep their defaults.
- **Cell type.** "Astrocytes," "Microglia," or "Neurons." Affects the output folder name and the displayed name in the log.
- **Image type.** "Confocal" or "Epifluorescence." Determines which preprocessing is applied before segmentation (grain index computation for confocal mode, LoG with gamma correction for epifluorescence mode).
- **Deconvolution.** Enables/disables 3D deconvolution in the preprocessing protocol.
- **XY (µm/pixel), Z (µm/slice).** Define the physical scale on which all volumetric and linear measurements depend.
- **Cells: NA, λ, ε.** Optical parameters for deconvolution of the cell channel (numerical aperture, emission wavelength, Wiener regularization).
- **DAPI: NA, λ, ε.** Analogous optical parameters for deconvolution of the DAPI channel (if used).
- **Granules: NA, λ, ε.** Analogous optical parameters for deconvolution of the granule channel (if used).
- **Soma autodetection.** Enables automatic soma detection (this function is currently in Beta testing and works inefficiently).
- **Manual cell filtering.** Enables the launch of the manual correction dialog after automatic segmentation.
- **Metrics…** Opens a window for selecting the metrics that will be included in the final CSV file.

To start the analysis, click **"▶ Start analysis."** To exit — **"Exit."** Language switching (RU/EN) is available in the upper right corner.

---

## IV. WORKFLOW

### 1. Loading and preparing the image
1. Click "Browse…" and select a file or folder. Ideally, the metadata (voxel scale and optical parameters required for deconvolution) will be extracted from the file automatically. If the metadata in the file is missing or incomplete, the program uses default values that can be adjusted manually in the settings panel before launching.
2. In the automatically opened dialog with channel thumbnails, assign the channels to be analyzed (for cells, nuclei, and intracellular granules), then click **Apply**.
3. If the metadata in step 1 was not loaded or was loaded incorrectly — click "Read from file" or enter the metadata manually.
4. Set the remaining parameters — make sure the **Image type** matches your data and select the **Cell type**.
5. Click **▶ Start analysis**.
6. If the image is larger than 512 pixels on any side, the program offers in a separate dialog to downsize it so that the larger side becomes ≤ 512. This significantly speeds up processing. If you agree, the program automatically recalculates the real voxel scale so that all subsequent measurements remain in correct physical units.
7. If the image has a physical size larger than 350 pixels on any side, the program offers in a separate dialog to perform manual cropping. This significantly speeds up processing. If you agree, the next dialog displays an interactive cropping window. On the MIP projection, you can select a rectangular region or specify X/Y coordinates to cut off empty edges or artifacts at the periphery of the image. The **Apply** button crops the stack across all channels. The **Skip** button keeps the original size.
8. Wait for deconvolution to be applied to each channel — simplified for the cell and granule channels (background suppression via a broad blurred PSF, unsharp mask, median filter), and full for the DAPI channel. After that, a separate slice selection dialog is displayed for all channels. If the stack after deconvolution contains clearly junk Z-slices (e.g., empty ones at the beginning or end), they can be disabled to improve analysis accuracy.
9. In the opened dialog with a color MIP over Z (R = lower third, G = middle third, B = upper third), mark all somas. A yellow crosshair follows the cursor on the MIP for precise targeting; XY coordinates are displayed in the upper left corner. Clicking **"Done"** completes the marking. After marking, a confirmation window appears: you can scroll through Z and verify that all points are placed on the correct somas and that all somas have been marked. Clicking **"Start segmentation"** continues the process. **Controls:**
  - **LMB** — place a soma point at this XY (Z is selected automatically based on the maximum signal in a window around the point).
  - **Ctrl+LMB** — delete the nearest placed point.
  - **Shift+LMB** — change the type of the nearest soma (regular / "donut").
  - **RMB** — open a magnifier (enlarged Z-slice locally around the point). In the magnifier, the Z slider scrolls through slices, LMB places a soma precisely at the selected XYZ, OK closes the magnifier.

### 2. Segmentation
10. At the next stage, a soma is automatically grown in 3D from each point, bounded by an ellipsoid around the center. Soma boundaries are determined by the intensity level relative to the starting point. In the vicinity of each soma, a nucleus is searched for in the DAPI channel using adaptive thresholding (triangle threshold) and persistence analysis to separate clumped nuclei. If DAPI is not used, nuclei are not segmented. Next, from each soma, a weight map is computed using Dijkstra's algorithm, where the transition cost is inversely proportional to the brightness, structureness, and filamentness of the local structure (with a penalty for intensity dips). At the final stage, granules are segmented within each cell's mask using the Otsu method (if the corresponding option is enabled by the user).
11. After each segmentation stage, a quality check window opens for each structure (somas, nuclei, cells, granules). If the quality is unsatisfactory, you can indicate the problem:
  - **Merged cells** — return to soma correction (add/remove somas).
  - **Artifacts** — return to Z-slice removal.
  - **Truncated processes** — enable "soft" growth mode.
  - **Extra objects** — manually remove segments by clicking.
- If the quality is satisfactory — simply click "Yes" and the program continues.
12. For each cell, 3D volumes (soma, nucleus, whole cell), MIP projection areas, perimeters, intensities, elongation, solidity, and nuclear-cytoplasmic ratio are automatically computed. The number of processes and their characteristics are separately assessed using the ellipsoid scanning method. If granules are selected, their sizes, counts, shape, and spatial distribution are computed.
13. The result is saved to {base}_measurements.csv and appended to the cumulative all_measurements.csv. **Visualizations** are also generated:
- **3D stacks:** raw signal, DAPI, soma/nucleus/cell masks, RGB cell mask, overlay stack (colored masks over grayscale signal).
- **Z-slice mosaics in PNG:** segmentation of somas, nuclei, cells.
- **Final MIP visualizations:** colored cell contours on gray background and on black background.
- **Deconvolved stacks with MIP images.**

---

## V. OUTPUT FILES

For the file sample.nd2 located in the directory {dir} with cell type "astrocytes," the following folder is created:

{dir}/sample_astrocytes_analysis/<br>
├── crop_info.txt — preprocessing log: allow_downsample, crop_box<br>
├── sample_measurements.csv — main metrics table (row = cell)<br>
├── 3d_stacks/<br>
│   ├── sample_astro_3d.tiff — cell signal<br>
│   ├── sample_dapi_3d.tiff — DAPI signal (if present)<br>
│   ├── sample_soma_ids.tiff — soma IDs<br>
│   ├── sample_nuclei_ids.tiff — nucleus IDs (if present)<br>
│   ├── sample_cells_ids.tiff — cell IDs<br>
│   ├── sample_cells_rgb.tiff — cell coloring on black<br>
│   └── sample_cells_overlay_rgb.tiff — cell coloring over the original image<br>
├── deconvolved/<br>
│   ├── sample_cells_deconvolved.tiff<br>
│   ├── sample_cells_deconvolved_mip.png<br>
│   ├── sample_dapi_deconvolved.tiff (only if a nucleus channel is present)<br>
│   ├── sample_dapi_deconvolved_mip.png (only if a nucleus channel is present)<br>
│   ├── sample_granules_deconvolved.tiff (only if a granule channel is present)<br>
│   └── sample_granules_deconvolved_mip.png (only if a granule channel is present)<br>
visualizations/<br>
    ├── sample_soma_slices.png — 2×5 Z mosaic, soma overlay<br>
    ├── sample_nuclei_slices.png — 2×5 Z mosaic, nucleus overlay<br>
    ├── sample_cells_slices.png — 2×5 Z mosaic, cell overlay<br>
    ├── sample_mask_astrocytes.tiff — cell labels (name hardcoded)<br>
    ├── sample_mask_somas.tiff — soma labels<br>
    ├── sample_mask_nuclei.tiff — nucleus labels<br>
    ├── sample_mask_astrocytes_rgb.tiff — cell coloring (duplicate of cells_rgb)<br>
    ├── sample_mip_colored_bg.png — grayscale MIP + colored masks<br>
    └── sample_mip_colored_black.png — black background + colored masks<br>


A cumulative folder is also created:

{dir}/all_results_astrocytes_3d/<br>
└── all_measurements.csv — morphometry information for cells from all files in {dir}<br>
└── summary_cell_counts.csv — cell count information from all files in {dir}<br>

---

## VI. COLUMNS IN .CSV FILES (sample_measurements.csv and all_measurements.csv)

- cell_id — cell number within this file<br>
- filename — name of the source file<br>
- centroid_x, centroid_y, centroid_z — coordinates of the cell's center of mass (in pixels/voxels)<br>
- cell_volume_um3 — whole-cell volume (soma + processes), µm³<br>
- soma_volume_um3 — soma volume, µm³<br>
- nucleus_volume_um3 — nucleus volume, µm³<br>
- cytoplasm_volume_um3 — cytoplasm volume, µm³<br>
- cell_surface_area_um2 — cell surface area, µm²<br>
- cell_sphericity — cell sphericity<br>
- cell_solidity — cell solidity<br>
- soma_solidity — soma solidity<br>
- nucleus_solidity — nucleus solidity<br>
- cell_elongation_3d — elongation (ratio of principal axes of the inertia tensor) of the cell<br>
- cell_eccentricity — eccentricity<br>
- cell_flatness — flatness<br>
- cell_convexity_3d — convexity (3D)<br>
- cell_extent_3d — bounding box fill fraction (3D)<br>
- cell_mean_intensity — mean signal intensity of the cell<br>
- soma_mean_intensity — mean signal intensity of the soma<br>
- nucleus_mean_intensity — mean signal intensity of the nucleus<br>
- cytoplasm_mean_intensity — mean signal intensity of the cytoplasm<br>
- nuclear_cytoplasmic_ratio — ratio of nucleus volume to cytoplasm volume (0 if no nucleus)<br>
- cell_nucleus_intensity_ratio — ratio of nucleus intensity to cytoplasm<br>
- cell_nucleus_centroid_distance_um — displacement of the nucleus centroid from the cell centroid, µm<br>
- cell_soma_centroid_distance_um — displacement of the soma centroid from the cell centroid, µm<br>
- process_max_count — maximum number of processes detected at any scanning radius<br>
- process_max_radius — radius at which this count is reached<br>
- cell_MIP_projection_perimeter_um — perimeter of the cell projection, µm<br>
- cell_perimeter_um — cell perimeter, µm (for 2D images)<br>
- soma_perimeter_um — soma perimeter, µm (for 2D images)<br>
- nucleus_perimeter_um — nucleus perimeter, µm (for 2D images)<br>
- granule_mean_intensity — mean granule intensity<br>
- granule_count — number of granules<br>
- granule_mean_size_um3 — mean granule size, µm³<br>
- granule_volume_um3 — total granule volume, µm³<br>
- granule_volume_fraction — fraction of granule volume in the cell volume<br>
- granule_spatial_dispersion_um — spatial dispersion of granules, µm<br>
- granule_nearest_neighbor_mean_um — mean nearest-neighbor distance, µm<br>
- granule_clark_evans_index — Clark–Evans index<br>
- granule_signal_to_background_ratio — granule SBR<br>

---

## VIII. WHAT TO DO WHEN RESULTS ARE POOR

| Problem | Solution |
| :--- | :--- |
| **Runs slowly** | Reduce the input size manually. Enable downsampling. Disable deconvolution if it is not critical. |
| **Dirty segmentation** | Make sure the correct *Image type* is selected. Try background correction in ImageJ. Remove out-of-focus slices at the selection stage. |
| **Cells merge** | Enable manual filtering. Indicate the problem "Somas need correction" and add/remove somas. |
| **Processes are not fully grown** | Enable manual filtering. Indicate the problem "Growth aggressiveness needs to be increased." |
| **Incorrect process attachment** | Enable manual filtering. Use manual removal of extra segments. |
| **Nuclei are not segmented** | Make sure the DAPI channel is correctly assigned in the channel selection dialog. |
| **File does not read** | Check whether nd2 (for .nd2), oirfile (for .oir), oiffile (for .oif), czifile (for .czi) are installed. |

---

## IX. LICENSE AND AUTHORSHIP

**Author** — Vyacheslav Ivanovich Alekseev (c), 2025–2026.

The program is distributed under the **GNU Affero General Public License v3.0 (AGPL-3.0)**. The full license text is available in the LICENSE file at the project root and at: https://www.gnu.org/licenses/agpl-3.0.html.

It is used in laboratory research for quantitative morphometry of nervous tissue cells.

If you use this program in your scientific work, please cite it as software (with a link to the repository and the publication DOI once available).

**WARNING!** Some functions of the program are currently in beta testing. Functionality may change without notice, and errors in segmentation and calculations are possible. Results obtained in the beta version are not recommended for use in publications without prior validation. If you are interested in beta testing — write to (alekseevvi1@my.msu.ru). We welcome feedback: bug reports, improvement suggestions, and examples of data on which the program performs unsatisfactorily.
