ASTRO-ANALYZER v0.1.0-beta.1

This program is designed for automatic segmentation and morphometry of astrocytes and microglial cells in 3D stacks and 2D images obtained by confocal and epifluorescence microscopy (tissue sections or cell cultures). The output of this program is a table of morphological parameters for each cell: soma volume, nucleus volume, whole-cell volume with processes, projection areas, mean intensities, elongation, nuclear-to-cytoplasmic ratio, process branching complexity, and other features.

I. PIPELINE

For each input file (.tif / .nd2 / .czi / .lif / .ome.tif stack), the program:

1. Reads the image and metadata (XY and Z scale, channels, optical parameters).
2. If requested by the user, performs 3D deconvolution to suppress out-of-focus haze.
3. Asks the user to place soma centers on the maximum intensity projection.
4. Automatically outlines the soma of each marked cell in 3D.
5. Automatically finds the nucleus of each cell in the DAPI channel.
6. Grows the whole cell (soma + processes) in the cell channel.
7. Allows manual correction or deletion of problematic cells.
8. Computes morphological parameters for each cell and saves them to CSV.
9. Generates visualizations of each cell's masks.
10. Aggregates results from all files in a folder into a cumulative CSV.

II. INSTALLATION

For correct operation, the program requires Python 3.8+ and the following packages:
pip install numpy pandas tifffile pillow matplotlib scipy scikit-image nd2 czifile readlif

To run the program, launch Astrocyte_analyzer_refactoring_GUI.py via:
python Astrocyte_analyzer_refactoring_GUI.py

III. INTERFACE

The main window of the program is divided into two parts:

- Left — workspace. In settings mode, the parameter panel is displayed here. After analysis starts, the work log appears here. You can freely switch between the settings screen and the log using the buttons "← Settings" (in the log) and "→ Back to log" (in the settings).
- Right — dialog panel. All pop-up panels appear here (channel selection, soma marking, manual filtering, etc.). Instructions for each panel are shown at its bottom.

The settings panel includes the following fields:

- File/folder. Offers interactive selection of a file or folder for analysis ("Browse…").
- Load mode. Offers selection of the operation mode — "Single file" (to analyze one stack / image) or "Folder" (to analyze a folder with multiple files).
- Metadata. The "Read from file" button scans the selected file and fills in the fields below (XY/Z scale, cell and DAPI optical parameters). If metadata is absent, the fields keep their defaults.
- Cell type. Offers selection of the analyzed cell type — "Astrocytes" or "Microglia". The choice affects the output folder name and the name displayed in the log.
- Image type. Offers selection of the image acquisition method — "Confocal" or "Epifluorescence". The choice determines which preprocessing is applied before segmentation. When "Confocal" is selected, preprocessing based on multi-scale grain index with slice-wise brightness normalization is used (suitable for confocal images). When "Epifluorescence" is selected, LoG preprocessing with gamma correction is used, which removes the broad diffuse halo and autofluorescence (suitable for epifluorescence images).
- Deconvolution. Enables or disables 3D deconvolution in the preprocessing pipeline. It is recommended to keep it enabled for confocal stacks with noticeable out-of-focus haze.
- XY (µm/pixel), Z (µm/slice). Define the physical scale on which all volumetric and linear measurements depend.
- Cells: NA, λ, ε. Define the optical parameters for cell-channel deconvolution (numerical aperture, emission wavelength, Wiener regularization).
- DAPI: NA, λ, ε. Define the optical parameters for DAPI-channel deconvolution.
- Manual cell filtering. Enables launching the manual correction dialog after automatic segmentation. If disabled, the program proceeds directly to measurements.

To start the analysis, click "▶ Run analysis" at the bottom of the screen. To exit the interface, click "Exit" at the bottom of the screen.

IV. WORKFLOW

1. Click "Browse…" and select a file or folder.
2. If needed, click "Read from file" to load metadata.
3. Make sure the Image type matches your data.
4. Select the Cell type.
5. Leave the other parameters at defaults or configure them manually.
6. Click "▶ Run analysis".
7. Complete the manual dialogs (see below).
8. When finished, view the results in the folder <file_name>_<cell_type>_analysis/.

V. FILE PROCESSING STAGES AFTER CLICKING "▶ Run analysis"

1. Image loading
- The file is opened, the axis order and number of channels are determined.
- Voxel scale and optical parameters required for deconvolution are extracted.
- If the image exceeds 512 pixels on either side, it is automatically downsampled so that the larger side becomes ≤ 512 (this speeds up processing). Important: after downsampling, the real voxel scale is automatically recalculated so that measurements remain in correct physical units.

2. Channel determination
- A channel selection dialog opens showing thumbnails (MIPs) of all channels found in the file. Default values in the dialog are taken from the image metadata: if the program confidently identifies the DAPI channel (by name, excitation, or emission wavelength), it will be preset as DAPI, and the remaining channel as cells. If metadata is absent or DAPI could not be determined, channel 0 — cells, channel 1 — DAPI is suggested by default. Channels can be reassigned manually in the dialog:
  LMB on the channel image — assign it as the cell channel.
  RMB on the channel image — assign it as the DAPI channel.
  LMB on the bottom button under the channel — include/exclude the channel from analysis.
  "No DAPI channel" checkbox — disable DAPI usage (nuclei will not be segmented).
- The "Apply" button fixes the selection. The "Cancel" button uses the default values without changes.

3. Deconvolution (if enabled)
- Simplified deconvolution of the cell channel is applied (background suppression via a wide blurred PSF, unsharp mask, median filter).
- DAPI-channel deconvolution is applied (PSF is calculated from the Airy radius based on NA and emission wavelength).
- For both channels, a slice-selection dialog is offered. If the deconvolved stack contains obviously garbage Z-slices (e.g., empty at the start or end), they can be disabled to increase analysis accuracy. The deconvolved stack is saved to the deconvolved folder.

4. Soma marking
- A dialog opens with a color-coded Z-MIP (R = lower third, G = middle third, B = upper third). In this dialog, the user is asked to manually place the centers of cells (somas). Small variations between operators in soma placement do not affect segmentation quality. Controls:
  LMB — place a soma point at this XY (Z is chosen automatically by the maximum signal in a window around the point).
  Ctrl+LMB — delete the nearest placed point.
  RMB — open a magnifier (enlarged Z-slice locally around the point). In the magnifier, the Z slider browses slices, LMB places the soma at the exact XYZ, OK closes the magnifier.
  A yellow crosshair follows the cursor on the MIP for precise targeting; XY coordinates are shown in the upper left corner.
- Clicking "Done" finishes marking.
- After marking, a confirmation window appears: you can browse Z to check that all points are on the correct somas and that all somas have been marked.
- Clicking "Start segmentation" continues the process.

5. Segmentation
- From each point, a soma is grown in 3D via BFS, bounded by an ellipsoid around the center. Soma boundaries are determined by the intensity level relative to the seed point.
- In the neighborhood of each soma, a nucleus is searched in the DAPI channel. Otsu thresholding on the local region + watershed for splitting merged nuclei is used. If DAPI is not used, nuclei are not segmented.
- Preprocessing depends on the image type (see above).
- A "structure" map (grain index / LoG) is built on the preprocessed signal.
- Pixels with high structure fall into the growth mask.
- From each soma, a weight map is computed using Dijkstra's algorithm on a weight map where the transition cost is inversely proportional to brightness, structure, and the filamentous nature of the local structure containing the voxel (with a penalty for brightness dips).
- Voxels are distributed between somas by the principle of "connectivity with a specific soma" (each connected component is assigned to the soma to which it is connected by a structural path). Isolated components without somas are discarded.

6. Manual filtering (optional)

- A MIP map of all found cells opens with a dialog allowing re-segmentation of problematic cells. Controls:
  LMB on a cell — mark as UNDERSIZED (processes need to be grown).
  RMB on a cell — mark as OVERSIZED (cell needs to be narrowed).
- The "Continue" button starts re-segmentation of the selected cells. The "Skip" button skips filtering.
- During re-segmentation, the same algorithm as in the main segmentation is applied to each problematic cell, but with different parameters.
- After re-segmentation, a review dialog opens to decide on the re-segmented cells. Controls:
  LMB on a cell — accept the new segmentation (green).
  RMB on a cell — delete the cell (red).
  No click (then pressing "Continue") — revert to the old segmentation.

7. Measurements
- For each cell, 3D volumes (soma, nucleus, whole cell), MIP projection areas, perimeters, intensities, elongation, solidity, nuclear-to-cytoplasmic ratio are computed. Separately, the number of processes and their characteristics are estimated using the ellipsoid scanning method.
- The result is saved to {base}_measurements.csv and appended to the cumulative all_measurements.csv.

8. Visualizations
- 3D stacks: raw signal, DAPI, soma/nuclei/cell masks, RGB cell mask, overlay stack (colored masks over the gray signal).
- Z-slice mosaics as PNG: soma, nuclei, cell segmentation.
- Final MIP visualizations: colored cell contours on gray background and on black background.
- Deconvolved stacks with MIP images.

VI. OUTPUT FILES

For a file sample.nd2 located in directory {dir} with cell type "astrocytes", the following folder is created:
{dir}/sample_astrocytes_analysis/
├── 3d_stacks/
│ ├── sample_astro_3d.tiff — raw cell signal (16-bit)
│ ├── sample_dapi_3d.tiff — DAPI (16-bit)
│ ├── sample_nuclei_ids.tiff — nucleus labels
│ ├── sample_astrocytes_ids.tiff — cell labels
│ ├── sample_astrocytes_rgb.tiff — RGB cell mask
│ └── sample_cells_overlay_rgb.tiff — colored masks over the raw grayscale image
├── deconvolved/
│ ├── sample_cells_deconvolved.tiff
│ ├── sample_cells_deconvolved_mip.png
│ ├── sample_dapi_deconvolved.tiff
│ └── sample_dapi_deconvolved_mip.png
├── visualizations/
│ ├── sample_soma_slices.png — soma segmentation mosaic
│ ├── sample_nuclei_slices.png — nuclei segmentation mosaic
│ ├── sample_cells_slices.png — cell segmentation mosaic
│ ├── sample_mask_astrocytes.tiff
│ ├── sample_mask_somas.tiff
│ ├── sample_mask_nuclei.tiff
│ ├── sample_mask_astrocytes_rgb.tiff
│ ├── sample_mip_colored_bg.png
│ └── sample_mip_colored_black.png
├── manual_filtering/
│ └── sample_manual_filtering.json — manual selection decisions
└── sample_measurements.csv — parameter table

A cumulative folder is also created:
{dir}/all_results_astrocytes_3d/
└── all_measurements.csv — information on all cells from all files in {dir}

VII. COLUMNS IN .CSV FILES (sample_measurements.csv and all_measurements.csv)

cell_id — cell number in this file
filename — source file name
centroid_x, centroid_y, centroid_z — cell centroid coordinates (in pixels/voxels)
cell_volume_um3 — whole cell volume (soma + processes), µm3
soma_volume_um3 — soma volume, µm3
nucleus_volume_um3 — nucleus volume, µm3
cell_area_um2 — whole cell area, µm2 (for 2D images)
soma_area_um2 — soma area, µm2 (for 2D images)
nucleus_area_um2 — nucleus area, µm2 (for 2D images)
cell_solidity — cell solidity (area to convex hull area ratio)
soma_solidity — soma solidity
nucleus_solidity — nucleus solidity
elongation — elongation (ratio of principal axes of the inertia tensor) of the cell
cell_mean_intensity — cell signal intensity
soma_mean_intensity — soma signal intensity
nucleus_mean_intensity — nucleus signal intensity
nuclear_cytoplasmic_ratio — nucleus-to-cytoplasm volume ratio (0 if no nucleus)
cell_MIP_projection_perimeter_um — cell projection perimeter, µm
soma_MIP_projection_perimeter_um — soma projection perimeter, µm
nucleus_MIP_projection_perimeter_um — nucleus projection perimeter, µm
cell_perimeter_um — cell perimeter, µm (for 2D images)
soma_perimeter_um — soma perimeter, µm (for 2D images)
nucleus_perimeter_um — nucleus perimeter, µm (for 2D images)
process_max_count — maximum number of processes detected at any scanning radius
process_max_radius — radius at which this number was reached
process_avg_mid — average number of processes in the middle radius zone
process_integral — integral of the "number of processes vs radius" profile
process_decay_ratio — ratio of the number of processes at the last radius to the maximum
process_branching_complexity — standard deviation of the profile derivative (branching measure)
process_half_max_radius — radius at which the number of processes drops by half
process_avg_size — average process size (voxels)
process_total_voxels — total number of voxels in all processes
process_to_cell_ratio — ratio of process volume to cell volume

VIII. WHAT TO DO WITH POOR RESULTS

Runs slowly:
- Reduce the input size manually.
- Disable deconvolution if it is not critical.

Dirty segmentation:
- Make sure the correct Image type is selected.
- Try background correction with other software (ImageJ).
- Try removing out-of-focus slices at the slice selection stage (especially important for epifluorescence images).

Cells merge together:
- Use manual filtering (enable it in settings). Mark merged cells as OVERSIZED (RMB).

Processes are not grown:
- Use manual filtering (enable it in settings). Mark such cells as UNDERSIZED (LMB).

Incorrect attachment of processes to another cell:
- Use manual filtering (enable it in settings). Mark the conflicting cells as OVERSIZED and UNDERSIZED respectively.

Nuclei are not segmented:
- Make sure the DAPI channel is correctly identified in the channel selection dialog.

File cannot be read:
- Check whether nd2 (for .nd2) or czifile (for .czi) is installed.

IX. LICENSE AND AUTHORSHIP

Author — Alekseev Vyacheslav Ivanovich (c), 2026.

The program is distributed under the GNU Affero General Public License v3.0 (AGPL-3.0). The full license text is available in the LICENSE file at the project root and at: https://www.gnu.org/licenses/agpl-3.0.html.

Used in laboratory research for quantitative morphometry of astrocytes.

If you use this program in your scientific work, please cite it as software (with a link to the repository and the DOI of the publication when it becomes available).

ATTENTION! The program is currently in beta testing. Functionality may change without notice; segmentation and calculation errors are possible. Results obtained in the beta version are not recommended for use in publications without prior validation. If you are interested in beta testing, please write to (alekseevvi1@my.msu.ru). We welcome feedback: bug reports, improvement suggestions, and examples of data on which the program performs unsatisfactorily.