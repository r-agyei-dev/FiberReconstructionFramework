## Project Overview

This project presents a complete computational workflow for automated **2D fiber segmentation and 3D fiber reconstruction** from volumetric tomographic datasets.

The framework progresses systematically from slice-wise feature extraction to full volumetric reconstruction and orthogonal validation, forming a robust and scalable pipeline for microstructural fiber characterization.

The analysis includes:

* Section-wise memory-efficient processing of large grayscale volumes
* Automated 2D segmentation of fiber cross-sections
* Ellipse fitting and geometric feature extraction
* Computation of centroid coordinates and orientation angles
* Stacking of segmented ellipses into crude 3D fiber morphologies
* Detection and separation of conjoined fiber structures
* Re-stitching of segmented fiber volumes
* Mutual orthogonal comparison of reconstructions along x, y, and z directions
* Export of reconstruction results in `.mat`, `.h5`, and `.xdmf` formats

The utility of this reconstruction framework is demonstrated using tomographic datasets (e.g., S1 and S2 at initial and final loading states), illustrating how volumetric microstructural information can be extracted, quantified, and reconstructed into analyzable 3D morphologies.

This repository is organized to reflect professional computational research standards commonly used in materials science laboratories, computational mechanics research groups, and high-performance imaging workflows.

---

## Computational Methodology

The framework follows a structured multi-stage methodology designed for robustness, scalability, and reproducibility.

### 1. Volume Sectioning and Preprocessing

High-resolution grayscale tomographic volumes are partitioned into equal sections to ensure efficient memory usage. Each section is processed independently and later concatenated to reconstruct the full dataset. This approach enables large-scale processing without exceeding system memory limitations.

### 2. 2D Fiber Segmentation

For a selected principal axis (x, y, or z), slice-wise segmentation is performed on the grayscale data. The segmentation workflow consists of:

* Thresholding and region isolation
* Morphological filtering to reduce noise and artifacts
* Connected component labeling
* Ellipse fitting to each uniquely segmented fiber cross-section

From each segmented region, geometric descriptors are extracted, including:

* Centroid coordinates
* Major and minor axis lengths
* In-plane orientation
* Out-of-plane orientation

These descriptors serve as the fundamental inputs for 3D reconstruction.

### 3. 3D Fiber Reconstruction

The 3D reconstruction stage transforms independent 2D ellipses into continuous volumetric fiber morphologies through:

* Stacking consecutive ellipses based on spatial proximity and geometric consistency
* Formation of crude 3D fiber volumes
* Identification of merged or conjoined fibers
* Segmentation of incorrectly merged structures
* Re-stitching of separated fiber segments to recover correct topology

This iterative process ensures geometrically consistent and physically meaningful fiber morphologies.

### 4. Mutual Orthogonal Comparison

To enhance reconstruction reliability, segmentation and reconstruction are performed independently along all three principal axes. The reconstructed volumes are compared orthogonally to:

* Validate fiber continuity
* Reduce axis-dependent reconstruction bias
* Improve morphological accuracy

This cross-validation step significantly increases robustness in complex fiber networks.

### 5. Data Export and Output Structure

Final outputs are stored in structured formats for analysis and visualization:

* `.mat` files containing centroid coordinates and pixel indices
* `.h5` files containing full 3D reconstructed volumes
* `.xdmf` files for visualization in ParaView or similar platforms

Condensed datasets are also generated containing fiber-level descriptors such as:

* Start and end coordinates
* Fiber length
* In-plane and out-of-plane angles
* Unique fiber label identifiers

---

## Implementation Design

The framework is implemented entirely in MATLAB using a modular architecture.

Each major computational stage is controlled by a primary driver script:

* `MAIN_2D_SEGM_MFILE.m`
* `MAIN_3D_RECON_MFILE.m`
* `Compare_MAIN_mfile.m`

Supporting sub-functions handle image processing operations, parallel execution, feature extraction, volume stitching, and data export.

The modular design ensures:

* Scalability to large datasets
* Compatibility with parallel computing environments
* Reproducibility of segmentation and reconstruction steps
* Extensibility for future algorithmic development

Parallelization is supported via MATLAB’s Parallel Computing Toolbox, allowing users to specify the number of workers for improved performance.

Memory efficiency is maintained through section-wise processing and structured intermediate data handling, enabling reconstruction of high-resolution tomographic datasets without requiring high-performance computing clusters.

---

Together, these components form a complete supervised iterative fiber reconstruction framework suitable for advanced microstructural analysis and mechanical behavior investigations.
