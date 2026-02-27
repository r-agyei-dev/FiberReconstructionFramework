# Methodology

## Overview

This document describes the methodology for automated 2D fiber segmentation and 3D fiber reconstruction from tomographic volumes.

The framework follows a structured multi-stage approach:

1. Section-wise preprocessing of grayscale tomographic volumes
2. Slice-wise 2D fiber segmentation and feature extraction
3. Stacking of 2D ellipses into 3D fiber morphologies
4. Detection and separation of conjoined fibers
5. Orthogonal comparison to validate reconstruction

---

## 1. Volume Preprocessing

High-resolution grayscale volumes are partitioned into equal sections to reduce memory usage. Each section is processed independently, and results are concatenated for full-volume reconstruction.

Key considerations:

- Volume partitioning ensures scalability to large datasets
- Section-wise processing avoids memory overflow
- Sections maintain original spatial coordinates for accurate reconstruction

---

## 2. 2D Fiber Segmentation

For each section along a selected axis (x, y, or z):

- Thresholding isolates potential fiber regions
- Morphological operations remove noise and small artifacts
- Connected component labeling identifies unique regions
- Ellipse fitting is applied to each segmented region

Extracted features include:

- Centroid coordinates
- Major and minor axis lengths
- In-plane and out-of-plane orientations

---

## 3. 3D Fiber Reconstruction

2D ellipse features are stacked to form 3D fibers:

- Consecutive ellipses are matched using spatial proximity and orientation similarity
- Crude 3D fiber volumes are generated
- Merged or conjoined fibers are detected and segmented
- Re-stitching corrects topology of separated fibers

---

## 4. Mutual Orthogonal Comparison

Reconstructions along x, y, and z axes are compared to:

- Validate fiber continuity
- Reduce axis-dependent bias
- Improve morphological accuracy

---

## 5. Data Export

Final outputs:

- `.mat` files: centroid coordinates, pixel indices
- `.h5` files: full volumetric reconstructions
- `.xdmf` files: visualization in ParaView or compatible tools
- Condensed datasets: fiber length, start/end coordinates, orientation, label IDs
