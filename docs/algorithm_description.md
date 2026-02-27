# Algorithm Description

This document provides a detailed description of the algorithms used in the fiber reconstruction framework.

## 1. 2D Segmentation Algorithm

**Input:** Grayscale 2D slice  
**Output:** Ellipse parameters for each fiber

Steps:

1. Threshold image to isolate potential fiber regions
2. Apply morphological operations to clean artifacts
3. Label connected components
4. Fit ellipses to each component
5. Extract geometric descriptors:
   - Centroid (x, y)
   - Major/minor axis lengths
   - Orientation (in-plane, out-of-plane)

---

## 2. 3D Reconstruction Algorithm

**Input:** Stacked 2D ellipses  
**Output:** 3D fiber volumes

Pseudocode:
for each axis in [x, y, z]:
for each section:
segment 2D fibers
extract ellipse features
end
stack 2D ellipses into 3D fibers
detect conjoined fibers
separate merged fibers
re-stitch segments
end
compare reconstructions across axes
export final volumes

---

## 3. Mutual Orthogonal Comparison

- Compare reconstructed volumes from all three axes
- Identify fibers with inconsistent topology
- Adjust segmentation or stitching if necessary
