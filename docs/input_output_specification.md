# Input and Output Specification

## Input Data

- Grayscale tomographic volume
- Supported formats: `.mat` or `.h5`
- Required variables:
  - Volume data array
  - Voxel size (if applicable)
- Axis convention: 1 = x-axis, 2 = y-axis, 3 = z-axis

---

## Output Data

### 1. Centroid Data

- File: `S1_3D_Centroid_at_Initial.mat`
- Content: [x, y, z] coordinates of ellipse centroids

### 2. Condensed Fiber Data

- File: `S1_3D_Condensed_Data_at_Initial.mat`
- Columns:
  1. Start x-coordinate
  2. Start y-coordinate
  3. Start z-coordinate
  4. End x-coordinate
  5. End y-coordinate
  6. End z-coordinate
  7. Out-of-plane angle
  8. In-plane angle
  9. Fiber length
  10. Fiber label ID

### 3. Pixel Indices

- File: `S1_Linear_IDX_at_Initial.mat`
- Content: linear indices of all fiber pixels

### 4. Full 3D Reconstruction

- File: `S1_REC_at_Initial.h5`
- Content: volumetric fiber reconstruction

### 5. Visualization

- File: `S1_REC_at_Initial.xdmf`
- For visualization in ParaView or compatible software
