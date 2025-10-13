# Parallel Frequency-Multiplexed Aberration Measurement (PFAM) for Widefield Fluorescence Microscopy

This is the public repository for **PFAM** — *Parallel frequency-multiplexed aberration measurement for widefield fluorescence microscopy*.  
The original preprint is available on **bioRxiv**: **https://www.biorxiv.org/content/10.1101/2025.10.11.681535v1**.
**Data access:** The associated data can be provided via Google Drive upon request ([here](https://drive.google.com/file/d/1HjPRzHwD6yLF7Xw8yzognaLx_FoOjXlm/view?usp=sharing)).

---

## Abstract

Widefield fluorescence microscopy is widely used for imaging at subcellular resolution, but its performance in complex samples is degraded by optical aberrations. Because aberrations can vary spatially across the field of view (FOV), accurate aberration measurement and correction at multiple FOV locations are essential for achieving high-quality imaging over large areas. Here, we introduce **parallel frequency-multiplexed aberration measurement (PFAM)** to perform massively parallel aberration measurements across an extended FOV. We validated PFAM using fluorescent beads and demonstrated simultaneous measurement and effective correction of spatially varying aberrations at **125** FOV locations. To address challenges of wavefront sensing in complex samples, we further developed **PFAM-SIFT** by integrating structured illumination, achieving robust aberration measurement in both brain slices and the mouse cortex *in vivo*. Together, PFAM and PFAM-SIFT provide accurate and scalable wavefront-sensing solutions for widefield imaging, enabling simultaneous measurement and correction of spatially varying aberrations in complex biological samples.

---

## MATLAB Dependencies

- **MATLAB R2023b** (or later)
- **Curve Fitting Toolbox**
- **Image Processing Toolbox**
- **Optimization Toolbox**
- **Statistics and Machine Learning Toolbox**

> Fast TIFF reading uses the MATLAB File Exchange function cited below.

---

## Data

- Place each dataset in the corresponding **code** subfolder as a `data/` directory (e.g., `Code1_PFAM_single_ROI/data/`, `Code2_PFAM_multiple_ROIs/data/`) **before running**.  
- Geometric configuration and segment indices for the deformable mirror (DM) **IrisAO PTT489 segmented DM** and the **Zernike polynomial** definitions are provided in the `DM_data/` directory.

---

## Code Overview

Each `code*` folder contains self-contained scripts and (upon request) example data for a specific PFAM/PFAM-SIFT component.

- **Code1 – PFAM (single ROI)**  
  Generates tip/tilt interference maps from frequency-multiplexed measurements for a **single** ROI.

- **Code2 – PFAM (multiple ROIs)**  
  Extends Code1 to **multiple** ROIs, enabling parallel measurement across the FOV.

- **Code3 – PFAM-SIFT (single ROI)**  
  Implements PFAM with **structured illumination and Fourier Transform** (SIFT) for a single ROI.  
  Includes a **simulation** of structured-illumination patterns and sidebands.

- **Code4 – PFAM-SIFT (multiple ROIs)**  
  Multi-ROI PFAM-SIFT for robust aberration measurement in complex samples (brain slices and *in vivo*).

- **Code5 – Phase Retrieval**  
  **Gerchberg–Saxton** phase retrieval from 3D image stacks of sub-diffraction fluorescence beads.

- **Code6 – Wavefront Reconstruction**  
  Computes **piston** values from measured **tip/tilt** per DM segment using the **Zonal Matrix Iterative Method**.  
  Includes both **3-group** and **2-group** DM configurations.

**Processing pipeline (Codes 1–4):**
1. Load camera frames acquired during tip/tilt scanning.  
2. Extract fluorescence signal at designated locations per frame.  
3. Generate **tip/tilt interference maps** for each DM segment.  
4. Extract the **peak** of each map (per segment) for wavefront estimation.

---

## External References & Utilities

- **Wavefront reconstruction method**  
  Panagopoulou, S. I., & Neal, D. P. *Zonal matrix iterative method for wavefront reconstruction from gradient measurements*. **Journal of Refractive Surgery**, **21**(6), S563–S569 (2005).

- **Fast TIFF reader (MATLAB File Exchange)**  
  Yoon-Oh Tak (2025). *Multipage TIFF stack*.  
  https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack  
  (Retrieved October 13, 2025.)

---

## How to Run

1. Ensure MATLAB dependencies are installed.  
2. Place the appropriate dataset in each target code folder’s `data/` directory (request Google Drive access if needed).  
3. Open MATLAB in that folder and run the main script(s).

---

## Citation

If you use PFAM/PFAM-SIFT or code from this repository, please cite:

**Paper (bioRxiv)**  
Kim, H., Kang, I., Natan, R., & Ji, N. (2025). *Parallel frequency-multiplexed aberration measurement for widefield fluorescence microscopy*. **bioRxiv**. https://doi.org/10.1101/2025.10.11.681535

**BibTeX**
```bibtex
@article{kim2025pfam_bioRxiv,
  author  = {Kim, Hyeonggeon and Kang, Iksung and Natan, Ryan and Ji, Na},
  title   = {Parallel frequency-multiplexed aberration measurement for widefield fluorescence microscopy},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.10.11.681535}
}
