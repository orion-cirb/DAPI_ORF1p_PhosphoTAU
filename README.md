#  Hoechst_ORF1p_pTAU

* **Developed for:** Rania
* **Team:** Fuchs
* **Date:** November 2023
* **Software:** Fiji

### Images description

3D images of taken on a spinning-disk with a x40 objective

3 channels:
  1. *CSU_405:* Hoechst nuclei
  2. *CSU_561:* ORF1p
  2. *CSU_642:* pTAU cells


### Plugin description

* Detect Hoechst nuclei and pTAU cells with Cellpose
* Compute their colocalization and determine whether each nucleus is pTAU+ or pTAU-
* Compute background noise of ORF1p channel
* For each nucleus, give volume and background-corrected ORF1p intensity

### Dependencies

* **3DImageSuite** Fiji plugin
* **Cellpose** conda environment + *cyto* and *cyto2_pTAU* (homemade) models

### Version history

Version 1 released on November 28, 2023.

