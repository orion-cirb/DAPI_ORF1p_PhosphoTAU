#  Hoechst_ORF1p_pTAU

* **Developed for:** Rania
* **Team:** Fuchs
* **Date:** November 2023
* **Software:** Fiji

### Images description

3D images taken on a spinning-disk microscope with a x60 objective.

3 channels:
  1. *CSU_405:* Hoechst nuclei
  2. *CSU_561:* ORF1p
  2. *CSU_642:* pTAU cells (optional)


### Plugin description

* Detect Hoechst nuclei with Cellpose
* If channel provided, detect pTAU cells with Cellpose, compute their colocalization with nuclei and determine whether each nucleus is pTAU+ or pTAU-
* Compute background noise of ORF1p channel
* For each nucleus, give volume and background-corrected ORF1p intensity

### Dependencies

* **3DImageSuite** Fiji plugin
* **Cellpose** conda environment + *cyto* and *cyto2_pTAU* (homemade) models

### Version history

Version 1 released on November 30, 2023.

