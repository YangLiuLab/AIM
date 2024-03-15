# AIM
Adaptive Intersection Maximization (AIM) is a high-speed drift correction lgorithm for single molecule localization microscopy. 

The details are presented in our paper entitled "Towards drift-free high-throughput nanoscopy through adaptive intersection maximization".

All the codes under \DME_RCC are from https://github.com/qnano/drift-estimation published in Jelmer Cnossen, Tao Ju Cui, Chirlmin Joo, and Carlas Smith, "Drift correction in localization microscopy using entropy minimization," Opt. Express 29, 27961-27974 (2021).

## Hardware requirement: 
AIM requires only a standard computer. 
RCC and DME require a minimal of 32 GB RAM for the big datasets from large field of view system (e.g., 2048 x 2048).

## Software requirement:
The provided codes have been tested on MATLAB version 2020b to 2023a on Windows 10 Operating System.

## Installation:
Users can direacly download the codes and run the demo code on MATLAB. 
Users need to replace the file name when processing users' own datasets.

## Demo datasets:
We provided four experimental datasets (Origami_PAINT, Microtublue_3d, Tissue_colon and CTCF_MCF10A_DRB_6h) and one simulated dataset (simulationSMLM) in MATLAB mat format available at <a href = "https://doi.org/10.5061/dryad.2v6wwpzw3" title = "[Dryad](https://doi.org/10.5061/dryad.2v6wwpzw3)"> Dryad</a>. Please download these dataset and put them in the <b>Data</b> folder.

## Example files
We provide four MATLAB codes as examples to demonstrate how to use AIM. 
<p>
<b>example_ExperimentalData.m: </b>This code performs drift correction with AIM on 2D or 3D localization coordinates of experimental data. Sample experimental data are avaialble at <a href = "https://doi.org/10.5061/dryad.2v6wwpzw3" title = "[Dryad](https://doi.org/10.5061/dryad.2v6wwpzw3)"> Dryad</a>.
</p>
<p>
<b> example_code_2D.m </b>: This code compares the performance of drift correction for AIM, RCC and DME using 2D localization coordinates for experimental dataset of DNA origami (Origami_PAINT.mat) or simulated data (simulationSMLM.mat) available at <a href = "https://doi.org/10.5061/dryad.2v6wwpzw3" title = "[Dryad](https://doi.org/10.5061/dryad.2v6wwpzw3)"> Dryad</a>.
</p>
<p>
<b> example_code_3D.m </b>: This code compares the performance of drift correction for AIM, RCC and DME using 3D localization coordinates of experimental data of simulated data or or experimental data of microtubules (Microtublue_3d.mat) available at <a href = "https://doi.org/10.5061/dryad.2v6wwpzw3" title = "[Dryad](https://doi.org/10.5061/dryad.2v6wwpzw3)"> Dryad</a>.
  </p>
<p>
<b> example_code_FigureS1.m: </b> This code is used to reproduce Supplementary Figure S1, which shows drift tracking precision under a wide range of image sizes from 128×128 pixels to 2048×2048 pixels.
</p>

## Other files
<p>
<b> simulationSMLM.m: </b> This code is used to generate the simulated SMLM dataset from DNA origami structures used in Figure 2 in the main text.
</p>
<p>
<b> Load_ThunderSTORM.m: </b> This code is used to provide a MATLAB function to read the localization dataset (csv files) from the commonly used ThunderSTORM software.
</p>
