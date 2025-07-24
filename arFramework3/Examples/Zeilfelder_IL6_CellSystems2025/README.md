# Diclofenac and acetaminophen dim the acute phase response but amplify expression of the iron regulator hepcidin in liver cancer cells

## Overview

In this work, several ODE-based dynamic models were developed to simulate the interactions between drug-induced 
signaling pathways in liver cancer cells. These models help elucidate the effects of the non-opioid analgesics
diclofenac (*DCF*) and acetaminophen (*APAP*) on the acute phase response and the regulation of the iron hormone hepcidin.

For more details read here: https://doi.org/XXXX

## Repository Structure
- `Data/`: Experimental input data required for the simulations.
- `Models/`: Contains the different mathematical models used in the study.
- `Parameters/`: Stores parameter files for each model, defining the numerical values used during simulations.
- `Results/`: Output files generated from running the simulations, organized by model.

Ready-to-run scripts for running simulations associated with each model are included, named according to the model of interest 
(e.g., `SetupControlIL6.m` for the IL6 model, `SetupPostL1.m` for the Post-L1 model).

## Description of Mathematical Models

- IL-6 Core Model - `control_IL6.def`
This model simulates the core *IL-6* signaling pathway. It serves as the base model to capture the dynamics of 
*IL-6*-mediated responses in liver cancer cells.

- L1 Model - `L1_model.def`
This model incorporates the effects of drugs (*DCF* and *APAP*) on *IL-6* signaling, penalizing, and regularizing 
these effects.

- Post-L1 Model - `PostL1_model.def`
This model is derived from the IL-6 core model after applying L1 regularization to identify drug-modified parameters. 

- BMP Signal Transduction Module - `control_SMAD.def`
A minimalistic model that captures the *SMAD* signaling in response to bone morphogenetic protein 2 (*BMP2*), 
or *DCF*/*APAP* stimulation. 

- Integrative Model - `combined_model.def`
This is the combined model that integrates both *IL-6* signaling and *SMAD* signaling. 
It models the complex responses when combinations of *IL-6*, *BMP2*, *DCF*, and *APAP* are present. 
This model helps in understanding the crosstalk between these pathways.

- PHH Model - `PHH_model.def`
A primary human hepatocyte model (PMM) which includes drug effects identified from
L1 regularization in HepG2 cells and included parameters reflecting differences between primary
human hepatocytes and HepG2 cells. This model has the same structure as the integrative mondel (`combined_model.def`)

## Prerequisites
The following installations are required for successful use of the scripts:
- MATLAB (version 2019 or higher)
- The modeling software [Data2Dynamics](https://github.com/Data2Dynamics/d2d).

## Example
To simulate and generate figures, run the following scripts one after the other:
`   SetupControlIL6;
    loadParameters;
    arPaperFigures(1, 0, 0, 0);
`

To specify the filename and save generated figures use: `arPaperFigures(1, 0, 0, 0, 'FigureName');` 

## Citation
If you use the models, data, or scripts from this repository in your research, please cite the original paper:

XXXXX
 
## Contact
XXXX
