# Dose-Dependent Retarded Transient Function (RTF) Approach 
## Example Files for Dose-Dependent RTF Modeling using Data2Dynamics

The dose-dependent RTF approach is an alternative to ODE modeling for signaling pathways. It produces intuitively interpretable parameters that describe signal dynamics and enables predictive modeling of time- and dose-dependencies, even if only individual cellular components are quantified. Moreover, it provides an intuitive way to characterize signaling differences between biological conditions or cell types.

This folde contains files that exemplify the usage of the dose-dependent RTF approach in Data2Dynamics (d2d). The provided scripts and data were used to produce the results presented in the paper *Dynamic Modeling of Signaling Pathways When ODEs Are Not Feasible* by Rachel et al. The data consists of time- and dose-dependent measurements from inflammasome activation in bone marrow-derived macrophages (BMDMs) treated with nigericin sodium salt.

**Note**: The example focuses only on the sustained part of the RTF due to the characteristics of the data used. It can be easily expanded to include the transient part of the RTF by making minor modifications to the model definition files.

**Subfolders:**
* **TimeDependent:** Script uses a single-dose RTF to model each dose separately. Only wild-type data is used.
* **DoseDependent:** Script uses the dose-dependent RTF to model all doses together. Only wild-type data is used.
* **ConditionDependent:** Script uses the dose-dependent RTF with condition-dependent expansion to model all doses together across wild-type and knockout data.
* **Functions:** Additional utility functions for the model.


**Cite:**
* [A New Approximation Approach for Transient Differential Equation Models](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2020.00070) - Kreutz C. *Front. Phys.*, 8, 70, 2020.
* [Dynamic Modeling of Signaling Pathways When ODEs Are Not Feasible](https://) - Rachel T., et al. *Bioinformatics*.
