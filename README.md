# Surrogate Modeling for Irreversible Electroporation
This code tries to build a surrogate model for the quantities of interest related to irreversible electroporation. In particular, we try to model the relative tumor ablation in the case of irreversible electroporation for colorectal metastasis in the liver.

# Getting Started
Set the variable `phases` in [`main.m`] and run the script.

## Prerequisites
- UQLab - The Framework for Uncertainty Quantification (www.uqlab.com).
- COMSOL Multiphysics v5.6 or higher.
- MATLAB version R2017a or higher.

# Extension to other applications
To reuse the code for different applications, change the functions `solve_COMSOL_EPlist`, `solve_COMSOL_EP`, and `createInputs` in [`UtilFuncs.m`].
