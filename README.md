# ViscoFit-HPC-ZTransform
 High-Powered Cluster Computing Analysis of Cell Viscoelasticity using the Z-Transform Inversion Technique

## Description
This repository is designed to be either built using the built-in Matlab C-Compiler (mcc) and then executed via Slurm scripts on a linux-based high-powered computing cluster, or executed directly by that cluster using the "matlab" environment argument (no compiling necessary; this is how we executed our analysis). The methodologies built into the ViscoFit and ViscoFitZ classes are derived from the principles and best practices described by Uluutku, Lopez-Guerra and Solares [(found here)](https://www.beilstein-journals.org/bjnano/articles/12/79).

This approach utlizes the Z-Transform (a type of integral transform, similar to the modified Fourier) to acquire implied viscoelastic quantities from force-indentation experiments *without* the assumption of an inherent material model. As such, it represents a marked improvement over the state-of-the-art by relaxing the restrictive assumptions required by traditional continuum mechanics material descriptions. Essentially, it steps away from parameterizing models (with sometimes unknown or ambiguous applicability) and towards more direct measurement/observation of viscoelastic properties. The function that executes this analysis is analyze_map_zTransform(), but an alternative version (fit_map_zTransform()) will extract the same information via Z-Transform and then attempt to parameterize one of the common generalized rheological models (Generalized Maxwell, Generalized Kelvin-Voigt, or PLR) using that data. Note that the fitting functionality has *not* been tested thoroughly, and is mostly based upon legacy code from a previous publication [(found here)](https://www.nature.com/articles/s42003-021-02959-5).

This repository is experimental---specific choices have been made to modify the Z-Transform recommendations according the needs of a particular research project involving high-resolution viscoelastic imaging of living cells via [Quantitative Imaging (QI)](https://www.bruker.com/en/products-and-solutions/microscopes/bioafm/bioafm-accessories/qi-advanced-mode.html). Users are cautioned to investigate the fundamental assumptions made in subfunctions and class descriptions to ensure they meet their own requirements before attempting to fork and utilize this code.

To modify this repository for your own use-case, users are directed to the LoadAFMData() function which performs the loading, correction, and initial configuration before performing the integral transform. Here, users will have to ensure that a case exists to handle the data format they use. Either the file extension will be used, or a trigger within the filename to identify which data loading procedure to use. Alternatively, users can create a new copy of analyze_map_zTransform() and modify it to load files directly, though the data streams will not be corrected by the code in the same manner as our datasets.

Should users require clarification for any code contained within this repository, they are encouraged to reach out via github or email.

## Requirements
- Matlab (>= 2021a), and...
  - Parallel Processing Toolbox
  - Signal Processing Toolbox
  - Statistics & Machine Learning Toolbox
  - Image Processing Toolbox
  - Matlab Compiler
- Your Preferred IDE (for modifying SLRUM sripts; I use VS Code)

## Acknowledgements
Some of the underlying code to import the exceedingly complex binary format of QI files has been drawn from publically available repositories. The author acknowledgements are left in the function descriptions. In every other circumstance, the code was created by [Cameron Parvini](https://scholar.google.com/citations?hl=en&user=NCzlZ4YAAAAJ).
