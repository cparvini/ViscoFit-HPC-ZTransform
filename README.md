# ViscoFit-HPC-ZTransform
 High-Powered Cluster Computing Analysis of Cell Viscoelasticity using the Z-Transform Inversion Technique

## Description
This repository is designed to be built using the built-in Matlab C-Compiler (mcc) and then executed via Slurm scripts on a high-powered computing cluster. The methodologies built into the ViscoFit and ViscoFitZ classes are derived from the principles and best practices described by Uluutku, Lopez-Guerra and Solares [(found here)](https://www.beilstein-journals.org/bjnano/articles/12/79). 

This approach utlizes the Z-Transform to acquire implied viscoelastic quantities from force-indentation experiments *without* the assumption of an inherent material model. As such, it represents a marked improvement over the state-of-the-art by relaxing the restrictive assumptions required by traditional continuum mechanics material descriptions. In effect, it steps away from parameterizing models with sometimes unknown or ambiguous applicability and towards more direct measurement of viscoelastic properties.

This repository is largely experimental---specific choices have been made to modify the Z-Transform recommendations according the needs of a particular research project involving high-resolution viscoelastic imaging of living cells via [Quantitative Imaging (QI)](https://www.bruker.com/en/products-and-solutions/microscopes/bioafm/bioafm-accessories/qi-advanced-mode.html). Users are cautioned to investigate the fundamental assumptions made in subfunctions and class descriptions to ensure they meet their own requirements before attempting to fork and utilize this code.

## Acknowledgements
Some of the underlying code to import the exceedingly complex binary format of QI files has been drawn from publically available repositories. The author acknowledgements are left in the function descriptions. In every other circumstance, the code was created by [Cameron Parvini](https://scholar.google.com/citations?hl=en&user=NCzlZ4YAAAAJ).
