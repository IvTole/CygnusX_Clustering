# CygnusX_Clustering

We use data from the Nobeyama 45m radiotelescope towards CygnusX region. The objective is to apply a clustering technique (python package astrodendro) to detect c18o in the region.

## Repository Structure
```
.
├── data/         # Original data cubes in .fits format (for 12co, 13co and c18o (J=1-0)).
├── fitsfiles/    # Output fits files from the script. At the moment, this includes mom8 and smoothed cubes files.
├── modules/      # Python scripts for running the program and modules. Main script : run_clustering.py
├── plots/        # Output plots from the script.
├── catalog/      # Output catalog (csv) from the script. This includes catalogs from astrodendro, fit parameters from Gaussian fit, and final catalog of clumps with physical parameters.
├── mask/         # Output mask (numpy array) from the script, where clusters are identified.
├── environment/  # Conda environment for python.
└── README.md     # Repository documentation

```