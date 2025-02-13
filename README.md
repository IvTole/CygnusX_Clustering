# CygnusX_Clustering

We use data from the Nobeyama 45m radiotelescope towards CygnusX region. The objective is to apply a clustering technique (python package astrodendro) to detect c18o clumps in the region.

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

## Conda environment setup

Inside directory **environment/** there is a file named **environment.yml**. This file is used to set up a dedicated Conda environment with all the necessary dependencies for running the code in this repository.

To create the environment, first ensure you have **Anaconda** or **Miniconda** installed on your system. You can download it from [Anaconda's official website](https://www.anaconda.com/download). Then, open a terminal and run the following command:

`conda env create -f environment.yml`

```bash
conda env create -f environment/environment.yml
```




