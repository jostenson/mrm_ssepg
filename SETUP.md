# SETUP #

=============================

Ostenson, Smith, Does, Damon. "Slice-selective extended phase graphs in gradient-crushed, transient-state free precession sequences: an application to magnetic resonance fingerprinting". Magnetic Resonance in Medicine

Requirements:
-------------
* MATLAB  - for reconstruction, processing, and figure generation
* C compiler - for generating MEX files from contributing code
* CUDA compiler - for generating ptx file for GPU dictionary processing
* Berkeley Advanced Reconstruction Toolbox (BART) - for coil combination and image reconstruction (https://mrirecon.github.io/bart/)


Tested Configuration:
---------------------
* Windows 10 Enterprise (v1803) and Subsystem for Linux (Ubuntu 16.04) for simulation, MRF recon, analysis, and figure generation
* Linux (Ubuntu 18.04) for MRF pEPG and ssEPG dictionary generation
* MATLAB R2017b (v9.3.0) in WSL used for image reconstruction and subsequent T1 and T2 fitting
* MATLAB R2018a/b (v9.4.0) in Windows 10 and Linux used for simulation, analysis, figure generation, and ssEPG and pEPG dictionary generation
* Berkeley Advanced Reconstruction Toolbox (v0.4.01)
* NVIDIA CUDA compiler (v9.2.148)

* 12-core Intel Xeon CPU E5-2687W v4
* 128 GB RAM

Installation Options:
---------------------
* Download the zip file of this repository to your local machine
* OR clone the git repository to your local machine
* OR fork to your own repository and then clone the git repository to your local machine


* Download contributing sample density compensation code (https://zenodo.org/record/401057/files/sdc3_nrz_11aug.zip) into `./contrib/`
* Compile necessary mex files in respective contributing code directories:
    -(e.g. `mex -compatibleArrayDims sdc3_MAT.c`)
* Compile ptx file in `./src/`
    -(e.g. `nvcc -ptx ssepg_functions.cu`)
* Download contributing code for figure capture:
    -https://github.com/altmany/export_fig/archive/f0af704d84608f5a69e3de82581869e7b6161d4f.zip

Usage:
------
* Download input data (see Releases)
* Go to `./scripts/` and modify the setup parameters in the first section of `run_all.m`:
    -modify `bart_path` to point to BART path
    -modify other paths to other contributing code (if different from default)
* Run `run_all.m`
* Alternatively, select the scripts of interest to run from `run_all.m`

* Code is not optimized for speed, the total processing time ~days if all processing is performed
* Figure generation quality is display dependent

Folder Structure:
--------

* `./scripts/` - contains all scripts for the different experiments
* `./src/` - contains all functions for executing the scripts (less the contributing code and BART)
* `./contrib/` - the downloaded contributing code
* `./data_in/` - the data input directory
* `./data_out/` - the reconstruction and processing output directory
* `./figures/` - the figure output directory

