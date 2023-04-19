# MMTPy: A Python Implementation of the Microbiome Modeling Toolbox (MMT) MgPipe Pipeline

Original Author: [Izzy Goodchild-Michelman, Zomorrodi Lab](https://github.com/zomorrodilab/izzy-gm/commits?author=zomorrodilab)

Original MATLAB implementation: [Microbiome Modeling Toolbox 2.0 by Heinken et al. (2022)](https://pubmed.ncbi.nlm.nih.gov/35157025/)

Latest version updated [here](https://github.com/kevinliu-bmb/MMTPy) by Kevin Liu

## Changes

### Bug Fixes

* Fixed error with reading in VMH diet with tab-delimiter.
* Fixed diet bounds being set as positive.

### Improvements

* Added optimization function ```opt_comm_gem.py```.
* Combined core pipeline into single file ```compy.py```.
* Added clearer printed checkpoints when running pipeline.
* Removed string manipulation redundancies.
* Cleaned comments and syntactic redundancies.
* Follows PEP 8 styling guidelines.
* Parallel processes for optimizing each multi-species model (experimental).

## Pending Implementations

* Automatic incorporation of metabolites and reactions needed for AGORA 2.01 model growth (emulating adaptVMHDietToAGORA() by [Heinken et al. (2022)](https://pubmed.ncbi.nlm.nih.gov/35157025/)).
* Parallelized model construction.
* Implementation of different optimization solvers. Currently only uses Gurobi.
