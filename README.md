# Description
Code for MCMC parameter estimation of dark energy models from cosmographical data.

# Prerequisites (vesion usued in development)
* Python3 (3.7)
* C compiler (gcc Apple LLVM version 10.0.1)
* SWIG (3.0.12)

## Python modules
* numpy (1.18.1)
* scipy (1.4.1)
* configparser (4.0.2)
* emcee (3.0.2)
* tqdm (4.43.0)
* h5py (2.10.0)
* getdist (1.1.0)

# Installation
Compilation is required in order to build a python interface for HyRec (written in C).
1. Git clone source file via `git clone https://github.com/toyokazu-sekiguchi/dark-energy-mcmc.git`.
2. Go to `dark-energy-mcmc/HyRec`.
3. Edit `INC_PY` and `LIB_PY` in `Makefile` appropriately, so that `Python.h` and `libpython3.*` are incorporated correctly. Then `make pyrec`.
4. Go back to the parent directory. 

# Usage
There are two stages in the analysis: MCMC and postprocessing.

## Stage 1: MCMC (including calculation for fiducial model)

### Basic usage:
`python3 driver.py params.ini`

This first calculates evolution in the fiducial model and MCMC run afterwords. MCMC chains are written in HDF5 format. MCMC tries to restart from previous chains if they exist.

### Description of input parameter file
`params.ini` specifies a variety of parameters and consists of five sections:
* [OUTPUT]
  - `root`: Chacters specifying output prefix.
* [COSMOLOGY]
  - `obh2`, `odmh2`, `odeh2`: Density parameters $\Omega_i h^2$ of baryon, dark matter (assuming no decay) and dark energy.
  - `wtype`: This specifies a dark energy parameterization. 
    - aaa
    - bbb
w[0] = -1
w[1] = -1
w[2] = -1
w[3] = -1
  - `nnu`: Effective number of neutrinos. The total number of neutrinos are enhanced by this factor (temperature is fixed to the standard value i.e. $T_\nu = (4/11)^{1/3} T_\gamma$.
  - `mnu`: Sum of neutrino mass in units of eV.
  - `neutrino_hierarchy`: Flag for neutrino mass hierarchy. 1 for normal, 0 for degenerate and -1 for inverted ones.
* [DDM SETUP] 
Background evolution is computed based on iteration method using the following parameters.
  - `num_a`: Number of scale factor bin. Empirically `30` is recommended.  
  - `max_it`: Maximum iteration number. This is not relevant because usually convergence is achieved within one iteration.
  - `tol`: Tolerance parameter for the error from the true evolution measured by the current energy density of massless decay-product. Setting to `1e-10` works.
* [LIKELIHOODS]
  - `use_BAO`,`use_H0`, `use_CMB`: Flags for whether data is incorporated in likelihood calculation. They should be either `true` of `false`
* [MCMC]
  - `ob`, `odm`, `ol`, `decay_rate`, `mratio`, `nnu`, `mnu`: When each parameter is varied in the parameter estimation, four numbers should be given in order: lower limit, upper limit, slope of prior distribution, initial fluctuations. When left as blank, corresponding parameter is fixed to the fiducial value. Commas ',' should be used to separate each item. Prior distribution of a each parameter $x$ is assumed to be in proportional to $\theta(x-x_{\rm min})\theta(x_{\rm max}-x)x^{n_x}$, where $x_{\rm min}$, $x_{\rm max}$, $n_x$ are the first three components in each line. 
  - `nwalkers`: Number of walkers in affine invariant MCMC sampler. This should be at least twice the number of varied parameters.
  - `nsteps`: Number of steps for MCMC analysis
  - `parallel`: If `true`, parallelization is implemented in the MCMC calculation.

### Role of each python file:
* `const.py`: Definition of units and constants
* `mdd.py`: Calculation of cosmological background evolution. 
* `likelihoods.py`: Calculation of likelihood function incorporating recent BAO (arXiv:1607.03155, arXiv:1801.03062, arXiv:1702.00176), direct Hubble measurement (arXiv:2001.03624) and CMB $\theta_*$ (arXiv:1807.06209).
* `mcmc.py`: MCMC analysis based on Affine Invariant MCMC sampler (emcee). Parallelization is supported based on the multiprocessing python module. Restart functionarity is supported.
* `driver.py`: Main function.

## Stage 2: Postprocessing

### Basic usage:
`python3 post.py dist.ini`

This analyses MCMC chain(s) produced in Step 1 and obtain parameter constraints as well as triangle plot of posteior distributions.

### Description of parameter files
`dist.ini` specifies 
* [POSTPROCESS]
  - `postroot`: This specifies output prefix.
  - `paramnames`: Array of parameter names varied in chains. Commas separate items.
  - `paramlabels`: Array of parameter labels in LaTeX format. They are adopted in plotting. Commas separate items.
  - `dparamnames`: Array of derived parameter names in chains. Commas separate items. A Derived paraemters are either of `H0`, `Age`, `odm0`, `odm1`, `odm2`.
  - `dparamnames`: Array of derived parameter labels in chains. They are adopted in plotting. Commas separate items.
  - `chains`: Array of chain file(s) to be analysed. Commas separate items.
  - `chainlabels`: Array of chain label(s). They are adopted in plotting. Commas separate items.
  
# Notes
* Flatness is assumed.
* Neutrinos are assumed to consist of three mass eigenstates.
* 4He abundance $Y_p(\omega_b, N_\nu)$ is fitted with a look-up table in `BBN.dat`, which is taken from CLASS, which are originally obtained using the PArthENoPE code (http://parthenope.na.infn.it).
* Recombination history is computed based on HyReC (https://pages.jh.edu/~yalihai1/hyrec/hyrec.html). Hyrec in our code is modified from the original one so that massive neutrinos are incorporated and interface to Python is realized by SWIG.

# Version history
* March 30th, 2020
  - Some derived parameters are calculated and written in chains at run-time of the MCMC stage. Slopes of prior distributions of parameters can be now specified.
  - Postprocessing is updated. Some minor bugs are fixed. Derived parameters are now incorporated.
* March 4nd, 2020
  - Restart functionarity is now supported.
  - HyRec wrapper is modified. In the previous version, segmentation faults occur when HyRec is located on a path which includes symbolic links. All the files for look-up tables are now referred by their physical pathes.
* March 2nd, 2020
  - Initial release.

# To-do list 
- [ ] Multiple plots with different sets of (primary and derived) parameters in different scales (e.g. linear, log, etc.).
