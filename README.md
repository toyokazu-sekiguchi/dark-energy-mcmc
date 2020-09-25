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
  - `wtype`: Integer specifying a parameterization of dark energy equation of state (EoS). Currently the following parameterizations are supported.
    - `wtype = 0`: The CPL parameterization $w(a) = w_0+w_a(1-a)$.
    - `100> wtype >= 10`: Binned w(a). Suppose i and j are the tens and ones places of `wtype` (i.e. `wtype` = 10*i+j), $1>=a>0$ is divided uniformly in $i^\sqrt{a}$ into $j+1$ points. That is, provided fixed $j$, larger $i$ gives a finer binning around $a=1$.
  - `w[0]`, `w[1]`, ...: EoS parameters. 
    - `w[0]`=$w_0$ and `w[1]`=$w_a$ when `wtype = 0`.
    - Array (`w[0]`,w[1],..., `w[j]`) gives the binned $w(a)$. `w[0]` (`w[j]`) corresponds to the bin with largest (smallest) $a$.    
  - `nnu`: Effective number of neutrinos. The total number of neutrinos are enhanced by this factor (temperature is fixed to the standard value i.e. $T_\nu = (4/11)^{1/3} T_\gamma$.
  - `mnu`: Sum of neutrino mass in units of eV.
  - `neutrino_hierarchy`: Flag for neutrino mass hierarchy. 1 for normal, 0 for degenerate and -1 for inverted ones.
* [LIKELIHOODS]
  - `use_BAO`,`use_H0`, `use_CMB`, `use_SNeIa`: Flags for whether data is incorporated in likelihood calculation. They should be either `true` of `false`.
  - `data_BAO`: Choice of BAO data; see `data/bao.dataset` for data available.
* [MCMC]
  - `obh2`, `odmh2`, `odeh2`,`w[0]`,..., `nnu`, `mnu`: When each parameter is varied in the parameter estimation, four numbers should be given in order: lower limit, upper limit, slope of prior distribution, initial fluctuations. When left as blank, corresponding parameter is fixed to the fiducial value. Commas ',' should be used to separate each item. Prior distribution of a each parameter $x$ is assumed to be in proportional to $\theta(x-x_{\rm min})\theta(x_{\rm max}-x)x^{n_x}$, where $x_{\rm min}$, $x_{\rm max}$, $n_x$ are the first three components in each line. 
  - `nwalkers`: Number of walkers in affine invariant MCMC sampler. This should be at least twice the number of varied parameters.
  - `nsteps`: Number of steps for MCMC analysis
  - `parallel`: If `true`, parallelization is implemented in the MCMC calculation.

### Role of each python file:
* `const.py`: Definition of units and constants
* `background.py`: Calculation of cosmological background evolution. 
* `likelihoods.py`: Calculation of likelihood function incorporating recent BAO (arXiv:1607.03155, arXiv:1801.03062, arXiv:1702.00176), direct Hubble measurement (arXiv:2001.03624) and CMB $\theta_*$ (arXiv:1807.06209).
* `SN.py`: Likelihood code for type Ia supernovae, taken from CosmoMC (https://github.com/cmbant/CosmoMC).
* `mcmc.py`: MCMC analysis based on Affine Invariant MCMC sampler (emcee). Parallelization is supported based on the multiprocessing python module. Restart functionarity is supported.
* `driver.py`: Main function.

### Files in `/data` directories:
* `bbn.dat`: Lookup table of BBN $Y_p$.
* `Pantheon/*`: Pantheon SNeIa data.
* `jla*`: JLA SNeIa data.

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
* April 8th, 2020
  - Initial release.
* September 25th, 2020
  - BAO likelihood function is replaced; now based on Alam et al. (arXiv:2007.08991).

# To-do list
- [ ] Visualization of reconstructed EoS as function of $a$.
