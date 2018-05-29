# BayeZirChron.c

C version of Bayesian zircon eruption age estimation code

## Installation from command line

Installation from the command line requires a working C compiler. The default [makefile](src/Makefile) assumes [gcc](https://gcc.gnu.org) (or an alias) is available. 
On linux/unix/bsd this is likely already true; on Mac OS the necessary tools for compiling C source can be installed by typing `xcode-select --install` at the command line.

```bash
# Download
git clone https://github.com/brenhinkeller/BayeZirChron.c.git

# Move to folder containing source code
cd BayeZirChron.c/src/

# Compile
make serial
```

To compile the parallel code used for synthetic distribution tests, additionally run:

```bash
make parallel
```
or 

```bash
mpicc -std=c11 -O3 -o tzircrystestaccuracy tzircrystestaccuracy.c
```
Compiling and running this parallel version additionally requires a working installation of MPI (either [Open MPI](https://www.open-mpi.org) or [MPICH](https://www.mpich.org)) A sample batchfile is provided in the example folder: [runTest.pbs](examples/synthetic%20dataset%20tests/runTest.pbs)

## Usage

A range of [examples](examples/), including the application Bayesian zircon eruption age estimation code to literature datasets, is provided. A Matlab script to run all literature examples is provided in [examples/literature dataset tests/RunLiteratureExamples.m](examples/literature%20dataset%20tests/RunLiteratureExamples.m)

Basic command-line usage follows the pattern:
```bash
tzircrystmetropolis <nsteps> distribution.tsv  sample.tsv > output.tsv
```
where <nsteps> is the length of Markov chain to run, `distribution.tsv` is an ascii file describing a zircon saturation distribution, `sample.tsv` is a  tab or comma-separated file containing a list of individual zircon ages and 2-sigma uncertainties; the resulting stationary distributions will be written to `output.tsv`. For example
```bash
tzircrystmetropolis 10000 MeltsTZircDistribtuion.tsv  sample.tsv > output.tsv
```

The results of synthetic dataset tests, which compare traditional weighted-mean and youngest-zircon interpretations to Bayesian eruption age estimates (using several different crystallization distributions)  are provided in [examples/synthetic dataset tests/](examples/synthetic%20dataset%20tests/). Figures can be re-plotted using a Matlab script [tzircrystestaccuracyPlots.m](examples/synthetic%20dataset%20tests/tzircrystestaccuracyPlots.m).

To reproduce the datafiles provided in this folder, compile the parallel code [tzircrystestaccuracy.c](src/tzircrystestaccuracy.c) as described above (using mpicc), and run using

```bash
mpiexec -np <number-of-tasks> ./tzircrystestaccuracy  <sims-per-task> <nsteps> <dt/sigma>  Distribution.tsv  > results.tsv
```
where <number-of-tasks> is the number of MPI tasks to run (typically you want this to be equal to the number of CPU cores or hardware threads you are running on),  <sims-per-task> is the number of simulations (at each N) to run per MPI task, <nsteps> is the length of Markov chain to run, <dt/sigma> is the crystallization timescale in units of sigma (analytical uncertainty), pulling synthetic data from a distribution specified in an ascii file `Distribution.tsv	
for example:
```bash
mpiexec -np 16 ./tzircrystestaccuracy 4 10000 1 MeltsTZircDistribtuion.tsv > eruptionestimates1.tsv
```
to run on 16 cores with 4 simulations per taks, each 10000 MCMC steps long, with a dt/sigma of 1 and using the crystallization distribution found in [MeltsTZircDistribtuion.tsv](distributions/MeltsTZircDistribtuion.tsv).  To run on a cluster, you will need a batch file suited to your cluster's workload manager. For instance, the example batchfile [runTest.pbs](examples/synthetic%20dataset%20tests/runTest.pbs):

```bash
#!/bin/bash
#PBS -l nodes=20:ppn=16,walltime=00:40:00

module load openmpi
mpiexec ./tzircrystestaccuracy 4 10000 1 MeltsTZircDistribtuion.tsv > eruptionestimates1.tsv
```

runs 4 simulations per task on 20 nodes with 16 cores each, for a total of 1280 simulations at each N



