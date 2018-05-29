# BayeZirChron.c

C version of Bayesian zircon eruption age estimation code

#### Installation from command line

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
Compiling and running this parallel version additionally requires a working installation of MPI (either [Open MPI](https://www.open-mpi.org) or [MPICH](https://www.mpich.org)) A sample batchfile is provided in the example folder: [runTest.pbs](examples/synthetic\ dataset\ tests/runTest.pbs)