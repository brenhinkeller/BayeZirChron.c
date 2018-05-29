# BayeZirChron.c

C version of Bayesian zircon eruption age estimation code

#### Installation from command line

Installation from the command line requires a working C compiler. The default makefile assumes GCC is available. 

```bash
# Download
git clone https://github.com/brenhinkeller/BayeZirChron.c.git

# Move to folder containing source code
cd BayeZirChron.c/src/ 	

# Compile	
make serial 		
```

To compile parallel code used for synthetic distribution tests, additionally run:

```bash
make parallel
```
or 

```bash
mpicc -std=c11 -O3 -o tzircrystestaccuracy tzircrystestaccuracy.c
```
Compiling and running this parallel version additionally requires a working installation of MPI.