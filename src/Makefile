serial:
	gcc -std=c11 -O3 -o tzircrystmetropolis tzircrystmetropolis.c
	gcc -std=c11 -O3 -o tzircrystmetropolis1sigma tzircrystmetropolis1sigma.c
	gcc -std=c11 -O3 -o tzircrystimage tzircrystimage.c

parallel:
	mpicc -std=c11 -O3 -o tzircrystestaccuracy tzircrystestaccuracy.c
	mpicc -std=c11 -O3 -o tzircrystestaccuracyuniform tzircrystestaccuracyuniform.c
	mpicc -std=c11 -O3 -o tzircrystestaccuracydist tzircrystestaccuracydist.c

all: 	serial parallel
