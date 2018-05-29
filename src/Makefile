serial:
	gcc -std=c11 -O3 -o tzircrystmetropolis tzircrystmetropolis.c
	gcc -std=c11 -O3 -o tzircrystmetropolis1sigma tzircrystmetropolis1sigma.c
	gcc -std=c11 -O3 -o tzircrystimage tzircrystimage.c

parallel:
	mpicc -std=c11 -O3 -o tzircrystestaccuracy tzircrystestaccuracy.c


all: 	serial parallel

