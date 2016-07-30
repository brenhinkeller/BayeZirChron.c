/******************************************************************************
 * FILE: tzircrystimage.c
 * DESCRIPTION: Creates a 2-d matrix for the likelihood of zircon crystallization
 * models as a function of initial and final zircon crystallization time
 * AUTHOR: C. Brenhin Keller
 ******************************************************************************/


#include <stdio.h>
#include <math.h>
#include <time.h>
#include "arrays.h"

#include "pcg_variants.h"
#include "gauss.h"


void drawFromDistribution(pcg32_random_t* rng, const double* dist, const uint32_t distrows, double* x, const uint32_t xrows){
	const double dist_yrandmax = UINT32_MAX / maxArray(dist,distrows);
	const double dist_xscale = (double)(distrows-1); 
	const double dist_xrandmax = UINT32_MAX / dist_xscale;
	double rx, ry, y;
 
	for (int i=0; i<xrows; i++){
		x[i] = -1;
		while (x[i] == -1){
			// Pick random x value
			rx = pcg32_random_r(rng) / dist_xrandmax;
			// Interpolate corresponding distribution value
			y = interp1i(dist,rx);
			// See if x value is accepted
			ry = pcg32_random_r(rng) / dist_yrandmax;
			if (y > ry){
				x[i] = rx / dist_xscale;

			}
		}
	}
}

void generateSyntheticZirconDataset(pcg32_random_t* rng, const double* dist, const uint32_t distrows, const double tmin, const double tmax, const double* uncert, double* synzirc, const uint32_t datarows){
	double dt = fabs(tmax - tmin);
	double r;

	drawFromDistribution(rng, dist, distrows, synzirc, datarows);
	sort_doubles_descending(synzirc, datarows);
	for(int i=0; i<datarows; i++){
		r = pcg_gaussian_ziggurat(rng, uncert[i]);
		synzirc[i] = (1-synzirc[i]) * dt + tmin + r;

	}
}

double compareZirconPopulations(const double* data, const double* uncert, const double* synzirc, const uint32_t rows){
	double loglikelihood = 0;
	for (int i=0; i<rows; i++){
		loglikelihood += log10( 1/(uncert[i] * sqrt(2*M_PI)) * exp( - (synzirc[i]-data[i])*(synzirc[i]-data[i]) / (2*uncert[i]*uncert[i]) ));
	}
	return loglikelihood / (double)rows;
}


double testAgeModel(pcg32_random_t* rng, const double* dist, const uint32_t distrows, const double* data, const double* uncert, double* synzirc, const uint32_t datarows, const uint32_t nsims, const double tmin, const double tmax){
	double loglikelihood = 0;
	for (int i=0; i<nsims; i++){
		generateSyntheticZirconDataset(rng, dist, distrows, tmin, tmax, uncert, synzirc, datarows);
		sort_doubles(synzirc, datarows);
		loglikelihood += compareZirconPopulations(data, uncert, synzirc, datarows);
	}
	return loglikelihood / (double)nsims;
}

double checkZirconLikelihood(const double* restrict dist, const uint32_t distrows, const double* restrict data, const double* restrict uncert, const uint32_t datarows, const double tmin, const double tmax){
	double distx, likelihood, loglikelihood = 0;
	const double dt = fabs(tmax-tmin);

	for (int j=0; j<datarows; j++){
		likelihood = 0;
		for (int i=0; i<distrows; i++){
			distx = tmax - dt*i/(distrows-1);
			likelihood += dist[i] / (distrows * uncert[j] * sqrt(2*M_PI)) * exp( -(distx-data[j])*(distx-data[j]) / (2*uncert[j]*uncert[j]) );
		}
		loglikelihood += log10(likelihood);
	}
	return loglikelihood;
}

int main(int argc, char **argv){
	uint32_t distrows, distcolumns;
	uint32_t datarows, datacolumns;
	uint32_t i, j, k;

	//Check input arguments
	if (argc != 4) {
		fprintf(stderr,"USAGE: %s <nsims> <distribution.tsv> <zircondata.tsv>\n", argv[0]);
		exit(1);
	}

	// Get number of simulations per MPI task from command-line argument
	const uint32_t nsims = (uint32_t)atoi(argv[1]);

	// Import data	
	const double* dist = csvparseflat(argv[2],'\t', &distrows, &distcolumns);	
	const double* data = csvparseflat(argv[3],'\t', &datarows, &datacolumns);
	const double tmin = minArray(data, datarows);
	const double tmax = maxArray(data, datarows);
	const double dt = tmax - tmin;
	double* synzirc = malloc(datarows * sizeof(double));
	for (i=0; i<datarows; i++) {synzirc[i]=0;}

//	printf("distribution size: %i, %i\n", distrows, distcolumns);
//	printf("data size: %i, %i\n", datarows, datacolumns);

	pcg32_random_t rng;
	pcg32_srandom_r(&rng,time(NULL), clock());


	double tl, tu;


	printf("NaN\t");
	for (tu = tmax-dt; tu < tmax+dt;  tu += dt/50.0){
		printf("%g\t", tu);
	}
	printf("\n");

	// Print matrix image, with x and y scales
	for (tl = tmin-dt; tl < tmin+dt;  tl += dt/50.0){
		printf("%g\t", tl);

		for (tu = tmax-dt; tu < tmax+dt;  tu += dt/50.0){
			if(tu>tl){
//				printf("%g\t", testAgeModel(&rng, dist, distrows, data, &data[datarows*1], synzirc, datarows, nsims, tl, tu));
				printf("%g\t", checkZirconLikelihood(dist, distrows, data, &data[datarows*1], datarows, tl, tu));
			} else {
				printf("NaN\t");
			}
		}
		printf("\n");
	}

return 0;
}
