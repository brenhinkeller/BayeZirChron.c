/******************************************************************************
 * FILE: tzircrystmetroplis.c
 * DESCRIPTION: Runs metropolis walker to find optima of zircon age distribution
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
		loglikelihood += log10( 1 / (uncert[i] * sqrt(2*M_PI)) * exp( - (synzirc[i]-data[i])*(synzirc[i]-data[i]) / (2*uncert[i]*uncert[i]) ));
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


int main(int argc, char **argv){
	uint32_t distrows, distcolumns;
	uint32_t datarows, datacolumns;
	uint32_t i, j, k;

	//Check input arguments
	if (argc != 5) {
		fprintf(stderr,"USAGE: %s <nsims> <nsteps> <distribution.tsv> <zircondata.tsv>\n", argv[0]);
		exit(1);
	}

	// Get number of simulations per MPI task from command-line argument
	const uint32_t nsims = (uint32_t)abs(atoi(argv[1]));
	// Get number of steps for metropolis walker to take	
	const uint32_t nsteps = (uint32_t)abs(atoi(argv[2]));

	// Import data	
	const double* dist = csvparseflat(argv[3],'\t', &distrows, &distcolumns);	
	const double* data = csvparseflat(argv[4],'\t', &datarows, &datacolumns);
	const double tmin_obs = minArray(data, datarows);
	const double tmax_obs = maxArray(data, datarows);
	const double dt = tmax_obs - tmin_obs;
	double* synzirc = malloc(datarows * sizeof(double));
	for (i=0; i<datarows; i++) {synzirc[i]=0;}

//	printf("distribution size: %i, %i\n", distrows, distcolumns);
//	printf("data size: %i, %i\n", datarows, datacolumns);

	pcg32_random_t rng;
	pcg32_srandom_r(&rng,time(NULL), clock());


	double tmin, tmax, theta, tmin_proposed, tmax_proposed, theta_proposed;
	double r, tmin_step=dt/10, tmax_step=dt/10;

	tmin = tmin_obs;
	tmax = tmax_obs;
	tmin_proposed = tmin_obs;
	tmax_proposed = tmax_obs;
	theta = testAgeModel(&rng, dist, distrows, data, &data[datarows*1], synzirc, datarows, nsims, tmin_proposed, tmax_proposed);

	double p = 10;

	for (i=0; i<nsteps; i++){
		
		// Uniformly adjust either the upper or lower bound age
		r = pcg32_random_r(&rng)/(double)UINT32_MAX;
		if (r<0.5){
			tmin_proposed = tmin + pcg_gaussian_ziggurat(&rng, tmin_step);
			tmax_proposed = tmax;
		} else {
			tmin_proposed = tmin;
			tmax_proposed = tmax + pcg_gaussian_ziggurat(&rng, tmax_step);
		}

		// Calculate log likelihood for new proposal
		theta_proposed = testAgeModel(&rng, dist, distrows, data, &data[datarows*1], synzirc, datarows, nsims, tmin_proposed, tmax_proposed);

		printf("%f\t%f\t", theta_proposed, theta);

		// Decide to accept or reject the proposal
		r = pcg32_random_r(&rng)/(double)UINT32_MAX;
		if (r < pow(10,theta_proposed-theta) && tmin_proposed < tmax_proposed){
			tmin = tmin_proposed;
			tmax = tmax_proposed;
			theta = theta_proposed;
		}

		printf("%f\t%f\n", tmin, tmax);
	}

return 0;
}
