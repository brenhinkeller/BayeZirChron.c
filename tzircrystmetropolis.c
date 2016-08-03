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

double stirling_lgamma(const double x){
	// Stirling's Approximation for gamma function ( generalized log[(n-1)!] )
	return (x - 0.5)*log(x) - x + 0.5*log(2*M_PI) + 1/(12*x) + 1/(360 * x * x * x); // + 1/(1260 * x * x * x * x * x) - ...
}

double checkZirconLikelihood(const double* restrict dist, const uint32_t distrows, const double* restrict data, const double* restrict uncert, const uint32_t datarows, const double tmin, const double tmax){
	double Zf, ix, wm, wsigma, mswd, distx, likelihood, loglikelihood = 0;
	const double dt = fabs(tmax-tmin);
	const double dist_xscale = (double)(distrows-1);
	const double dist_yave = nanmean(dist, distrows);

	for (int j=0; j<datarows; j++){
		// If possible, prevent aliasing problems by interpolation
		if ((uncert[j] < dt/dist_xscale) && data[j] > tmin && data[j] < tmax){
			// Find (double) index
			ix = (tmax - data[j])/dt*dist_xscale;
			// Interpolate corresponding distribution value
			likelihood = interp1i(dist,ix)/dt;
		// Otherwise, sum contributions from Gaussians at each point in distribution
		} else {
			likelihood = 0;
			for (int i=0; i<distrows; i++){
				distx = tmax - dt*i/(distrows-1);
				likelihood += dist[i] / (dist_yave * distrows * uncert[j] * sqrt(2*M_PI)) * exp( - (distx-data[j])*(distx-data[j]) / (2*uncert[j]*uncert[j]) );
			}
		}
		loglikelihood += log10(likelihood);
	}


	wmean(data, uncert, datarows, &wm, &wsigma, &mswd);
	if (datarows == 1){
		Zf = 1;
	} else if (mswd*sqrt(datarows)>1000){
		Zf = 0;
	} else {
		const double f = (double)datarows-1;
//		Zf = exp(f/2*log(f/2) - stirling_lgamma(f/2) + (f/2-1)*log(mswd) - f/2*mswd); // MSWD distribution from Wendt and Carl
		Zf = exp((f/2-1)*log(mswd) - f/2*(mswd-1)); // Height of MSWD distribution relative to height at mswd = 1;
	}

	return loglikelihood - log10(dt/wsigma) * Zf;
}



int main(int argc, char **argv){
	uint32_t distrows, distcolumns;
	uint32_t datarows, datacolumns;
	uint32_t i, j, k;

	//Check input arguments
	if (argc != 6) {
		fprintf(stderr,"USAGE: %s <nsims> <nsteps> <stepfactor> <distribution.tsv> <zircondata.tsv>\n", argv[0]);
		exit(1);
	}

	// Get number of simulations per MPI task from command-line argument
	const uint32_t nsims = (uint32_t)abs(atoi(argv[1]));
	// Get number of steps for metropolis walker to take	
	const uint32_t nsteps = (uint32_t)abs(atoi(argv[2]));
	// sigma = stepfactor * last-step	
	const double stepfactor = atof(argv[3]);


	// Import data	
	const double* dist = csvparseflat(argv[4],'\t', &distrows, &distcolumns);	
	const double* data = csvparseflat(argv[5],'\t', &datarows, &datacolumns);
	const double tmin_obs = minArray(data, datarows);
	const double tmax_obs = maxArray(data, datarows);
	const double dt = tmax_obs - tmin_obs + data[datarows*1+0] + data[datarows*1 + datarows-1];
	double* synzirc = malloc(datarows * sizeof(double));
	for (i=0; i<datarows; i++) {synzirc[i]=0;}

//	printf("distribution size: %i, %i\n", distrows, distcolumns);
//	printf("data size: %i, %i\n", datarows, datacolumns);

	pcg32_random_t rng;
	pcg32_srandom_r(&rng,time(NULL), clock());


	double tmin, tmax, theta, tmin_proposed, tmax_proposed, theta_proposed;
	double r, tmin_step=dt/(double)datarows, tmax_step=dt/(double)datarows;

	tmin = tmin_obs;
	tmax = tmax_obs;
	tmin_proposed = tmin_obs;
	tmax_proposed = tmax_obs;
//	theta = testAgeModel(&rng, dist, distrows, data, &data[datarows*1], synzirc, datarows, nsims, tmin_proposed, tmax_proposed);
	theta =  checkZirconLikelihood(dist, distrows, data, &data[datarows*1], datarows, tmin_proposed, tmax_proposed);

	double p = 10;

	int transitions=0;

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

		if (tmin_proposed>tmax_proposed){
			r = tmin_proposed;
			tmin_proposed = tmax_proposed;
			tmax_proposed = r;
		}

		// Calculate log likelihood for new proposal
//		theta_proposed = testAgeModel(&rng, dist, distrows, data, &data[datarows*1], synzirc, datarows, nsims, tmin_proposed, tmax_proposed);
		theta_proposed =  checkZirconLikelihood(dist, distrows, data, &data[datarows*1], datarows, tmin_proposed, tmax_proposed);

		// Decide to accept or reject the proposal
		r = pcg32_random_r(&rng)/(double)UINT32_MAX;
		if (r < pow(10,theta_proposed-theta)){// && tmin_proposed < tmax_proposed){ // don't let tmax go below tmin
			if (tmin_proposed != tmin){
				tmin_step = fabs(tmin_proposed-tmin)*stepfactor;
			} 
			if (tmax_proposed != tmax){
				tmax_step = fabs(tmax_proposed-tmax)*stepfactor;
			}

			tmin = tmin_proposed;
			tmax = tmax_proposed;
			theta = theta_proposed;
		}

		printf("%f\t%f\t%f\n", tmin, tmax, theta);
	}

return 0;
}
