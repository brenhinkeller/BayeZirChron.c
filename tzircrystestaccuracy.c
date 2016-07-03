/******************************************************************************
 * FILE: tzircrystestaccuracy.c
 * DESCRIPTION: Creates multiple synthetic zircon datasets of varying N, and 
 * compares accuracy of weighted mean, youngest-zircon, and Bayesian 
 * metropolis walker estimates for time of final zircon crystallization.
 * AUTHOR: C. Brenhin Keller
 ******************************************************************************/


#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "arrays.h"

#include "pcg_variants.h"
#include "gauss.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ROOT 0

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
	double tm = min(tmax,tmin);
	double r;

	drawFromDistribution(rng, dist, distrows, synzirc, datarows);
	sort_doubles_descending(synzirc, datarows);
	for(int i=0; i<datarows; i++){
		r = pcg_gaussian_ziggurat(rng, uncert[i]);
		synzirc[i] = (1-synzirc[i]) * dt + tm + r;

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


int findMetropolisEstimate(pcg32_random_t* rng, const double* dist, const uint32_t distrows, const double* data, const double* uncert, double* synzirc, const uint32_t datarows, const uint32_t nsims, const uint32_t nsteps, double* const restrict mu, double* const restrict sigma){
	const uint32_t burnin = nsteps/10;
	const double tmin_obs = minArray(data, datarows);
	const double tmax_obs = maxArray(data, datarows);
	const double dt = tmax_obs - tmin_obs + uncert[0] + uncert[datarows-1];


	double tmin, tmax, theta, tmin_proposed, tmax_proposed, theta_proposed;
	double r, tmin_step=dt/(double)datarows, tmax_step=dt/(double)datarows;
	uint32_t i;

	double tmins[nsteps+burnin];

	
	tmin = tmin_obs;
	tmax = tmax_obs;
	tmin_proposed = tmin_obs;
	tmax_proposed = tmax_obs;
	theta = testAgeModel(rng, dist, distrows, data, uncert, synzirc, datarows, nsims, tmin_proposed, tmax_proposed);


	for (i=0; i<nsteps+burnin; i++){

//		// Uniformly adjust either the upper or lower bound age
//		r = pcg32_random_r(rng)/(double)UINT32_MAX;
//		if (r<0.5){
//			tmin_proposed = tmin + pcg_gaussian_ziggurat(rng, tmin_step);
//			if (tmin_proposed > tmax)
//				tmax_proposed = tmax;
//		} else {
//			tmin_proposed = tmin;
//			tmax_proposed = tmax + pcg_gaussian_ziggurat(rng, tmax_step);
//		}


		// Uniformly adjust upper and lower bound age
		tmin_proposed = tmin + pcg_gaussian_ziggurat(rng, tmin_step);
		tmax_proposed = tmax + pcg_gaussian_ziggurat(rng, tmax_step);

		// Swap upper and lower bound if they get reversed
		if (tmin_proposed>tmax_proposed){
			r = tmin_proposed;
			tmin_proposed = tmax_proposed;
			tmax_proposed = r;
		}

		// Calculate log likelihood for new proposal
		theta_proposed = testAgeModel(rng, dist, distrows, data, uncert, synzirc, datarows, nsims, tmin_proposed, tmax_proposed);

		// Decide to accept or reject the proposal
		r = pcg32_random_r(rng)/(double)UINT32_MAX;
		if (r < pow(10,theta_proposed-theta)){
			if (tmin_proposed != tmin){
				tmin_step = fabs(tmin_proposed-tmin)*2.718;
			} 
			if (tmax_proposed != tmax){
				tmax_step = fabs(tmax_proposed-tmax)*2.718;
			}
			tmin = tmin_proposed;
			tmax = tmax_proposed;
			theta = theta_proposed;
		}

		tmins[i] = tmin;
	}

return Offset_nanstd(&tmins[burnin], nsteps, mu, sigma);
}


int main(int argc, char **argv){
	uint32_t distrows, distcolumns;
	uint32_t N, i, j, k;
	int world_size, world_rank, rc;

	//Check input arguments
	if (argc != 4) {
		fprintf(stderr,"USAGE: %s <nsims> <nsteps> <distribution.tsv>\n", argv[0]);
		exit(1);
	}

	// Start MPI
	rc = MPI_Init(&argc,&argv); 
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n"); MPI_Abort(MPI_COMM_WORLD, rc);
	}

	// Get world size (number of MPI processes) and world rank (# of this process)
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

	if (world_rank==ROOT){
		printf("N\tWeighted Mean\tMSWD Test\tYoungest Zircon\tMetropolis Estimate\tMetropolis+MSWD\n");
	}

	// Get number of simulations per MPI task from command-line argument
	const uint32_t nsims = (uint32_t)abs(atoi(argv[1]));
	// Get number of steps for metropolis walker to take	
	const uint32_t nsteps = (uint32_t)abs(atoi(argv[2]));

	// Import data	
	const double* dist = csvparseflat(argv[3],'\t', &distrows, &distcolumns);	

	// Declare various variables
	pcg32_random_t rng;
	pcg32_srandom_r(&rng,time(NULL), world_rank);

	const uint32_t Nmin = 1;
	const uint32_t Nmax = 1024;
	uint32_t Ns[] = {1,2,3,4,6,8,11,16,23,32,45,64,91,128,181,256,362,512,724,1024};
	uint32_t nNs = 20;

	double tmin_true = 100;
	double relagerange = 0.1/100.0;
	double reluncert = 0.1/100.0;
	double absuncert = reluncert*tmin_true;
	double simspertask = 2;
	double tmax_true = tmin_true*(1+relagerange);


	double* synzirc = malloc(Nmax * sizeof(double));
	double* data = malloc(Nmax * sizeof(double));
	double* uncert = malloc(Nmax * sizeof(double));
	double wx, wsigma, mswd; // For weighted mean
	double tmin_metropolis_mswd, tmin_metropolis_est, tmin_obs, tmin_mswd_test; // Means
	double tmin_metropolis_mswd_sigma, tmin_metropolis_est_sigma, tmin_obs_sigma, tmin_mswd_test_sigma; // Standard deviations

	for (i=0; i<nNs; i++){
		N = Ns[i];

		for (j=0; j<N; j++){
			uncert[j] = tmin_true * reluncert;
		}

		for (j=0; j<simspertask; j++){
			generateSyntheticZirconDataset(&rng, dist, distrows, tmin_true, tmax_true, uncert, data, N);

			findMetropolisEstimate(&rng, dist, distrows, data, uncert, synzirc, N, nsims, nsteps, &tmin_metropolis_est, &tmin_metropolis_est_sigma);
			
			sort_doubles(data, N);
			tmin_mswd_test = data[0];
			tmin_mswd_test_sigma = absuncert;
			for (k=2; k<N+1; k++){
				wmean(data, uncert, k, &wx, &wsigma, &mswd);
				if (mswd > 1.0 + 2.0 * sqrt(2.0/(double)(k-1))){
					break;
				}
				tmin_mswd_test = wx;
				tmin_mswd_test_sigma = wsigma;
			}

			wmean(data, uncert, N, &wx, &wsigma, &mswd);

			tmin_obs = minArray(data, N);
			tmin_obs_sigma = absuncert;


			if (mswd < 1.0 + 2.0 * sqrt(2.0/(double)(N-1))){
				tmin_metropolis_mswd = wx;
				tmin_metropolis_mswd_sigma = wsigma;
			} else {
				tmin_metropolis_mswd = tmin_metropolis_est;
				tmin_metropolis_mswd_sigma = tmin_metropolis_est_sigma;
			}

			printf("%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", N, fabs(wx-tmin_true)/absuncert, fabs(tmin_mswd_test-tmin_true)/absuncert, fabs(tmin_obs-tmin_true)/absuncert, fabs(tmin_metropolis_est-tmin_true)/absuncert, fabs(tmin_metropolis_mswd-tmin_true)/absuncert, fabs(wx-tmin_true)/wsigma, fabs(tmin_mswd_test-tmin_true)/tmin_mswd_test_sigma, fabs(tmin_obs-tmin_true)/tmin_obs_sigma, fabs(tmin_metropolis_est-tmin_true)/tmin_metropolis_est_sigma, fabs(tmin_metropolis_mswd-tmin_true)/tmin_metropolis_mswd_sigma);
		}
	}
	
MPI_Finalize();
return 0;
}
