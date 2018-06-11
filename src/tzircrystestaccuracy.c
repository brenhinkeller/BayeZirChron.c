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
		while (x[i]==-1){
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
	const double dt = fabs(tmax - tmin);
	const double tm = min(tmax,tmin);
	double r;

	drawFromDistribution(rng, dist, distrows, synzirc, datarows);
//	sort_doubles_descending(synzirc, datarows);
	for(int i=0; i<datarows; i++){
		r = pcg_gaussian_ziggurat(rng, uncert[i]);
		while (isnan(r)){
			r = pcg_gaussian_ziggurat(rng, uncert[i]);
		}

		synzirc[i] = (1-synzirc[i]) * dt + tm + r;
	}
}


double stirling_lgamma(const double x){
	if (x<1){
		return 0;
	} else {
		// Stirling's Approximation for gamma function ( generalized log[(n-1)!] )
		return (x - 0.5)*log(x) - x + 0.5*log(2*M_PI) + 1/(12*x) + 1/(360 * x * x * x); // + 1/(1260 * x * x * x * x * x) - ...
	}
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
			likelihood = interp1i(dist,ix)/(dt * dist_yave);
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

	awmean(data, uncert, datarows, &wm, &wsigma, &mswd);
	if (datarows == 1 || mswd<1){
		Zf = 1;
		mswd = 1;
	} else if (mswd*sqrt(datarows)>1000){
		Zf = 0;
	} else {
		const double f = (double)datarows-1;
		//Zf = exp(f/2*log(f/2) - stirling_lgamma(f/2) + (f/2-1)*log(mswd) - f/2*mswd); // Distribution of the MSWD for normally-distributed data, from Wendt and Carl 1991
		Zf = exp((f/2-1)*log(mswd) - f/2*(mswd-1)); // Height of MSWD distribution relative to height at mswd = 1; from Wendt and Carl 1991
	}

	// At low N (low sample numbers), favor the weighted mean interpretation at high Zf
 	// (MSWD close to 1) and the youngest-zircon interpretation at low Zf (MSWD far from one).
 	// This is helpful to for preventing instability at low N, with the funcional form
 	// used here derived from training against synthetic datasets.
	return loglikelihood  - ( log10((fabs(tmin - wm)+wsigma)/wsigma)*Zf + log10((fabs(tmax - wm)+wsigma)/wsigma)*Zf + log10((fabs(tmin - data[0])+uncert[0])/uncert[0])*(1-Zf) + log10((fabs(tmax - data[datarows-1])+uncert[datarows-1])/uncert[datarows-1])*(1-Zf) ) * (2/log10(1+datarows));
}



int findMetropolisEstimate(pcg32_random_t* rng, const double* dist, const uint32_t distrows, const double* data, const double* uncert, const uint32_t datarows, const uint32_t nsteps, double* const restrict mu, double* const restrict sigma){
	const uint32_t burnin = nsteps/2;
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
	theta =  checkZirconLikelihood(dist, distrows, data, uncert, datarows, tmin_proposed, tmax_proposed);

	for (i=0; i<nsteps+burnin; i++){

		// Uniformly adjust either the upper or lower bound age
		r = pcg32_random_r(rng)/(double)UINT32_MAX;
		if (r<0.5){
			tmin_proposed = tmin + pcg_gaussian_ziggurat(rng, tmin_step);
			tmax_proposed = tmax;
		} else {
			tmin_proposed = tmin;
			tmax_proposed = tmax + pcg_gaussian_ziggurat(rng, tmax_step);
		}

		// // Uniformly adjust upper and lower bound age
		// tmin_proposed = tmin + pcg_gaussian_ziggurat(rng, tmin_step);
		// tmax_proposed = tmax + pcg_gaussian_ziggurat(rng, tmax_step);

		// Swap upper and lower bound if they get reversed
		if (tmin_proposed>tmax_proposed){
			r = tmin_proposed;
			tmin_proposed = tmax_proposed;
			tmax_proposed = r;
		}

		// Calculate log likelihood for new proposal
		theta_proposed = checkZirconLikelihood(dist, distrows, data, uncert, datarows, tmin_proposed, tmax_proposed);

		// Decide to accept or reject the proposal
		r = pcg32_random_r(rng)/(double)UINT32_MAX;
		if (r < pow(10,theta_proposed-theta)){
			if (tmin_proposed != tmin){
				tmin_step = fabs(tmin_proposed-tmin)*2.9;
			}
			if (tmax_proposed != tmax){
				tmax_step = fabs(tmax_proposed-tmax)*2.9;
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
	if (argc != 5) {
		fprintf(stderr,"USAGE: %s <sims-per-task> <nsteps> <dt/sigma> <distribution.tsv>\n", argv[0]);
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
		printf("N\tMSWD\tWeighted Mean err\tMSWD Test err\tYoungest Zircon err\tMetropolis Estimate err\tMetropolis+MSWD err\tWeighted Mean err/sigma\tMSWD Test err/sigma\tYoungest Zircon err/sigma\tMetropolis Estimate err/sigma\tMetropolis+MSWD err/sigma\n");
	}

	// Get number of simulations per MPI task from command-line argument
	const uint32_t simspertask = (uint32_t)abs(atoi(argv[1]));
	// Get number of steps for metropolis walker to take
	const uint32_t nsteps = (uint32_t)abs(atoi(argv[2]));
	// Get dt/sigma value to run at
	const double dt_sigma = fabs(atof(argv[3]));


	// Import data
	const double* dist = csvparseflat(argv[4],'\t', &distrows, &distcolumns);

	double* uniformdist = malloc(distrows*sizeof(double));
	for (i=0;i<distrows;i++){
		uniformdist[i]=1;
	}

	// Declare various variables
	pcg32_random_t rng;
	pcg32_srandom_r(&rng,time(NULL), world_rank);

	const uint32_t Nmin = 1;
	const uint32_t Nmax = 1024;
	uint32_t Ns[] = {1,2,3,4,6,8,11,16,23,32,45,64,91,128,181,256,362,512,724,1024};
	uint32_t nNs = 20;

	double tmin_true = 100;
	double reluncert = 0.1/100.0;
	double relagerange = reluncert * dt_sigma;
	double absuncert = reluncert*tmin_true;
	double tmax_true = tmin_true*(1+relagerange);

	double* data = malloc(Nmax * sizeof(double));
	double* uncert = malloc(Nmax * sizeof(double));
	double wx, wsigma, mswd, minuncert; // For weighted mean
	double tmin_metropolis_uniform, tmin_metropolis_est, tmin_obs, tmin_mswd_test; // Means
	double tmin_metropolis_uniform_sigma, tmin_metropolis_est_sigma, tmin_obs_sigma, tmin_mswd_test_sigma; // Standard deviations

	for (i=0; i<nNs; i++){
		N = Ns[i];

		for (j=0; j<N; j++){
			uncert[j] = tmin_true * reluncert;
		}

		for (j=0; j<simspertask; j++){
			generateSyntheticZirconDataset(&rng, dist, distrows, tmin_true, tmax_true, uncert, data, N);

			// Find metropolis estimate of eruption age
			findMetropolisEstimate(&rng, dist, distrows, data, uncert, N, nsteps, &tmin_metropolis_est, &tmin_metropolis_est_sigma);
			tmin_metropolis_est_sigma = 1.253*tmin_metropolis_est_sigma; //Convert from mean absolute deviation to standard deviation
			//Check that metroplis uncertainty hasn't fallen below theoretical minimum
			minuncert = nanmean(uncert,N)/sqrt(N);
			if (tmin_metropolis_est_sigma < minuncert){
				tmin_metropolis_est_sigma = minuncert;
			}


			// Find metropolis estimate of eruption age with uniform (noninformative) prior distribution
			findMetropolisEstimate(&rng, uniformdist, distrows, data, uncert, N, nsteps, &tmin_metropolis_uniform, &tmin_metropolis_uniform_sigma);
			tmin_metropolis_uniform_sigma = 1.253*tmin_metropolis_uniform_sigma; //Convert from mean absolute deviation to standard deviation
			//Check that metroplis uncertainty hasn't fallen below theoretical minimum
			minuncert = nanmean(uncert,N)/sqrt(N);
			if (tmin_metropolis_uniform_sigma < minuncert){
				tmin_metropolis_uniform_sigma = minuncert;
			}

			// Find weighted mean of youngest n zircons with acceptable MSWD ("MSWD test")
			sort_doubles(data, N);
			tmin_mswd_test = data[0];
			tmin_mswd_test_sigma = absuncert;
			for (k=2; k<N+1; k++){
				awmean(data, uncert, k, &wx, &wsigma, &mswd);
				if (mswd > 1.0 + 2.0 * sqrt(2.0/(double)(k-1))){
					break;
				}
				tmin_mswd_test = wx;
				tmin_mswd_test_sigma = wsigma;
			}

			// Find weighted mean
			awmean(data, uncert, N, &wx, &wsigma, &mswd);

			// Find youngest zircon
			tmin_obs = minArray(data, N);
			tmin_obs_sigma = absuncert;

			printf("%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", N, mswd, (wx-tmin_true)/absuncert, (tmin_mswd_test-tmin_true)/absuncert, (tmin_obs-tmin_true)/absuncert, (tmin_metropolis_est-tmin_true)/absuncert, (tmin_metropolis_uniform-tmin_true)/absuncert, (wx-tmin_true)/wsigma*1.253, (tmin_mswd_test-tmin_true)/tmin_mswd_test_sigma*1.253, (tmin_obs-tmin_true)/tmin_obs_sigma*1.253, (tmin_metropolis_est-tmin_true)/tmin_metropolis_est_sigma*1.253, (tmin_metropolis_uniform-tmin_true)/tmin_metropolis_uniform_sigma)*1.253;
			//N.B.: factor of 1.253 converts from mean absolute deviation to standard deviation
		}
	}

MPI_Finalize();
return 0;
}
