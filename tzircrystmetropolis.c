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
 	// Draw N = xrows random numbers from the distribution 'dist'
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

 double checkZirconLikelihood(const double* restrict dist, const uint32_t distrows, const double* restrict data, const double* restrict uncert, const uint32_t datarows, const double tmin, const double tmax){
 	// Return the log-likelihood of drawing the dataset 'data' from the distribution 'dist'
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

 	// Calculate a weighted mean and examine our MSWD
 	awmean(data, uncert, datarows, &wm, &wsigma, &mswd);
 	if (datarows == 1 || mswd<1){
 		Zf = 1;
 	} else if (mswd*sqrt(datarows)>1000){
 		Zf = 0;
 	} else {
 		const double f = (double)datarows-1;
 		Zf = exp((f/2-1)*log(mswd) - f/2*(mswd-1)); // Height of MSWD distribution relative to height at mswd = 1;
 	}

 	// At low N (low sample numbers), favor the weighted mean interpretation at high Zf
 	// (MSWD close to 1) and the youngest-zircon interpretation at low Zf (MSWD far from one).
 	// This is helpful to for preventing instability at low N, with the funcional form
 	// used here derived from training against synthetic datasets.
 	return loglikelihood  - ( log10((fabs(tmin - wm)+wsigma)/wsigma)*Zf + log10((fabs(tmax - wm)+wsigma)/wsigma)*Zf + log10((fabs(tmin - data[0])+uncert[0])/uncert[0])*(1-Zf) + log10((fabs(tmax - data[datarows-1])+uncert[datarows-1])/uncert[datarows-1])*(1-Zf) ) * (2/log10(1+datarows));
 }


 int main(int argc, char **argv){
 	// Run a Metropolis sampler to quantify zircon saturation and eruption ages
 	// Invoke from command line as: tzircrystmetroplis1sigma  <nsteps> <distribution.tsv> <zircondata.tsv>
 	// Where 'nsteps; is the number of steps in the Markov chain, 'distribution.tsv' is a column vector containing
 	// a prior probability distribution for zircon crystallisation, and 'zircondata.tsv' is the raw zircon
 	// age spectrum and 1-sigma uncertainty in comma- or tab-separated columns

 	uint32_t distrows, distcolumns;
 	uint32_t datarows, datacolumns;
 	uint32_t i, j, k;

 	//Check input arguments
 	if (argc != 4) {
 		fprintf(stderr,"USAGE: %s <nsteps> <distribution.tsv> <zircondata.tsv>\n", argv[0]);
 		exit(1);
 	}


 	// Get number of steps for metropolis walker to take
 	const uint32_t nsteps = (uint32_t)abs(atoi(argv[1]));
 	// sigma = stepfactor * last-step
 	const double stepfactor = 2.9; // standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%


 	// Import data -- try both comma and tab-delimited files
 	double* dist = csvparseflat(argv[2],'\t', &distrows, &distcolumns); // distribution has only one column, delimiter doesn't matter
 	double* data = csvparseflat(argv[3],'\t', &datarows, &datacolumns);
 	if (datacolumns<2){
 		free(data);
 		data = csvparseflat(argv[3],',', &datarows, &datacolumns);
 	}

 	for (i=0; i<datarows; i++){
 		data[datarows*1+i]=data[datarows*1+i]/2; // Convert from 2-sigma to 1-sigma uncertainties
 	}
 	const double tmin_obs = minArray(data, datarows);
 	const double tmax_obs = maxArray(data, datarows);
 	const double dt = tmax_obs - tmin_obs + data[datarows*1+0] + data[datarows*1 + datarows-1];

	// Print statements for debugging, if necessary
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

 	theta =  checkZirconLikelihood(dist, distrows, data, &data[datarows*1], datarows, tmin_proposed, tmax_proposed);

 	double p = 10;

 	int transitions=0;

 	// Step through each of the N steps in the Markov chain
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
