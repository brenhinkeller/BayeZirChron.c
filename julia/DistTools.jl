    # How far away from the mean (in units of sigma) should we expect proportion
    # F of the samples to fall
    function NormQuantile(F)
        return sqrt(2)*erfinv(2*F-1)
    end

    # Bootstrap a KDE of the pre-eruptive (or pre-deposition) zircon distribution
    # shape from a 2-d array of sample ages using a KDE of stacked sample data
    function BootstrapCrystDistributionKDE(data::Array{Float64})
        # Load all data points and scale from 0 to 1
        allscaled = Array{Float64,1}();
        for i=1:size(data,2)
            scaled = data[:,i]-minimum(data[:,i]);
            max = maximum(scaled);
            if max > 0
                scaled = scaled./max;
            end
            allscaled = [allscaled; scaled]
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7);
        t = kd.x .> -0.05;
        return kd.density[t];
    end

    # Draw random numbers from a distribution specified by a vector of points
    # defining the PDF curve
    function drawFromDistribution(dist::Array{Float64}, n::Int)
        # Draw n random numbers from the distribution 'dist'
        x = Array{Float64}(n);
        dist_ymax = maximum(dist);
        dist_xmax = length(dist)-1.0;

        for i=1:n;
            while true
                # Pick random x value
                rx = rand() * dist_xmax;
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand() * dist_ymax;
                if (y > ry)
                    x[i] = rx / dist_xmax;
                    break
                end
            end
        end
        return x
    end

    # Calculate a weigted mean, including MSWD, but without MSWD correction to uncertainty
    function awmean(x, sigma)
        n = length(x);
        s1 = 0.0; s2 = 0.0; s3 = 0.0;

        if n==1
            wx = x[1];
            mswd = 0;
            wsigma = sigma[1];
        else
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i]);
                s2 += 1 / (sigma[i]*sigma[i]);
            end
            wx = s1/s2;

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i]);
            end
            mswd = s3 / (n-1);
            wsigma = sqrt(1.0/s2);
        end
        return wx, wsigma, mswd
    end

    # Return the log-likelihood of a proposed crystallization distribution, with adjustments to prevent runaway at low N
    function checkCrystLogLikelihood(dist::Array{Float64}, data::Array{Float64}, uncert::Array{Float64}, tmin::Float64, tmax::Float64)
        # Define some frequently used variables
        loglikelihood = 0.0;
        datarows = length(data);
        distrows = length(dist);
        dist_xscale = distrows-1.00000000001;
        dist_yave = mean(dist);
        dt = abs(tmax-tmin);
        # Cycle through each datum in data array
        for j=1:datarows
            # If possible, prevent aliasing problems by interpolation
            if (uncert[j] < dt/dist_xscale) && data[j] > tmin && data[j] < tmax
                # Find (double) index
                ix = (data[j] - tmin) / dt * dist_xscale + 1;
                # Interpolate corresponding distribution value
                f = floor(Int,ix)
                likelihood = (dist[f+1]*(ix-f) + dist[f]*(1-(ix-f))) / (dt * dist_yave);
                # Otherwise, sum contributions from Gaussians at each point in distribution
            else
                likelihood = 0;
                for i=1:distrows
                    distx = tmin + dt*(i-1)/dist_xscale; # time-position of distribution point
                    # Likelihood curve follows a Gaussian PDF. Note: dt cancels
                    likelihood += dist[i] / (dist_yave * distrows * uncert[j] * sqrt(2*pi)) *
                    exp( - (distx-data[j])*(distx-data[j]) / (2*uncert[j]*uncert[j]) );
                end
            end
            loglikelihood += log(likelihood);
        end
        # Calculate a weighted mean and examine our MSWD
        wm, wsigma, mswd = awmean(data, uncert);
        if datarows == 1 || mswd < 1
            Zf = 1;
        elseif mswd*sqrt(datarows) > 1000
            Zf = 0;
        else
            f = datarows-1.0;
            # Height of MSWD distribution relative to height at MSWD = 1 (see Wendt and Carl, 1991, Chemical geology)
            Zf = exp((f/2-1)*log(mswd) - f/2*(mswd-1));
        end
        # To prevent instability / runaway of the MCMC for small datasets (low N),
        # favor the weighted mean interpretation at high Zf (MSWD close to 1) and
        # the youngest-zircon interpretation at low Zf (MSWD far from one). The
        # penalty factors used here are determined by training against synthetic datasets.
        return loglikelihood - (2/log(1+datarows)) * (                  # Scaling factor that decreases with log number of data points (i.e., no penalty at high N)
        log((abs(tmin - wm)+wsigma)/wsigma)*Zf +                        # Penalty for proposing tmin too far from the weighted mean at low MSWD (High Zf)
        log((abs(tmax - wm)+wsigma)/wsigma)*Zf +                        # Penalty for proposing tmax too far from the weighted mean at low MSWD (High Zf)
        log((abs(tmin - data[1])+uncert[1])/uncert[1])*(1-Zf) +         # Penalty for proposing tmin too far from youngest zircon at high MSWD (low Zf)
        log((abs(tmax - data[end])+uncert[end])/uncert[end])*(1-Zf) );  # Penalty for proposing tmax too far from oldest zircon at high MSWD (low Zf)
    end


    # Run a Metropolis sampler to estimate the extrema of a finite-range distribution from samples drawn
    # from that distribution -- e.g., estimate zircon saturation and eruption ages from a distribution of
    # zircon crystallization ages.
    function crystMinMetropolis(nsteps::Int,dist::Array{Float64},data::Array{Float64},uncert::Array{Float64},tminDist::Array{Float64,1})
        # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
        stepfactor = 2.9;
        # Sort the data array from youngest to oldest
        sI = sortperm(data);
        data = data[sI]; # Sort data
        uncert = uncert[sI]; # Sort uncertainty
        # These quantities will be used more than once
        datarows = length(data);
        tmin_obs = minimum(data);
        tmax_obs = maximum(data);
        # Step sigma for Gaussian proposal distributions
        dt = tmax_obs - tmin_obs + uncert[1] + uncert[end];
        tmin_step = dt / datarows;
        tmax_step = dt / datarows;
        # Use oldest and youngest zircons for initial proposal
        tmin = tmin_obs - uncert[1];
        tmax = tmax_obs + uncert[end];
        tmin_proposed = tmin;
        tmax_proposed = tmax;
        # Log likelihood of initial proposal
        ll =  checkCrystLogLikelihood(dist, data, uncert, tmin, tmax);
        ll_proposed = ll;
        # Step through each of the N steps in the Markov chain
        for i=1:nsteps
            tmin_proposed = copy(tmin);
            tmax_proposed = copy(tmax);
            # Adjust either upper or lower bound
            if rand()<0.5
                tmin_proposed += tmin_step*randn();
            else
                tmax_proposed += tmax_step*randn();
            end
            # Flip bounds if reversed
            if (tmin_proposed>tmax_proposed)
                r = tmin_proposed;
                tmin_proposed = tmax_proposed;
                tmax_proposed = r;
            end
            # Calculate log likelihood for new proposal
            ll_proposed =  checkCrystLogLikelihood(dist, data, uncert, tmin_proposed, tmax_proposed);
            # Decide to accept or reject the proposal
            if rand() < exp(ll_proposed-ll)
                if tmin_proposed != tmin
                    tmin_step = abs(tmin_proposed-tmin)*stepfactor;
                end
                if tmax_proposed != tmax
                    tmax_step = abs(tmax_proposed-tmax)*stepfactor;
                end

                tmin = copy(tmin_proposed);
                tmax = copy(tmax_proposed);
                ll = copy(ll_proposed);
            end
            tminDist[i] = tmin;
        end
        return 0
    end


## --- Function for parallel workers to evaluage

    function BootstrapEachN(MeltsVolcanicZirconDistribution,Ns,nsteps,burnin,dt_sigma)
        tminDist = Array{Float64,1}(nsteps);
        AgeEst = Array{Float64}(length(Ns));
        AgeEst_sigma = Array{Float64}(length(Ns));
        for j=1:length(Ns)
            N = Ns[j]
            # Draw new set of ages where true minimum == 0 and analytical sigma == 1
            ages = drawFromDistribution(MeltsVolcanicZirconDistribution,N).*dt_sigma + randn(N);
            uncert = ones(N)

            # Maximum extent of expected analytical tail (beyond eruption/deposition)
            maxTailLength = mean(uncert) * NormQuantile(1-1/(1+length(ages)));
            included = (ages-minimum(ages)) .>= maxTailLength;

            # Bootstrapped crystallization distribution, excluding maximum analytical tail
            if sum(included) > 5
                dist = BootstrapCrystDistributionKDE(ages[included]);
            else
                # Avoid edge cases at n = 0 and n = 2;
                # Default to n = 1 instead, which yields a half-normal distribution
                dist = BootstrapCrystDistributionKDE([0.]);
            end

            # Run MCMC to estimate saturation and eruption/deposition age distributions
            # (tminDist,~,~,~) = crystMinMaxMetropolis(nsteps,dist,ages,uncert);
            crystMinMetropolis(nsteps,dist,ages,uncert,tminDist);

            AgeEst[j] = mean(tminDist[burnin:end]);
    	    AgeEst_sigma[j] = std(tminDist[burnin:end]);

            if AgeEst_sigma[j] < 1/sqrt(N)
                # Ensure metroplis sigma is not below maximum theoretical precision
		# This can happen if markov chain gets stuck -- not usually a problem
		# and would be discovered by manual inspection of stationary distribution
		# but can lead to outliers when running many simulations in parallel
                AgeEst_sigma[j] = 1/sqrt(N);
	    end

            # print(i,", ", N, ": ", AgeEst_dist[k], " +/- ", AgeEst_sigma_dist[k], "\n")
        end
        return (AgeEst, AgeEst_sigma)
    end


## --- Distributions

MeltsVolcanicZirconDistribution = reverse([0.000126311537071135,0.138678719074772,0.277231126612473,0.415780749218191,0.554244033653465,0.692478823982408,0.828041603755384,0.958465872031531,1.08363596175037,1.20336837717604,1.31086441368344,1.39918689258733,1.46297882113932,1.50876500985997,1.53936076666456,1.55766152364649,1.56721340549175,1.57055368507083,1.57017608482977,1.56375900729345,1.55673917794689,1.54528204105819,1.53446510415179,1.51888561033322,1.50413627708598,1.48878532213817,1.47093626533172,1.45367531631901,1.43722401954538,1.41884953441590,1.40020284908090,1.38217173159262,1.36491496000665,1.34781339259003,1.32918673600548,1.31064062341720,1.29264559362098,1.27506639070296,1.25789091287708,1.24124868413795,1.22415364593791,1.20619847793269,1.18871940548722,1.17166291592868,1.15507107892548,1.13905814837756,1.12329852032616,1.10796490309432,1.09296239013487,1.07797063027672,1.06271754003706,1.04720311991123,1.03171712350443,1.01658666090717,1.00179914992514,0.987288063694791,0.973203816081034,0.959396520209570,0.945904001861541,0.932723654268405,0.919949056056459,0.907472933652423,0.895233778569715,0.883320760619869,0.871740160686158,0.860440380738532,0.849189937637448,0.836826698333265,0.824476857182168,0.812345940043929,0.800432420550727,0.788729129991484,0.777279410925614,0.766113420680074,0.755111910752756,0.744481086182191,0.734045338222697,0.723854066251978,0.713944076764060,0.704174797556201,0.694697443892042,0.685378323765953,0.676278206907816,0.667385068616828,0.658707640374667,0.650208560617337,0.641977872085481,0.633875778399037,0.625836373705725,0.617898448430796,0.610049502743681,0.602380243936629,0.594892672217814,0.587607589519472,0.580506613961880,0.573451694478375,0.566417381408414,0.559389964999820,0.552362549285121,0.545335133570427]);
