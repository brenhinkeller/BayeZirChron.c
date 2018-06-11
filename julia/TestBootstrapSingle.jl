## --- Load packages

    using Plots; gr();
    using KernelDensity: kde

    # Functions we'll be using here
    include("DistTools.jl")

## ---

    dt_sigma = 2;	# Timescale relative to analytical uncertainty
    nsteps = 500000;	# Length of Markov chain
    burnin = 10000;	# Number of steps to discard at beginning of Markov chain

    N = 10
    # Draw new set of ages where true minimum == 0 and analytical sigma == 1
    ages = drawFromDistribution(MeltsVolcanicZirconDistribution,N).*dt_sigma + randn(N);
    uncert = ones(N)

    # Maximum extent of expected analytical tail (beyond eruption/deposition)
    maxTailLength = mean(uncert) * NormQuantile(1-1/(length(ages) + 1));
    included = (ages-minimum(ages)) .>= maxTailLength;

    # Bootstrapped crystallization distribution, excluding maximum analytical tail
    if sum(included) > 5
        dist = BootstrapCrystDistributionKDE(ages[included]);
    else
        # Avoid edge cases at n = 0 and n = 2;
        # Default to n = 1 instead, which yields a half-normal distribution
        dist = BootstrapCrystDistributionKDE([0.]);
    end
    plot(dist, label="bootstrapped", ylabel="density", legend=:topleft)
    plot!(linspace(0,80,100),MeltsVolcanicZirconDistribution,label="original")

## ---

    # Run MCMC to estimate saturation and eruption/deposition age distributions
    # (tminDist,~,~,~) = crystMinMaxMetropolis(nsteps,dist,ages,uncert);
    tminDist = Array{Float64,1}(nsteps);
    @time crystMinMetropolis(nsteps,dist,ages,uncert,tminDist);

    AgeEst = mean(tminDist[burnin:end]);
    AgeEst_sigma = std(tminDist[burnin:end]);
    print("$AgeEst +/- $AgeEst_sigma Ma")
