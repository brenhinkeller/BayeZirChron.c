## --- Load packages

    using Plots; gr();

## --- Call parallel workers

    # # On single system
    # np = 3;
    # addprocs(np)

    # On cluster
    np = 320;
    using ClusterManagers
    addprocs_slurm(np,time="00:30:00");

## --- Generate a synthetic dataset

    @everywhere begin
        # For bootstrapping
        using KernelDensity: kde

        # Functions we'll be using here
        include("DistTools.jl")

        # Run the test for the following number of analyses
        Ns = [1,2,3,4,6,8,11,16,23,32,45,64,91,128,181,256,362,512,724,1024];

        dt_sigma = 1;	# Timescale relative to analytical uncertainty
        nsteps = 10000;	# Length of Markov chain
        burnin = 1000;	# Number of steps to discard at beginning of Markov chain

    end

    NSims = np * 4; # Run eight simulations per processor

    const jobs = RemoteChannel(()->Channel{Int}(NSims));
    const results = RemoteChannel(()->Channel{Tuple}(NSims));

    @everywhere function bootstrap_task(jobs, results) # define work function everywhere
        while true
            job_id = take!(jobs)
            out = BootstrapEachN(MeltsVolcanicZirconDistribution,Ns,nsteps,burnin,dt_sigma);
            put!(results, out)
        end
    end

    function make_jobs(n)
        for i in 1:n
            put!(jobs, i)
        end
    end

    @schedule make_jobs(NSims); # feed the jobs channel

    for p in workers() # start tasks on the workers to process requests in parallel
        @async remote_do(bootstrap_task, p, jobs, results)
    end

    AgeEst_dist = Array{Float64}(NSims*length(Ns));
    AgeEst_sigma_dist = Array{Float64}(NSims*length(Ns));
    N_dist = Array{Float64}(NSims*length(Ns));

    n = NSims;
    @elapsed while n > 0 # print out results
        AgeEst, AgeEst_sigma = take!(results)
        thisrange = (n-1)*length(Ns)+(1:length(Ns));
        AgeEst_dist[thisrange] = AgeEst;
        AgeEst_sigma_dist[thisrange] = AgeEst_sigma;
        N_dist[thisrange] = Ns;
        n = n - 1;
    end

    A = Array{Any,2}(NSims*length(Ns)+1,3)
    A[1,:] = ["N","AgeEst","AgeEst_sigma"];
    A[2:end,:] =  hcat(N_dist,AgeEst_dist,AgeEst_sigma_dist);
    writedlm(string("BootstrappedEruptionAge_dt_sigma_", dt_sigma, ".csv"), A, ',')

## ---
    MAD = Array{Float64}(size(Ns));
    MAD_ExpectedMAD = Array{Float64}(size(Ns));

    for i=1:length(Ns)
        t = N_dist .== Ns[i];
        MAD[i] = mean(abs.(AgeEst_dist[t])); # Note: expected age = 0
        MAD_ExpectedMAD[i] = mean(abs.(AgeEst_dist[t])./(1.253 * AgeEst_sigma_dist[t]));
    end

    f1 = plot(Ns,MAD,xlims=[1,1000],xscale=:log10)
    ylims!(0,3)
    savefig(f1,string("AbsoluteError_dt_sigma", dt_sigma, ".pdf"))

    f2 = plot(Ns,MAD_ExpectedMAD,xlims=[1,1000],xscale=:log10)
    ylims!(0,6)
    savefig(f2,string("Absolute_ExpectedError_dt_sigma", dt_sigma, ".pdf"))


## ---
