# A suite for MCMC by Gibb's sampling
#  While tailored to fitting experimental cone flash responses, the 
#  routines are flexible and should be possible to adapt to other
#  applications.
using Random, LinearAlgebra, DataFrames, CSV

#%% gibbsdata
"""
A structure to specify settings for the gibbs sampler. The associated
MersenneTwister object which controls the random number generation is 
handled separately, since that object is of mutable type.
"""
struct gibbsdata
	# Number of free parameters
	dim::Int64
	
	# Logicals for which model params to vary
	prmvary::Dict{Symbol,Bool}

	# Ranges that can be used to absolutely bound model parameters
	#  These ranges simply bound the prior's support and needn't imply
	#  the prior is uniform. 
	#  Should be same keys as prmvary
	prmrg::Dict{Symbol,Array{Float64,1}}
	
	# Dictionary keys for parameters
	prmkeys::Array{Symbol,1}

end

#%% gibbsdatamat
"""
Encode the data specific to the model for cone flash response
nsmp: Number of mcmc samples desired
"""
function gibbsdatamat()
	# Number of model parameters
	dim = 178;

	#-----
	# Disease parameters
	#  Incubation period
	d_E = [2.,6.]; flagd_E = false; # days
	#  Infectious period
	d_I = [4.,7.]; flagd_I = false; # days

	prmrg=Dict{Symbol,Array{Float64,1}}(:d_E=>d_E,:d_I=>d_I);
	prmvary=Dict{Symbol,Bool}(:d_E=>flagd_E,:d_I=>flagd_I);

	#-----
	# Epidemic transmission
	#  r0
	r0 = [.5,2.]; flagr0 = true;
	prmrg[:r0] = r0; prmvary[:r0] = flagr0;

	#-----
	# Vaccination
	#  Susceptibility
	α = [0.,1.]; flagα = true;
	#  Contagiousness
	ω = [0.,1.]; flagω = true;
	#  Fraction of pop not willing to be vax'd
	vh = [0.,1.]; flagvh = false;
	prmrg[:α] = α; prmrg[:ω] = ω; prmrg[:vh] = vh;
	prmvary[:α] = flagα; prmvary[:ω] = flagω; prmvary[:vh] = flagvh;

	#-----
	# Time discretization
	#  Runga-kutta step
	δt = [.5e-1,1.]; flagδt = false;
	prmrg[:δt] = δt;
	prmvary[:δt] = flagδt;

	#-----
	# Auxilliary parameters
	rptλ = [1.,4.]; flagrptλ = true;
	bayσ = [25.,60.]; flagbayσ = true;
	vι0 = [0.,1.]; flagvι0 = true;
	prmrg[:rptλ] = rptλ; prmrg[:bayσ] = bayσ; prmrg[:vι0] = vι0;
	prmvary[:rptλ] = flagrptλ; prmvary[:bayσ] = flagbayσ; prmvary[:vι0] = flagvι0;

	#-----
	# Aggregate dictionary keys
	prmkeys = Array{Symbol,1}(undef,length(keys(prmvary)));
	pos = 1;

	for key in keys(prmvary)
		prmkeys[pos] = key;
		pos += 1;
	end

	#-----
	# Write data to gibbsdata
	return gibbsdata(dim,prmvary,prmrg,prmkeys)	
end

#%% gibbsmodelrun
"""
A wrapper for the gibbs routine to run the model for which parameters are
estimated
"""
function gibbsmodelrun(sheet::data,mydep::Dict{Symbol,Vector{Float64}},frc_M::Matrix{Float64})
	# Simulate the model
	tpts,ypts = odesolver(sheet,mydep,frc_M);

	# Return outputs
	return tpts,ypts

end

#%% gibbsmodelerr
"""
Compute the error between the model prediction and reported experimental data

measurements:: A dictionary whose key's give (time,experiment) arrays that record time series for data. 
	       time and experiment are in two col. Your sheet.tspan SHOULD NOT go past ODH csv
"""
function gibbsmodelerr(sheet::data,myaux::Dict{Symbol,Float64},mydep::Dict{Symbol,Vector{Float64}},
		       frc_M::Matrix{Float64},
		       measurements::Dict{String,Vector{Float64}})
	# Run the model for given choice of parameters
	tpts,yptsorig = gibbsmodelrun(sheet,mydep,frc_M)
	
	# For comparing to ODH, aggregate across Unvax, Vax, Unwill
	ypts = reshape(yptsorig,(27,5,length(tpts)));
	ypts = ypts[1:9,:,:] + ypts[10:18,:,:] + ypts[19:27,:,:];

	# Quadrature compute daily new infections
	#  Approx as ∫ₐᵇfdx ≈ ∑ᵢfᵢΔx = (b-a)*∑ᵢfᵢ/n
	ti = Int64(floor(sheet.tspan[1])); tf = Int64(ceil(sheet.tspan[2]));
	Etot = ypts[:,2,:];

	dailyI = Matrix{Float64}(undef,tf-ti,9);
	for i=ti:(tf-1)	
		for j=1:9
			val1 = myinterp(tpts,Etot[j,:],Float64(i));
			val2 = myinterp(tpts,Etot[j,:],Float64(i+1));
			# Integrate over day so b-a = 1 by trapezoidal
			dailyI[i-(ti-1),j] = 1/sheet.d_E*.5*(val1+val2);
		end

	end
	dailyI *= (1/myaux[:rptλ]);

	# Construct the error dictionary
	SE = Matrix{Float64}(undef,tf-ti,9);
	keys = ["0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"];
	for i=1:9
		SE[:,i] = abs.(dailyI[:,i] - measurements[keys[i]]).^2 ;
	end

	return SE
end
#%% randdir!
"""
Sample a uniform random direction in n-dimensional Euclidean space. It 
mutates the rng.
"""
function randdir!(n::Int64; rng::MersenneTwister=MersenneTwister())
	# Sample n ind N(0,1) random vectors and normalize to unit length
	dir = randn(rng,n);

	L = 0;
	for i=1:n
		L += dir[i]^2;
	end
	L = sqrt(L);

	return dir/L
end

#%% randOn
"""
Sample a uniformly random rotation in n-dimensional Euclidean space. It
mutates the rng.
"""
function randOn!(n::Int64; rng::MersenneTwister=MersenneTwister())
	# Use the QR factorization of a matrix built from iid standard normals 
	N = randn(rng,n,n);
	FacQR = qr(N);

	return FacQR.Q
end

#%% gibbsprior
"""
Evaluate log of the unnormalized Bayesian prior for given choice of model parameters
"""
function gibbsprior(sheet::data,myaux::Dict{Symbol,Float64},mydep::Dict{Symbol,Vector{Float64}},
		    gibbssheet::gibbsdata)
	# Check if sheet values are inside constraint ranges
	for key in gibbssheet.prmkeys
		if gibbssheet.prmvary[key]
			if !in(key,keys(myaux))
				myval = getfield(sheet,key);
			else
				myval = myaux[key];
			end
			if (myval < gibbssheet.prmrg[key][1])|(myval > gibbssheet.prmrg[key][2])
				return -Inf
			end
		end
	end

	#-----
	# Custom checks that would need to be changed for generic Gibbs
	# Check that we have not exceeded Ohio population
	pop = mydep[:Sv0] + 
	      mydep[:E0] + mydep[:Ev0] + mydep[:Eu0] + 
	      mydep[:I0] + mydep[:Iv0] + mydep[:Iu0] + 
	      mydep[:R0] + mydep[:Rv0] + mydep[:Ru0] + 
	      mydep[:D0] + mydep[:Dv0] + mydep[:Du0] ;
	flag = (pop .> sheet.N);
	nnz = sum(flag);
	if nnz != 0	
		return -Inf
	end

	# Passed all checks
	return 0.
	
end

#%% gibbsinit
"""
Randomly initialize gibbs mcmc by a value for which the prior is nonzero. It only mutates the rng.
"""
function gibbsinit!(sheet::data,myaux::Dict{Symbol,Float64},gibbssheet::gibbsdata; rng::MersenneTwister=MersenneTwister())
	flagfd = false;
	
	# Revert sheet to dictionary
	mydat = datamat(sheet)
	
	candsheet = 0;
	candaux = 0;
	canddepmat = 0;
	while !flagfd

		# Randomly sample freely varied parameter values
		for key in gibbssheet.prmkeys
			if gibbssheet.prmvary[key]
				smp = rand(rng);
				if !in(key,keys(myaux))
					mydat[key] = (1-smp)*gibbssheet.prmrg[key][1]+smp*gibbssheet.prmrg[key][2];
				else
					myaux[key] = (1-smp)*gibbssheet.prmrg[key][1]+smp*gibbssheet.prmrg[key][2];
				end
			end
		end

		# Update dependent parameters
		candsheet = data(mydat);
		canddepmat = depmat(candsheet,myaux);
		
		# Check if admissible
		flagfd = (gibbsprior(candsheet,myaux,canddepmat,gibbssheet) != -Inf)

	end

	return candsheet,canddepmat,myaux

end

#%% gibbslikelihood
"""
Evaluate the log of the unnormalized Bayesian likelihood for given choice of model parameters.
Routine uses myinterp from the odesolver.jl library
"""
function gibbslikelihood(sheet::data,mydep::Dict{Symbol,Vector{Float64}},myaux::Dict{Symbol,Float64},
		SE::Matrix{Float64})

	# Assume daily reported cases error is normal iid across each age cohort
	val = 0
	for i=1:length(SE)
		val += -.5*SE[i]/myaux[:bayσ]^2 - log(myaux[:bayσ]);
	end

	return val
end

#%% gibbscondprp
"""
Evaluate the log of Metropolis within Gibbs proposal distribution for given choice of model parameters
"""
function gibbscondprp(sheet::data,mydep::Dict{Symbol,Vector{Float64}},myaux::Dict{Symbol,Float64},
		      SE::Matrix{Float64})

	# Use a proposal distribution like 1/((max(SE))/sigma^2+2.) (= 1/(max(z-score)^2+.2.))
	#  Note: that because of how used in condsmp! that this is actually multiplied 
	#        by the prior
	mymax = maximum(SE);
	
	val = -log(mymax/myaux[:bayσ]^2+2.);

	return val
end

#%% gibbscondsmp!
"""
Sample the 1d conditional probabilities of the Bayesian posterior like Metropolis 
within Gibbs for MCMC. Code only mutates the rng.

prmkey:: String that names the parameter conditioned on all other frozen values
"""
function gibbscondsmp!(sheet::data,myaux::Dict{Symbol,Float64},gibbssheet::gibbsdata,
		       frc_M::Matrix{Float64},measurements::Dict{String,Vector{Float64}},SE::Matrix{Float64},
		       prmkey::Symbol; rng::MersenneTwister=MersenneTwister())


	# Sample by Metropolis-Hastings within Gibbs
	#  Proposal distribution sampled by accept-reject
	#   Take an upper envelope before taking log. This upper bound is consistent with gibbscondprp
	uenv = 0.5;	

	canddatmat = datamat(sheet);
	prm0 = !in(prmkey,keys(myaux)) ? 
	          canddatmat[prmkey] : myaux[prmkey];
	candauxmat = deepcopy(myaux);
	candSE = deepcopy(SE);
	candsheet = deepcopy(sheet); 
	mydep = depmat(sheet,myaux); canddepmat = deepcopy(mydep);	
	prmpr = 0;
	flagfd = false;
	nrej = 0; prgbar = 0; δprgbar = 1000;
	#  Run accept-reject on the proposal distribution
	while !flagfd
		#  Uniformly sample envelope
		smp = rand(rng); prmpr = (1-smp)*gibbssheet.prmrg[prmkey][1]+smp*gibbssheet.prmrg[prmkey][2];
		smp = rand(rng); prmgr = smp*uenv;
		
		# Update for new random sample and evaluate model error for these parameters if required
		if !in(prmkey,keys(myaux))
			canddatmat[prmkey] = prmpr;
			candsheet = data(canddatmat);
			canddepmat = depmat(candsheet,candauxmat);
			
			# Fast check if admissible
			logρ = gibbsprior(candsheet,candauxmat,canddepmat,gibbssheet);
			if logρ == -Inf	
				continue
			end

			candSE = gibbsmodelerr(candsheet,candauxmat,canddepmat,frc_M,measurements);
		else	
			candauxmat[prmkey] = prmpr;
			
			# Fast check if admissible
			logρ = gibbsprior(candsheet,candauxmat,canddepmat,gibbssheet);
			if logρ == -Inf
				continue
			end

		end
			

		#  Uniformly restrict to subgraph
		logρ = gibbscondprp(candsheet,canddepmat,candauxmat,candSE)
		if log(prmgr) < logρ
			flagfd = true;
		else
			nrej += 1;
			if nrej >= prgbar + δprgbar
				prgbar = nrej;
				println(String(prmkey)*" has $prgbar many rejections in its Gibbs sampling...")
			end
		end
	end
	
	# Compute the log of the Metropolis-Hastings acceptance ratio
	#  Assymetric proposal means ratio is pi(y)/pi(x)*q(x|y)/q(y|x)
	#   For our accept-reject scheme, x and y are independent under q
	mhlogratio = gibbslikelihood(candsheet,canddepmat,candauxmat,candSE) - 
	                 gibbslikelihood(sheet,mydep,myaux,SE) +
		      gibbscondprp(sheet,mydep,myaux,SE) - 
		          gibbscondprp(candsheet,canddepmat,candauxmat,candSE);
	
	# Perform Metropolis-Hastings accept-reject
	if log(rand(rng)) < mhlogratio
		return prmpr,canddepmat,candauxmat,candSE
	else
		return prm0,mydep,myaux,SE
	end

end

#%% mysaverng
"""
Given a MersenneTwister rng, save its internal states to csv files which may be reloaded for
restarting future runs.
"""
function mysaverng(rng::MersenneTwister)
	# Loop over fieldnames of rng and save to their own csv's
	#  One of the fields was a structure of its own and needed slightly different handling
	for key in fieldnames(typeof(rng))
		if key != :state
			df = DataFrame(string(key)=>getfield(rng,key));
		else
			ram = getfield(rng,key);
			df = DataFrame(string(key)=>ram.val);
		end
		CSV.write("RNG"*string(key)*".csv",df);
	end
	return true
end

#%% myloadrng
"""
Given a MersenneTwister saved into csv's like my mysaverng function, restore a rng with these settings.
fname is the prefix that the CSV files begin with.
"""
function myloadrng(fname::String="RNG")
	fields = ["seed","state","vals","ints","idxF","idxI","adv","adv_jump","adv_vals","adv_ints"];
	
	# Seed
	DF = CSV.read(fname*"seed.csv",DataFrame);
	myseed = convert(Vector{UInt32},DF[!,1]);

	# State
	DF = CSV.read(fname*"state.csv",DataFrame);
	mystate = Random.DSFMT.DSFMT_state(convert(Vector{Int32},DF[!,1]));

	# vals
	DF = CSV.read(fname*"vals.csv",DataFrame);
	myvals = convert(Vector{Float64},DF[!,1]);

	# ints
	DF = CSV.read(fname*"ints.csv",DataFrame);
	myints = convert(Vector{UInt128},DF[!,1]);

	# idxF
	DF = CSV.read(fname*"idxF.csv",DataFrame);
	myidxF = convert(Int64,DF[!,1][1]);

	# idxI
        DF = CSV.read(fname*"idxI.csv",DataFrame);
	myidxI = convert(Int64,DF[!,1][1]);

	# adv
	DF = CSV.read(fname*"adv.csv",DataFrame);
	myadv = convert(Int64,DF[!,1][1]);

	# adv_jump
	DF = CSV.read(fname*"adv_jump.csv",DataFrame);
	myadv_jump = convert(BigInt,DF[!,1][1]);

	# adv_vals
	DF = CSV.read(fname*"adv_vals.csv",DataFrame);
	myadv_vals = convert(Int64,DF[!,1][1]);

	# adv_ints
	DF = CSV.read(fname*"adv_ints.csv",DataFrame,type=Int64);
	myadv_ints = convert(Int64,DF[!,1][1]);

	return MersenneTwister(myseed,mystate,myvals,myints,myidxF,myidxI,myadv,myadv_jump,myadv_vals,myadv_ints)
end

#%% gibbsmcmc
"""
Wrapper for the entire Gibbs sampling procedure

MHmix:: is optional argument that combines Gibbs sampler with Metropolis-Hastings like a mixture
        model to ensure ergodicity even though the prior may conceivably introduce disconectedness 
	in the posterior. The proposal distribution is a simple uniform draw from ranges in 
	gibbssheet. 
"""
function gibbsmcmc(nsmp::Int64; rng::MersenneTwister=MersenneTwister(), MHmix::Float64=0., frestart::String="START")
	@assert (nsmp >= 2) "Need at least two samples"
	
	# Extract csv file names from standard sheet for loading and writing dictionaries to arrays
	sheet = data(); csv_vac = sheet.csv_vac; csv_odh = sheet.csv_odh;
	
	# Load the ODH files
	cols = ["0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"];
	DF = CSV.read(sheet.csv_odh,DataFrame);
	#  Store ODH daily infections over truncating to ti->tf interval
	#   ti,tf are days after Jan 1, 2020
	ti = Int64(floor(sheet.tspan[1])); tf = Int64(ceil(sheet.tspan[2]));
	measurements = Dict{String,Vector{Float64}}();
	for i=1:9
		measurements[cols[i]] = mymvavg(DF[!,cols[i]][ti:tf-1]);
	end	
	#  Load vaccination data
	frc_M = (isempty(sheet.csv_vac)) ? [0. 0.; 1. 0.] : vaxld(); # vaxld adjusts to ti=0-based index

	# Define the gibbs data sheet
	gibbssheet = gibbsdatamat();

	# Initialize a starting value
	if frestart == "START"	
		sheet = data(); myaux = auxmat();

		initsheet, initdepmat, initauxmat = gibbsinit!(sheet,myaux,gibbssheet; rng=rng)
		initdatamat = datamat(initsheet);
		initSE = gibbsmodelerr(initsheet,initauxmat,initdepmat,frc_M,measurements);
	else 
		# Load the data from files
		DF_init = CSV.read(frestart,DataFrame,header=false);
		M_init = convert(Matrix{Float64},DF_init[:,:]);
		rng = myloadrng();
		
		# Restore the last run
		initdatamat,initauxmat = csvdat(M_init[:,end],csv_vac,csv_odh);	
		initsheet = data(initdatamat);
		initdepmat = depmat(initsheet,initauxmat);
		initSE = gibbsmodelerr(initsheet,initauxmat,initdepmat,frc_M,measurements);

		# Clear memory hogs
		DF_init = 0;
		M_init = 0;
	end	
	println("Found MCMC starting value...");

	# Run the gibbs sampler
	#  Initialize and store first sample
	initprm,pos = csvdat(initdatamat,initauxmat);
	nprm = length(initprm);
	Posterior = Array{Float64,2}(undef,nprm,nsmp+1); 
	Posterior[:,1] = initprm;

	Errors = Array{Float64,2}(undef,nsmp+1,1);
	Errors[1,1] = sum(initSE);
	
	#  For writing out progress
	prgbar = 0; δprgbar = .01; nMH = 0; nGibbs = 0;

	@inbounds for i=2:(nsmp+1)
		if rand(rng) < MHmix
			# MH sample by accept-reject under gibbscondprp
			candprm = 0; candsheet=0; canddepmat=0; candauxmat=0; candSE=0;
			flag_fd = false; nrej = 0.; mhprgbar = 0.; mhδprgbar = 1000.;
			while !flag_fd
				# Sample cand parameter values
				candprm = deepcopy(initdatamat); 
				candauxmat = deepcopy(initauxmat);
				for key in gibbssheet.prmkeys
					if gibbssheet.prmvary[key]
						smp = rand(rng); smp = (1-smp)*gibbssheet.prmrg[key][1]+smp*gibbssheet.prmrg[key][2]
						if in(key, keys(candprm))
							candprm[key] = smp;
						else
							candauxmat[key] = smp;
						end
					end
				end
		
				candsheet = data(candprm);
				canddepmat = depmat(candsheet,candauxmat);

				logρ = gibbsprior(candsheet,candauxmat,canddepmat,gibbssheet);	
				if (logρ == -Inf)
					continue
				end

				#  Take an upper envelope before taking log. This upper bound is consistent with gibbscondprp
				uenv = 0.5;
				prmgr = rand(rng)*uenv;

				candSE = gibbsmodelerr(candsheet,candauxmat,canddepmat,frc_M,measurements);
				logρ = gibbscondprp(candsheet,canddepmat,candauxmat,candSE);
				
				#   Uniformly restrict to subgraph 
				if log(prmgr) < logρ
					flag_fd = true;
				else
					nrej += 1
					if nrej >= mhprgbar + mhδprgbar
						mhprgbar = nrej;
						println("Pure MH has $mhprgbar many rejections in its accept-reject...")
					end
				end


			end
			
			# Compute log of MH acceptance factor
			#  Assymetric proposal means ratio is pi(y)/pi(x)*q(x|y)/q(y|x)
			#   For our accept-reject scheme, x and y are independent under q
			mhlogratio = gibbslikelihood(candsheet,canddepmat,candauxmat,candSE) -
			                 gibbslikelihood(initsheet,initdepmat,initauxmat,initSE) +
				     gibbscondprp(initsheet,initdepmat,initauxmat,initSE) -
				         gibbscondprp(candsheet,canddepmat,candauxmat,candSE);	    
						
			if log(rand(rng)) < mhlogratio
 				# Accept the proposed candidate
				initdatamat = candprm;
				initsheet = candsheet;
				initdepmat = canddepmat; 
				initauxmat = candauxmat;
				initSE = candSE;
			end
			nMH += 1;
		else
			# Gibbs sample
			for key in gibbssheet.prmkeys
				if gibbssheet.prmvary[key]
					prmpr,canddepmat,candauxmat,candSE = gibbscondsmp!(initsheet,initauxmat,gibbssheet,
						                                            frc_M,measurements,initSE,
						                                            key; rng=rng);
					if !in(key,keys(initauxmat))
						initdatamat[key] = prmpr;
					else
						initauxmat[key] = prmpr;
					end
	
					initsheet = data(initdatamat);
					initdepmat = canddepmat;
					initauxmat = candauxmat;
					initSE = candSE;
				
				end
			end	
			nGibbs += 1;
		end
		
		Posterior[:,i],pos = csvdat(initdatamat,initauxmat);
		Errors[i,1] = sum(initSE)/length(initSE);

		# Print progress
		myprg = i/(nsmp+1);
		if (myprg >= prgbar + δprgbar)
			prgbar = floor(myprg/δprgbar)*δprgbar;
			CSV.write("GibbsMCMC.csv", DataFrame(Posterior[:,1:i]), writeheader=false);
			CSV.write("ErrorsMCMC.csv", DataFrame(Errors[1:i,:]), writeheader=false);
			mysaverng(rng);
			println("$prgbar" *"/1 complete with MCMC samples ...")
		end
	end

# Save chain parameter values to csv
println("Number of MH samples: $nMH");
println("Number of Gibbs samples: $nGibbs");
CSV.write("GibbsMCMC.csv", DataFrame(Posterior), writeheader=false);
CSV.write("ErrorsMCMC.csv", DataFrame(Errors), writeheader=false);
mysaverng(rng);

end

#%% gibbscorr
"""
Given two sample vectors of 1d random variables, compute their correlation
"""
function gibbscorr(X::Array{Float64,1},Y::Array{Float64,1})
	@assert length(X) == length(Y) "Need a common number of samples"
	nmcmc = length(X);

	# Compute the means
	Xμ = sum(X)/nmcmc;
	Yμ = sum(Y)/nmcmc;

	# Compute the std dev
	Xvar = X .- Xμ;
	Xvar = Xvar.^2;
	Xσ = sum(Xvar)/nmcmc;
	Xσ = sqrt(Xσ);
	
	Yvar = Y .- Yμ;
	Yvar = Yvar.^2;
	Yσ = sum(Yvar)/nmcmc;
	Yσ = sqrt(Yσ);

	# Compute the correlation
	cor = sum( (X .- Xμ).*(Y .- Yμ) )/nmcmc/( Xσ*Yσ );

	return cor, Xμ, Yμ, Xσ, Yσ

end

#%% mergeoutputs
"""
Given consecutive gibbs sample runs collected in their own Run# folder,
merge these into a common arrays for use in the gibbs_postprocess.pynb
"""
function mergeoutputs(idbeg::Integer,idend::Integer;flag_write::Bool=false)
	M_mcmc = 0;
	M_err = 0;

	for k=idbeg:idend
		# Load a csv
		df_mcmc = CSV.read("Run"*string(k)*"/GibbsMCMC.csv",DataFrame,header=false);
		df_err = CSV.read("Run"*string(k)*"/ErrorsMCMC.csv",DataFrame,header=false);
		
		# Merge into the matrices
		if k == idbeg
			M_mcmc = convert(Matrix,df_mcmc[:,2:end]);
			M_err = convert(Vector,df_err[!,1][2:end]);
		else
			M_mcmc = [M_mcmc convert(Matrix,df_mcmc[:,2:end])];
			M_err  = [M_err; convert(Vector,df_err[!,1][2:end])];
		end

		println("Finished Run"*string(k)*"...")
	end

	if flag_write
		# Write csv's
		CSV.write("GibbsMCMC.csv",DataFrame(M_mcmc),writeheader=false);
		CSV.write("ErrorsMCMC.csv",DataFrame(reshape(M_err,(length(M_err),1))),writeheader=false);
	end

	return M_mcmc, M_err

end

#%% gibbstraj
"""
Compute trajectories for all the parameters output by gibbsmcmc
"""
function gibbstraj(idbeg::Int64,idend::Int64;
		   dwnsmp::Int64=50,flag_write::Bool=true,tspan=Vector{Float64}([NaN,NaN]))
	
	sheet = data();
	frc_M = (isempty(sheet.csv_vac)) ? [0. 0.; 1. 0.] : vaxld();
	# Load tspan from data if it was not specified
	if tspan == [NaN,NaN]	
		tspan = sheet.tspan;
	end

	# Load the MCMC smps
	#  Used to group trajs into suff small chunks so can be saved. Stores 
	#  number of samples in each Run.
	nmcmc = Vector{Int64}(undef,idend-(idbeg-1));
	shift = -(idbeg-1);
	for i=idbeg:idend
		df_mcmc = CSV.read("Run"*string(i)*"/GibbsMCMC.csv",DataFrame,header=false);
		nmcmc[i+shift] = size(df_mcmc)[2]-1; # Bc mergeoutputs discards overlapping first sample
	end

	M_mcmc,M_err = mergeoutputs(idbeg,idend);

	# Run the trajectories
	nsmp = size(M_mcmc,2);
	#  Define the axis
	taxis = [tval for tval in LinRange(tspan[1],tspan[2],dwnsmp+1)];
	#  Evaluate the models
	M_trajs = Matrix{Float64}(undef,9*dwnsmp,nsmp+1);
	M_trajs[:,1] = repeat(taxis[1:end-1],outer=(9,));
	pos = 0.; δprg = .005;
	for j=1:nsmp
		mydat,myaux = csvdat(M_mcmc[:,j],sheet.csv_vac,sheet.csv_odh);
		mydat[:tspan] = tspan;
		mysheet = data(mydat);
		mydep = depmat(mysheet,myaux);
		tpts,yptsorig = gibbsmodelrun(mysheet,mydep,frc_M);

		#-----
		# For comparing to ODH, aggregate across Unvax, Vax, Unwill
		ypts = reshape(yptsorig,(27,5,length(tpts)));
		ypts = ypts[1:9,:,:] + ypts[10:18,:,:] + ypts[19:27,:,:];
		# Quadrature compute daily new infections
		#  Approx as ∫ₐᵇfdx ≈ ∑ᵢfᵢΔx = (b-a)*∑ᵢfᵢ/n
		Etot = ypts[:,2,:]; # Age x time
		dailyI = Matrix{Float64}(undef,dwnsmp,9); # time x Age
		for i=1:dwnsmp
			for k=1:9
				val1 = myinterp(tpts,Etot[k,:],taxis[i]);
				val2 = myinterp(tpts,Etot[k,:],taxis[i+1]);
				# Integrate over day so b-a = 1 by trapezoidal
				dailyI[i,k] = 1/mydat[:d_E]*.5*(val1+val2);
			end	
		end
		M_trajs[:,j+1] = 1/myaux[:rptλ]*dailyI[:];
		#-----
		
		prg = j/nsmp;
		if prg >= pos + δprg
			inc = floor((prg-pos)/δprg);
			pos += inc*δprg;
			println("Fraction $pos through with samples ...");
		end

	end

	# Save to csv's
	if flag_write
		shift = -(idbeg-1);
		for i=idbeg:idend
			# Identify block of M_trajs with this Run's trials
			pos1 = (i != idbeg) ? 1+sum(nmcmc[1:i+shift-1]) : 2 ;
			pos2 = 1+sum(nmcmc[1:i+shift]);
		
			# Concatenate with first column being time axis
			M_ram = [M_trajs[:,1] M_trajs[:,pos1:pos2]];
			
			# Write to csv
			CSV.write("Run"*string(i)*"/GibbsTrajs.csv",DataFrame(M_ram),writeheader=false);
		end
	end

	return M_trajs
end
