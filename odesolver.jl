#%% myfindfirst
"""
A binary search routine for finding the first time point greater than the query point
among a given sequence of times. Writing because the findfirst routine is proving to
be expensive in Julia. Routine uses that tpts is ordered least to great. It's 
assumed teval is strictly contained in the interval partitioned by tpts.
"""
function myfindfirst(tpts::Array{Float64,1},teval::Float64)
	ntpts = length(tpts);
	
	# Find the smallest interval of type (,] containing point.
	idx = [1,ntpts];
	flag_fd = false;

	while !flag_fd
		mid = ceil(Int64,.5*idx[1]+.5*idx[2]);
		if mid == idx[2]
			flag_fd = true;
		elseif teval <= tpts[mid]
			idx[2] = mid;
		else
			idx[1] = mid;
		end
	end

	return idx[2]
end

#%% myinterp
"""
A simple 1d linear interpolation scheme to extend a discrete data set to an interval

vals: ntpts x 2 array of floats. First column is time, second is function value
teval: time point at which to evaluate
"""
function myinterp(tpts::Array{Float64,1},ypts::Array{Float64,1},teval::Float64)
	
	if teval <= tpts[1]
		val = ypts[1];
	elseif teval >= tpts[end]
		val = ypts[end];
	else
		pos = myfindfirst(tpts,teval);
		t1,t2 = tpts[pos-1:pos];
		s = (teval-t1)/(t2-t1);
		v1 = ypts[pos-1];
		v2 = ypts[pos];
		val = v1+s*(v2-v1);
	end

	return val

end

#%% vacprt
"""
Define the transfer rate from S unvacc to vacc
"""
function vacprt(sheet::data,t::Float64,S::Vector{Float64},
		Nv::Float64,frc_M::Matrix{Float64}=[0. 0. ; 0. 0.])
	
	if isempty(sheet.csv_vac)
		# Normalize distr_vac to probability mass for compartments
		# left with susceptibles
		distr_vac = (sheet.distr_vac).*(S .> 0);
		mass = sum(distr_vac);
		if mass == 0.
			return distr_vac
		end
		distr_vac = distr_vac/mass;

		# Compute effective vaccination rates
		Srate = sheet.vrate*distr_vac;

		# Determine if still vaccine available
		Srate = Srate.*(Nv <= sheet.vtot);

	else
		Srate = Vector{Float64}(undef,9);
		for i=1:9
			Srate[i] = myinterp(frc_M[:,1],frc_M[:,1+i],t);
		end

		# Determine if S left to vaccinate
		Srate = Srate.*(S .> 0.);
	end

	return Srate

end


#%% flowfield
"""
Define the flowfield to be used with the 4th order Runga-Kutta solver below
"""
function flowfield(sheet::data,mydep::Dict{Symbol,Vector{Float64}},
		   t::Float64,u::Vector{Float64},
		   frc_M::Matrix{Float64}=[0. 0. ; 0. 0.])

	# Reshape u from column vector into 3 arrays
	n = length(u);
	u = reshape(u,(27,Int64(n/27)));

	#  Separate into unvacc, vacc, and vacc hes
	#  X = [S E I R D]
	n = 9;
	X = u[1:n,:];
	Y = u[n+1:2n,:];
	Z = u[2n+1:3n,:];

	S = X[:,1]; E = X[:,2]; I = X[:,3]; R = X[:,4]; D = X[:,5];
	Sv = Y[:,1]; Ev = Y[:,2]; Iv = Y[:,3]; Rv = Y[:,4]; Dv = Y[:,5];
	Sx = Z[:,1]; Ex = Z[:,2]; Ix = Z[:,3]; Rx = Z[:,4]; Dx = Z[:,5];

	# Compute size of vaccinated population for possibile limit to vax
	Nv = sum(Y);

	# Compute ODE systems
	#  Unvaccinated
	DS = (-S.*mydep[:β]).*(sheet.C*(
			       (I+sheet.ω*Iv+Ix)./(sheet.N-D-Dv-Dx)
			       )
			    );
	DE = -DS - (1/sheet.d_E)*E;
	#   Now add in vaccination protocol
	DS = DS - vacprt(sheet,t,S,Nv,frc_M);	
	DI = (1/sheet.d_E)*E - (1/sheet.d_I)*I;
	DR = (1/sheet.d_I)*I.*(1 .- sheet.IFR);
	DD = (1/sheet.d_I)*I.*sheet.IFR;

	DX = [DS DE DI DR DD];

	#  Vaccinated
	DS = sheet.α*
	     (-Sv.*mydep[:β]).*(sheet.C*(
				 (I+sheet.ω*Iv+Ix)./(sheet.N-D-Dv-Dx)
				)
			     );
	DE = -DS - (1/sheet.d_E)*Ev;
	#   Now add in vaccination protocol
	DS = DS + vacprt(sheet,t,S,Nv,frc_M);	
	DI = (1/sheet.d_E)*Ev - (1/sheet.d_I)*Iv;
	DR = (1/sheet.d_I)*Iv.*(1 .- sheet.α*sheet.IFR);
	DD = (1/sheet.d_I)*Iv.*(sheet.α*sheet.IFR);

	DY = [DS DE DI DR DD];

	#  Unwilling to vaccinate
	DS = (-Sx.*mydep[:β]).*(sheet.C*(
				 (I+sheet.ω*Iv+Ix)./(sheet.N-D-Dv-Dx)
				)
			     );
	DE = -DS - (1/sheet.d_E)*Ex;
	DI = (1/sheet.d_E)*Ex - (1/sheet.d_I)*Ix;
	DR = (1/sheet.d_I)*Ix.*(1 .- sheet.α*sheet.IFR);
	DD = (1/sheet.d_I)*Ix.*(sheet.α*sheet.IFR);

	DZ = [DS DE DI DR DD];

	Du = [DX;DY;DZ];
	flow = Du[:];

	return flow
	
end

#%% runga
"""
This is the rungakutta part of the odesolver. It is separated to its own routine
to address the memory allocation issues that Julia can generate when working with
arrays
"""
function runga!(sheet::data,mydep::Dict{Symbol,Vector{Float64}},
		tnow::Float64,δt::Float64,ynow::Vector{Float64},
		frc_M::Matrix{Float64},
 	        K0::Array{Float64,1},K1::Array{Float64,1},K2::Array{Float64,1},K3::Array{Float64,1})
	
	K0[:] = δt*flowfield(sheet,mydep,tnow,ynow,frc_M);
	K1[:] = δt*flowfield(sheet,mydep,tnow+0.5*δt,ynow+0.5*K0,frc_M);
	K2[:] = δt*flowfield(sheet,mydep,tnow+0.5*δt,ynow+0.5*K1,frc_M);
	K3[:] = δt*flowfield(sheet,mydep,tnow+δt,ynow+K2,frc_M)
		
end


#%% odesolver
"""
An odesolver for the Bubar age-stratified model. 
Uses a 4th order Runga-Kutta solver taken from Hildebrand Introduction to Numerical Analysis
at (6.14.8)-(6.14.10)

Returns vector tpts and array ypts whose rows are sol at different time points.  All 27 x 5 
SEIR variables have been flattened into single column array like [:]. Use reshape to revert
to human readable format.
"""
function odesolver(sheet::data,mydep::Dict{Symbol,Vector{Float64}},
		   frc_M::Matrix{Float64}=[0. 0. ; 0. 0.])
	# Initialize initial data
	S = [mydep[:S0];mydep[:Sv0];mydep[:Su0]];
	E = [mydep[:E0];mydep[:Ev0];mydep[:Eu0]];
	I = [mydep[:I0];mydep[:Iv0];mydep[:Iu0]];
	R = [mydep[:R0];mydep[:Rv0];mydep[:Ru0]];
	D = [mydep[:D0];mydep[:Dv0];mydep[:Du0]];

	u0 = [S E I R D];
	u0 = u0[:];
	
	# Declare time span
	tspan = mydep[:tspan];
	
	# Initialize time axis
	ntpts = Int64(ceil((tspan[2]-tspan[1])/sheet.δt));
	tpts = [x for x in LinRange(tspan[1],tspan[2],ntpts)];
	tpts[end] = tspan[2];
	
	#  δt index goes to left end as tpts index
	δt = tpts[2:end]-tpts[1:end-1];

	nyval = length(u0);

	ypts = Array{Float64,2}(undef,nyval,ntpts);
	ypts[:,1] = u0;
	
	K0 = Array{Float64,1}(undef,nyval);
	K1 = Array{Float64,1}(undef,nyval);
	K2 = Array{Float64,1}(undef,nyval);
	K3 = Array{Float64,1}(undef,nyval);
	@inbounds for i=1:ntpts-1
		
		# Extract the time increment
		h = δt[i];
			
		# Define the vectorial slope values
		ynow = ypts[:,i];
		runga!(sheet,mydep,tpts[i],h,ynow,frc_M,K0,K1,K2,K3);	

		# Compute the new y-value by weighted average
		ypts[:,i+1] = ynow + 1/6*(K0+2*K1+2*K2+K3);	

	end
	
	# Transpose array
	ysol = Array{Float64,2}(undef,ntpts,nyval);
	@inbounds for j=1:ntpts
		for i=1:nyval
			ysol[j,i] = ypts[i,j];
		end
	end

	return tpts,ysol
end