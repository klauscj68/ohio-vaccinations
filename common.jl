# Routines for defining the model parameters needed to run the age
# stratified vaccination model
using CSV, DataFrames, Dates, LinearAlgebra

#%% datamat
"""
Declare a base set of parameter values for the model
"""
function datamat()
	mydata = Dict{Symbol,Any}()

	# Disease parameters
	#  Incubation period
	mydata[:d_E] = 3.; # days

	#  Infectious period
	mydata[:d_I] = 5.; # days

	# Epidemic transmission
	#  Relative susceptibility before r0 rescaling
	mydata[:β] = [0.4, 0.38, 0.79,
		      0.86, 0.8, 0.82,
		      0.88, 0.74, 0.74];
	#  Intended r0
	mydata[:r0] = 1.1620695433381327;

	#  Infection fatality rate
	mydata[:IFR] = 1e-2*[0.001, 0.003, 0.01,
			     0.04, 0.12, 0.40,
			     1.36, 4.55, 15.24];
	
	# Population parameters
	mydata[:N] = [695933. + 1/3*2233431,  # a<=5 + 1/3*5<=a<=19
		      2. /3*2233431,  # 2/3* 5<=a<=19
		      1554324.,  # 20<=a<=29,
		      1428941.,  # 30<=a<=39
		      691199. + 1/2*1544429,  # 40<=a<=44 + 1/2*45<=a<=54
		      1/2*1544429. + 834920,  # 1/2*45<=a<=54 + 55<=a<=59
		      765841. + 631247,  # sum 60's half decades
		      449394. + 331259,  # sum 70's half decades
		      228889. + 252072];  # sum 80's+

	# Vaccine parameters
	#  Vaccine susceptibility
	mydata[:α] = 0.22369101657958515;

	#  Vaccine contagiousness
	mydata[:ω] = 0.00028274830966965037;

	#  Total available vaccine
	mydata[:vtot] = Inf;

	#  Fraction of population unwilling to be vaccinated
	mydata[:vh] = 0.;

	#  Vaccination strategy
	#   Name of file if from csv. Empty denotes not.
	mydata[:csv_vac] = "All_by_Sept_7.csv";

	#   If not csv specify unnormalized age distr from dashboard
	mydata[:distr_vac] = [343799.0*.5,  # 0-9
			      343799.0*.5,   # 10-19
			      598957.,    # 20-29
			      673683.,   # 30-39
			      721343.,   # 40-49
			      930382.,   # 50-59
			      564371.0+527817.0,   # 60-69
			      437683.0+287414.0,   # 70-79
			      390603.0];    # 80+

	#   Constant vaccination rate if not from csv
	mydata[:vrate] = 24063.47; # 1997268 vax started Dec 15th - 
				   #                     March 8 (83 days)
	
	# Contact matrix
	C = CSV.read("contact_matrix_usa.csv",DataFrame);
	mydata[:C] = convert(Matrix{Float64},C[:,2:end]);
	
	# Initial conditions
	#  CSV to read initial conditions from, leave empty if not a csv
	mydata[:csv_odh] = "ODH_Data/ODH_0616.csv";

	#  Aggregate unvaccinated population
	mydata[:I0] = [  670.2065124033747
		         692.7934875966254
			  1405.0
			   1254.0
			    1124.0
			     1120.0
			       843.0
			         365.0
				   155.0];
	
	#  Aggregate initial exposed apartment before d_E normalization
	mydata[:E0] = [ 172.59169908553523
		        178.40830091446477
			 358.0
			  314.0
			   263.0
			    256.0
			     171.0
			       88.0
			         47.0];

	#  Aggregate deceased
	mydata[:D0] = [    2.950285454453595
		           3.049714545546405
			      22.0
			         96.0
				   252.0
				     915.0
				      2638.0
				       4959.0
				        9627.0];

	#  Aggregate recovered
	mydata[:R0] = [  61335.94288384782
		         63403.057116152166
			  176059.0
			   146471.0
			    140021.0
			     146513.0
			      112669.0
			        62862.0
				  41326.0];

	#  Aggregate vaccinated at time 0
	#   Although the key says Sv0, the code treats it as the aggregate
	#   of all vaccinated and then splits these across categories in 
	#   depmat
	mydata[:Sv0] = [  56325.05481156863
			  56325.05481156863
			   196255.86959108495
			    220740.79273425453
			     236357.1971584638
			      304851.48092750035
			       357869.2722465017
			        237587.2429454651
				 127986.03477359236];

	# ODE solver params
	#  Time span for simulation. Day is relative Jan 1, 2020
	mydata[:tspan] = [425., 508.]; # Run1: 515 Run2: 470

	#  Runga kutta time step
	mydata[:δt] = .25;

	return mydata

end

#%% auxmat
"""
Define a dictionary of auxiliary parameters used to redefine any params
with relationships to others
"""
function auxmat()
	mydata = Dict{Symbol,Float64}();
	
	# Factor used to scale ODH reported data 
	mydata[:rptλ] = 2.9359450218356584;	

	# Standard deviation for the likelihood 
	mydata[:bayσ] = 108.58090715085467;

	# Factor for an additional increase in reporting factor due to variants
	mydata[:rptinc] = 1.;

	# Weight for initial amount of infection among vaccinated
	mydata[:vι0] = 0.;

	# Factor used to additionally scale ODH reported exposure
	mydata[:rptλE] = 1.;

	# Factor used to additionally scale ODH reported infection
	mydata[:rptλI] = 1.

	return mydata
end

#%% csvdat
"""
Format the datmat and auxmat dictionaries into a standard vector
compatible for reading and writing csv files. One method converts
dictionaries into an array and the other converts an array into a 
dictionary.
"""
function csvdat(mydat::Dict{Symbol,Any},myaux::Dict{Symbol,Float64})
	datary = Vector{Float64}(undef,181);

	# Primary
	datary[1] = mydat[:d_E];
	datary[2] = mydat[:d_I];
	datary[3:11] = mydat[:β];
	datary[12] = mydat[:r0];
	datary[13:21] = mydat[:IFR];
	datary[22:30] = mydat[:N];
	datary[31] = mydat[:α];
	datary[32] = mydat[:ω];
	datary[33] = mydat[:vtot];
	datary[34] = mydat[:vh];
	datary[35] = NaN; # csv_vac
	datary[36:44] = mydat[:distr_vac];
	datary[45] = mydat[:vrate];
	datary[46:126] = mydat[:C][:];
	datary[127] = NaN; # csv_odh
	datary[128:136] = mydat[:I0];
	datary[137:145] = mydat[:E0];
	datary[146:154] = mydat[:D0];
	datary[155:163] = mydat[:R0];
	datary[164:172] = mydat[:Sv0];
	datary[173:174] = mydat[:tspan];
	datary[175] = mydat[:δt];

	# Auxiliary
	pos = 175; # Starting index for auxiliary params
	datary[pos+1] = myaux[:rptλ];
	datary[pos+2] = myaux[:bayσ];
	datary[pos+3] = myaux[:rptinc];
	datary[pos+4] = myaux[:vι0];
	datary[pos+5] = myaux[:rptλE];
	datary[pos+6] = myaux[:rptλI];

	return datary, pos
end
function csvdat(datary::Vector{Float64},csv_vac::String="",csv_odh::String="")
	mydat = Dict{Symbol,Any}();
	myaux = Dict{Symbol,Float64}();

	# Primary
	mydat[:d_E] = datary[1];
	mydat[:d_I] = datary[2];
	mydat[:β] = datary[3:11] ;
	mydat[:r0] = datary[12];
	mydat[:IFR] = datary[13:21];
	mydat[:N] = datary[22:30];
	mydat[:α] = datary[31];
	mydat[:ω] = datary[32];
	mydat[:vtot] = datary[33];
	mydat[:vh] = datary[34];
	mydat[:csv_vac] = csv_vac;
	mydat[:distr_vac] = datary[36:44];
	mydat[:vrate] = datary[45];
	mydat[:C] = reshape(datary[46:126],(9,9));
	mydat[:csv_odh] = csv_odh;
	mydat[:I0] = datary[128:136];
	mydat[:E0] = datary[137:145];
	mydat[:D0] = datary[146:154];
	mydat[:R0] = datary[155:163];
	mydat[:Sv0] = datary[164:172];
	mydat[:tspan] = datary[173:174];
	mydat[:δt] = datary[175];

	# Auxiliary
	pos = 175; # Starting index for auxiliary params
	myaux[:rptλ] = datary[pos+1];
	myaux[:bayσ] = datary[pos+2];
	myaux[:rptinc] = datary[pos+3];
	myaux[:vι0] = datary[pos+4];
	myaux[:rptλE] = datary[pos+5];
	myaux[:rptλI] = datary[pos+6];

	return mydat,myaux
end

# data
"""
Structure for encoding the independent model parameters in the Bubar
age-stratified model. Includes three inner constructors for defining
either by an array, by a dictionary, or directly by call to datamat()
"""
struct data
	# Disease progression
	d_E::Float64
	d_I::Float64

	# Epidemic transmissibility
	β::Vector{Float64}
	r0::Float64
	IFR::Vector{Float64}

	# Population
	N::Vector{Float64}

	# Vaccine 
	α::Float64
	ω::Float64
	vtot::Float64
	vh::Float64
	csv_vac::String
	distr_vac::Vector{Float64}	
	vrate::Float64

	# Contact matrix
	C::Matrix{Float64}

	# Initial conditions
	csv_odh::String
	I0::Vector{Float64}
	E0::Vector{Float64}
	R0::Vector{Float64}
	D0::Vector{Float64}
	Sv0::Vector{Float64}

	# ODE solver
	tspan::Vector{Float64}
	δt::Float64

	function data(datamat::Vector{Union{Float64,String}})
		d_E,d_I = datamat[1:2];
		β = datamat[3:11]; r0 = datamat[12]; 
			IFR = datamat[13:20];
		N = datamat[21:29];
		α,ω,vtot,vh = datamat[30:33];  
		csv_vac = datamat[34]; distr_vac = datamat[35:43];
			 vrate = datamat[44];
		C = reshape(datamat[45:125],(9,9));
		ram = reshape(datamat[126:170],(9,5)); 
			I0 = ram[:,1];
			E0 = ram[:,2]; 
			R0 = ram[:,3]; 
			D0 = ram[:,4]; 
			Sv0 = ram[:,5];
		tspan = datamat[171:172]; δt = datamat[173];
		csv_odh = datamat[174];

		return new(d_E,d_I,
			   β,r0,IFR,
			   N,
			   α,ω,vtot,vh,
			   csv_vac,distr_vac,vrate,
			   C,
			   csv_odh,I0,E0,R0,D0,Sv0,
			   tspan,δt)
	end

	function data(dict::Dict{Symbol,Any})

		return new(dict[:d_E],dict[:d_I],
			   dict[:β],dict[:r0],dict[:IFR],
			   dict[:N],
			   dict[:α],dict[:ω],dict[:vtot],dict[:vh],
			   dict[:csv_vac],dict[:distr_vac],dict[:vrate],
			   dict[:C],
			   dict[:csv_odh],dict[:I0],dict[:E0],dict[:R0],dict[:D0],dict[:Sv0],
			   dict[:tspan],dict[:δt])
	end

	function data()
		mydata = datamat();

		return data(mydata)
	end
end
"""
Revert to a dictionary using a data structure element
"""
function datamat(sheet::data)
# Restore a dictionary from the sheet, as dictionaries are mutable
	mydata = Dict{Symbol,Any}();
	for key in fieldnames(data)
		mydata[key] = getfield(sheet,key);
	end

	return mydata
end

#%% odhld
"""
Extract initial conditions from odh downloaded csv file.
You will then need to manually input this data into datamat() for example
"""
function odhld(sheet::data)
	println("Assumes files were output by parseODH.py with its naming"*
		" conventions ...")
	println("You should update the vaccination by age distribution in"*
		" distr_vac field")
	println("Manually input this data into the datamat file")
	
	# ODH Dashboard Vaccination
	odhvax = sheet.distr_vac; 
	odhvax = odhvax/sum(odhvax);


	fname = sheet.csv_odh;
	if isempty(fname)
		println("No csv specified for data");
		return
	end

	# Load the csv with ODH case data
	dfodh = CSV.read(fname,DataFrame);
	cols = ["0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"];

	# Locate the initial day in spreadsheet
	dates = dfodh[!,:time]; # Vector of dates
	day0 = Date("2020-01-01");
	dayi = day0 + Day(sheet.tspan[1]);
	pos = 0; flagfd = false;
	while !flagfd
		pos += 1;
		if dates[pos] == dayi
			break
		end
	end

	# Aggregate the initial infection (sum over last 5 days)
	#               initial exposed (unnormalized by d_E so can use
	#                                MCMC properly on d_E. d_E 
	#                                normalization done in depmat)
	#               total infection over pandemic (- deceased 
	#                                              = recovered)
	I0 = Vector{Float64}(undef,9);
	E0 = Vector{Float64}(undef,9);
	Cum0 = Vector{Float64}(undef,9);
	for k=1:9
		E0[k] = dfodh[pos,cols[k]];
		I0[k] = sum(dfodh[pos-4:pos,cols[k]]);
		Cum0[k] = sum(dfodh[1:pos-5,cols[k]]);
	end

	# Aggregate the total deceased
	dfdec = CSV.read(fname[1:end-4]*"_death.csv",DataFrame);
	dates = dfdec[!,:time];
	pos = 0; flagfd = false;
	while !flagfd
		pos += 1;
		if dates[pos] == dayi
			break
		end
	end
	
	# Aggregate the initial deceased
	D0 = Vector{Float64}(undef,9);
	for k=1:9
		D0[k] = sum(dfdec[1:pos-5,cols[k]]);
	end

	# Compute the initial recovered
	R0 = Cum0 - D0;

	# Compute all initial vaccinated
	dfvax = CSV.read(fname[1:end-4]*"_vax.csv",DataFrame);
	dates = dfvax[!,:time];
	pos = 0; flagfd = false;
	while !flagfd
		pos += 1;
		if dates[pos] == dayi
			break
		end
	end

	Sv0 = dfvax[pos,:cum_start]*odhvax;

	# Create a dataframe with the starting values
	df = DataFrame(:AGE=>cols,
		       :E0=>E0,:I0=>I0,:R0=>R0,:D0=>D0,:Sv0=>Sv0);
	println(df);
	return df
	
end
"""
When called with no arguments it loads the sheet defined by datamat
"""
function odhld()
	sheet = data();
	df = odhld(sheet);

	return df
end

#%% mymvavg
"""
Perform a moving average smoothing of a vector. This may be called by 
other libraries for computation of errors
"""
function mymvavg(v::Vector{Float64}; inc::Int64=7)
	n = length(v);
	inc = isodd(inc) ? inc : inc+1;
	half = Int64((inc-1)/2);
	
	vsmth = Vector{Float64}(undef,n);
	vpad = [repeat([v[1]],outer=(half,)); v; repeat([v[end]],outer=(half,))]
	for pos=1:n
		vsmth[pos] = sum(vpad[half+pos-half:half+pos+half])/inc;
	end

	return vsmth

end

#%% vaxld
"""
Import a vaccination schedule from a csv like what is output by Chance's
routines
"""
function vaxld(sheet::data)
	# Load the csv
	fname = sheet.csv_vac;
	if isempty(fname)
		return
	end

	df = CSV.read(fname,DataFrame);

	# Revert the dates to Jan 1, 2020 0-index and then to ti 0-index
	day0 = Date("2020-01-01");
	dates = [Date(x,"mm/dd/yy") + Year(2000) for x in df[!,:date]];
	ram = dates .- day0;
	taxis = [Float64(getfield(t,:value)) for t in ram];
	taxis = taxis .- sheet.tspan[1];

	# Build the M matrix
	M = Matrix{Float64}(undef,size(df)[1],10);
	M[:,1] = taxis;
	cols = ["0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"];
	for i=1:9
		M[:,1+i] = df[!,cols[i]];
	end

	return M
end
function vaxld()
	sheet = data();
	
	return vaxld(sheet)
end

#%% depmat
"""
Define a dictionary using the inputs of datamat() and auxmat() which
updates the values of parameters depending on others
"""
function depmat(sheet::data,auxmat::Dict{Symbol,Float64})	
		
	mydata = Dict{Symbol,Vector{Float64}}();
	
	# Initial conditions
	vaxdstr0 = sheet.Sv0./sheet.N;

	# E0 across groups
	agg = auxmat[:rptλE]*auxmat[:rptλ]*sheet.d_E*sheet.E0;
	mydata[:Ev0] = agg.*vaxdstr0*auxmat[:vι0];
	mydata[:E0] = (agg - mydata[:Ev0])*(1-sheet.vh);
	mydata[:Eu0] = (agg - mydata[:Ev0])*sheet.vh;

	# I0 across groups
	agg = auxmat[:rptλI]*auxmat[:rptλ]*sheet.I0;
	mydata[:Iv0] = agg.*vaxdstr0*auxmat[:vι0];
	mydata[:I0] = (agg - mydata[:Iv0])*(1-sheet.vh);
	mydata[:Iu0] = (agg - mydata[:Iv0])*sheet.vh;

	# R0 across groups
	agg = auxmat[:rptλ]*sheet.R0;
	mydata[:Rv0] = agg.*vaxdstr0*auxmat[:vι0];
	mydata[:R0] = (agg - mydata[:Rv0])*(1-sheet.vh);
	mydata[:Ru0] = (agg - mydata[:Rv0])*sheet.vh;

	# D0 across groups
	agg = sheet.D0;
	mydata[:Dv0] = agg.*vaxdstr0*auxmat[:vι0];
	mydata[:D0] = (agg - mydata[:Dv0])*(1-sheet.vh);
	mydata[:Du0] = (agg - mydata[:Dv0])*sheet.vh;

	# S0 across groups
	mydata[:Sv0] = sheet.Sv0 - mydata[:Ev0] - mydata[:Iv0] 
	                         - mydata[:Rv0] - mydata[:Dv0];
	agg = sheet.N - (mydata[:E0] + mydata[:Ev0] + mydata[:Eu0] +
			 mydata[:I0] + mydata[:Iv0] + mydata[:Iu0] + 
			 mydata[:R0] + mydata[:Rv0] + mydata[:Ru0] + 
			 mydata[:D0] + mydata[:Dv0] + mydata[:Du0] +
			 mydata[:Sv0]);
	mydata[:S0] = (1-sheet.vh)*agg;
	mydata[:Su0] = sheet.vh*agg;

	# Reproduction number
	M = sheet.β.*sheet.C*sheet.d_I;
	Mλ = eigvals(M);
	λmax = maximum(abs.(Mλ));
	mydata[:β] = sheet.r0/λmax*sheet.β;

	# tspan length of simulation
	mydata[:tspan] = [0.,sheet.tspan[2]-sheet.tspan[1]];

	return mydata
end
function depmat(mydata::Dict{Symbol,Any},auxmat::Dict{Symbol,Float64})
	sheet = data(mydata);

	return depmat(sheet,auxmat)
end
