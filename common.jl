# Routines for defining the model parameters needed to run the age
# stratified vaccination model
using CSV, DataFrames, Dates, LinearAlgebra, Plots

#%% datamat
"""
Declare a base set of parameter values for the model
"""
function datamat()
	mydata = Dict{Symbol,Any}()
	# Disease parameters
	#  Incubation period
	mydata[:d_E] = 4.; # days

	#  Infectious period
	mydata[:d_I] = 7.; # days

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
	mydata[:csv_vac] = "VaxCdcSep08_1D.csv";

	#   If not csv specify unnormalized age distr from dashboard
	mydata[:distr_vac] = [448430.0*.5,  # 0-9
			      448430.0*.5,   # 10-19
			      670676.0,    # 20-29
			      742184.0,   # 30-39
			      783253.0,   # 40-49
			      992897.0,   # 50-59
			      590849.0+546319.0,   # 60-69
			      450960.0+295134.0,   # 70-79
			      398997.0];    # 80+

	#   Constant vaccination rate if not from csv
	mydata[:vrate] = 24063.47; # 1997268 vax started Dec 15th - 
				   #                     March 8 (83 days)
	
	# Contact matrix
	C = CSV.read("contact_matrix_usa.csv",DataFrame);
	mydata[:C] = convert(Matrix{Float64},C[:,2:end]);
	
	# Initial conditions
	#  CSV to read initial conditions from, leave empty if not a csv
	mydata[:csv_odh] = "ODH_0909.csv";

	#  Aggregate unvaccinated population
	mydata[:I0] = [ 1026.6993381498505
		        1061.3006618501495
			 2097.0
			  1852.0
			   1693.0
			    1624.0
			     1245.0
			       549.0
			         242.0];
	
	#  Aggregate initial exposed apartment before d_E normalization
	mydata[:E0] = [ 176.52541302480668
		        182.47458697519332
			 361.0
			  312.0
			   264.0
			    257.0
			     177.0
			       86.0
			         45.0];

	#  Aggregate deceased
	mydata[:D0] = [    2.950285454453595
		           3.049714545546405
			      21.0
			         94.0
				   250.0
				     911.0
				      2627.0
				       4961.0
				        9628.0];

	#  Aggregate recovered
	mydata[:R0] = [61200.22975294295
		         63262.77024705704
			  175614.0
			   145975.0
			    139637.0
			     146195.0
			      112421.0
			        62756.0
				  41347.0];

	#  Aggregate vaccinated at time 0
	#   Although the key says Sv0, the code treats it as the aggregate
	#   of all vaccinated and then splits these across categories in 
	#   depmat
	mydata[:Sv0] = [  0.
			   1278.
			    132980.
			     122250.
			      125196.
			       137494.
			        325090.
				 296504.
				  182665. ];

	# ODE solver params
	#  Time span for simulation. Day is relative Jan 1, 2020
	mydata[:tspan] = [425., 610.]; 

	#  Runga kutta time step
	mydata[:δt] = .25;

	# Vax-Unvax interaction weight
	mydata[:ɾ] = 1.;

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
	mydata[:rptλI] = 1.;
	
	# Time point at which switch to change point parameters
	mydata[:Δpt] = Inf;

	# New change point r0 value
	mydata[:Δr0] = 2.
	
	# New change point α value
	mydata[:Δα] = .1;

	# New change point ω value
	mydata[:Δω] = .5;

	# New change point rptλ value
	mydata[:Δrptλ] = 2.;

	# New change point bayσ value
	mydata[:Δbayσ] = 108.58090715085467;

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
	datary = Vector{Float64}(undef,188);

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

	# Auxiliary and ɾ
	pos = 175; # Starting index for auxiliary params
	datary[pos+1] = myaux[:rptλ];
	datary[pos+2] = myaux[:bayσ];
	datary[pos+3] = myaux[:rptinc];
	datary[pos+4] = myaux[:vι0];
	datary[pos+5] = myaux[:rptλE];
	datary[pos+6] = myaux[:rptλI];
	datary[pos+7] = myaux[:Δpt];
	datary[pos+8] = myaux[:Δr0];
	datary[pos+9] = myaux[:Δα];
	datary[pos+10] = myaux[:Δω];
	datary[pos+11] = myaux[:Δrptλ];
	datary[pos+12] = myaux[:Δbayσ];
	
	datary[pos+13] = mydat[:ɾ];

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
	myaux[:Δpt] = datary[pos+7];
	myaux[:Δr0] = datary[pos+8];
	myaux[:Δα] = datary[pos+9];
	myaux[:Δω] = datary[pos+10];
	myaux[:Δrptλ] = datary[pos+11];
	myaux[:Δbayσ] = datary[pos+12];

	mydat[:ɾ] = datary[pos+13];

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

	# Interaction between vax and unvax weight
	ɾ::Float64

	function data(datamat::Vector{Union{Float64,String}})
		d_E,d_I = datamat[1:2];
		β = datamat[3:11]; r0 = datamat[12]; 
			IFR = datamat[13:21];
		N = datamat[22:30];
		α,ω,vtot,vh = datamat[31:34];  
		csv_vac = datamat[35]; distr_vac = datamat[36:44];
			 vrate = datamat[45];
		C = reshape(datamat[46:126],(9,9));
		csv_odh = datamat[127];
		ram = reshape(datamat[128:172],(9,5)); 
			I0 = ram[:,1];
			E0 = ram[:,2]; 
			R0 = ram[:,3]; 
			D0 = ram[:,4]; 
			Sv0 = ram[:,5];
		tspan = datamat[173:174]; δt = datamat[175];
		ɾ = datamat[188];

		return new(d_E,d_I,
			   β,r0,IFR,
			   N,
			   α,ω,vtot,vh,
			   csv_vac,distr_vac,vrate,
			   C,
			   csv_odh,I0,E0,R0,D0,Sv0,
			   tspan,δt,
			   ɾ)
	end

	function data(dict::Dict{Symbol,Any})

		return new(dict[:d_E],dict[:d_I],
			   dict[:β],dict[:r0],dict[:IFR],
			   dict[:N],
			   dict[:α],dict[:ω],dict[:vtot],dict[:vh],
			   dict[:csv_vac],dict[:distr_vac],dict[:vrate],
			   dict[:C],
			   dict[:csv_odh],dict[:I0],dict[:E0],dict[:R0],dict[:D0],dict[:Sv0],
			   dict[:tspan],dict[:δt],
			   dict[:ɾ])
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
Extract initial conditions from odh and cdc downloaded csv files.
You will then need to manually input this data into datamat() for example
"""
function odhld(sheet::data)
	println("Assumes files were output by parseODH.py with its naming"*
		" conventions and parsevax in common.jl ...");
	println("Manually input this data into the datamat file");


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

	# Aggregate the initial infection (sum over last 7 days)
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
		I0[k] = sum(dfodh[pos-6:pos,cols[k]]);
		Cum0[k] = sum(dfodh[1:pos-7,cols[k]]);
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
		D0[k] = sum(dfdec[1:pos-7,cols[k]]);
	end

	# Compute the initial recovered
	R0 = Cum0 - D0;

	# Compute all initial vaccinated
	dfvax = CSV.read(sheet.csv_vac,DataFrame);
	dates = dfvax[!,:Date];
	pos = 0; flagfd = false;
	while !flagfd
		pos += 1;
		if dates[pos] == dayi
			break
		end
	end
	
	Sv0 = sum(convert(Matrix,dfvax[1:pos,2:end]),dims=1)[:];

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
routines time 0 matches sheet.tspan[1]
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
	dates = df[!,:Date];	
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

	# Directly copy the change point parameters from auxmat
	mydata[:Δpt] = [auxmat[:Δpt]];
	mydata[:Δr0] = [auxmat[:Δr0]];
	mydata[:Δα] = [auxmat[:Δα]];
	mydata[:Δω] = [auxmat[:Δω]];
	mydata[:Δrptλ] = [auxmat[:Δrptλ]];

	return mydata
end
function depmat(mydata::Dict{Symbol,Any},auxmat::Dict{Symbol,Float64})
	sheet = data(mydata);

	return depmat(sheet,auxmat)
end

#%% Process cdc vaccine into time series age distribution
# parsevax
"""
Turn the cdc spreadsheet into age stratified data for vaccine delivery. 
Routine extrapolates vaccination rates into the future by taking the 
average rate over last two weeks and extrapolating

Data available at
https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-County/8xkx-amqh
"""
function parsevax(fname::String;
	          nfutday::Int64=60)
	# Read in CDC data and reduce to Ohio, sort by date
	df = CSV.read(fname,DataFrame);
	gdf = groupby(df,:Recip_State); df = DataFrame(gdf[("OH",)]);
	df[!,:Date] = Date.(df[!,:Date],"m/d/y"); sort!(df,:Date); 
	
	# CDC total doses by date should increase w/in each county so fill in missings
	# respecting mon inc w/in each county
	gdf = groupby(df,:Recip_County);
	for g in gdf
		for key in names(g)
			if nonmissingtype(eltype(g[!,key]))<:Real
				for i=1:size(g)[1]
					val = g[!,key][i];
					g[!,key][i] = ismissing(val) ? (i>1 ? g[!,key][i-1] : 0) : val;
				end
			end
		end
	end
	
	# Stack the subdataframes into one dataframe for future split by date
	dftemp = similar(df,0);
	for g in gdf
		dftemp = [dftemp;DataFrame(g)];
	end

	# Compute data columns marginalized over each date.
	#  Note: you can also use [:col1,:col2,:col3]=>((x,y,z)->x.*y-z)=>:new_col for ex
	#  in combine for other ways to compute columns
	gdf = groupby(df,:Date);
	df = combine(gdf,:Series_Complete_12Plus=>sum=>"2D12+",:Series_Complete_18Plus=>sum=>"2D18+",
		         :Series_Complete_65Plus=>sum=>"2D65+",
			 :Administered_Dose1_Recip_12Plus=>sum=>"1D12+",
			 :Administered_Dose1_Recip_18Plus=>sum=>"1D18+",
			 :Administered_Dose1_Recip_65Plus=>sum=>"1D65+");
	transform!(df,["1D12+","1D18+"]=>(-)=>"1D12-18",["1D18+","1D65+"]=>(-)=>"1D18-65",
		      ["2D12+","2D18+"]=>(-)=>"2D12-18",["2D18+","2D65+"]=>(-)=>"2D18-65");
	
	# Go from cumulative doses to daily doses administered
	for i=2:ncol(df)
		df[2:end,i] = df[2:end,i]-df[1:end-1,i];
	end
	
	# Make plots for inspecting extracted CDC vaccine schedule
	plot(df[!,:Date],[df[!,"1D12+"],df[!,"1D12-18"],df[!,"1D18-65"],df[!,"1D65+"]],
	     labels=["12+" "12-18" "18-65" "65+"],xlabel="date",ylabel="daily doses",
	     title="CDC: Administered First Vaccines");
	savefig("CDCvax_1dose.pdf");

	plot(df[!,:Date],[df[!,"2D12+"],df[!,"2D12-18"],df[!,"2D18-65"],df[!,"2D65+"]],
	     labels=["12+" "12-18" "18-65" "65+"],xlabel="date",ylabel="daily doses",
	     title="CDC: Administered Second Vaccines");
	savefig("CDCvax_2dose.pdf");

	# Shift dates by 14 days for allowing vaccine efficacy to take effect
	df[!,:Date] += Day(14);
	dayf = df[end,:Date];


	# Split up CDC like Ohio Population counts
	mydata = datamat();
	wt2065 = deepcopy(mydata[:N][3:7]); wt2065[end] *= .5;
	wt2065 *= 1/sum(wt2065);
	wt65p = deepcopy(mydata[:N][7:9]); wt65p[1] *= .5;
	wt65p *= 1/sum(wt65p);

	df1D = combine(df,"Date"=>(x->x)=>"Date",
			    "1D12-18"=>(x->x)=>"10-19","1D18-65"=>(x->wt2065[1]*x)=>"20-29",
			    "1D18-65"=>(x->wt2065[2]*x)=>"30-39","1D18-65"=>(x->wt2065[3]*x)=>"40-49",
			    "1D18-65"=>(x->wt2065[4]*x)=>"50-59",
			    ["1D18-65","1D65+"]=>((x,y)->wt2065[end]*x+wt65p[1]*y)=>"60-69",
			    "1D65+"=>(x->wt65p[2]*x)=>"70-79","1D65+"=>(x->wt65p[3]*x)=>"80+");
	df1D[!,"0-9"] = zeros(nrow(df1D));

	df2D = combine(df,"Date"=>(x->x)=>"Date",
			    "2D12-18"=>(x->x)=>"10-19","2D18-65"=>(x->wt2065[1]*x)=>"20-29",
			    "2D18-65"=>(x->wt2065[2]*x)=>"30-39","2D18-65"=>(x->wt2065[3]*x)=>"40-49",
			    "2D18-65"=>(x->wt2065[4]*x)=>"50-59",
			    ["2D18-65","2D65+"]=>((x,y)->wt2065[end]*x+wt65p[1]*y)=>"60-69",
			    "2D65+"=>(x->wt65p[2]*x)=>"70-79","2D65+"=>(x->wt65p[3]*x)=>"80+");
	df2D[!,"0-9"] = zeros(nrow(df2D));

	select!(df1D,[1,10,2,3,4,5,6,7,8,9]);
	select!(df2D,[1,10,2,3,4,5,6,7,8,9]);

	# Convert to integer values
	for i=2:ncol(df1D)
		df1D[!,i] = Int64.(floor.(df1D[!,i]));
		df2D[!,i] = Int64.(floor.(df2D[!,i]));
	end

	# Plot age-stratified vaccinations
	plot();
	for key in names(df1D)
		if key!="Date"
			plot!(df1D[!,:Date],df1D[!,key],labels=key,xlabel="date",ylabel="doses administered",
			      title="Single Dose");
		end
	end
	savefig("Vax1D_age.pdf");

	plot();
	for key in names(df2D)
		if key!="Date"
			plot!(df2D[!,:Date],df2D[!,key],labels=key,xlabel="date",ylabel="doses administered",
			      title="Second Dose");
		end
	end
	savefig("Vax2D_age.pdf");

	# Extrapolate into the future using constant rate over the
	# last two weeks
	vrate1D = Int64.(floor.(sum(convert(Matrix,df1D[end-13:end,2:end]),dims=1)/14));
	vrate2D = Int64.(floor.(sum(convert(Matrix,df2D[end-13:end,2:end]),dims=1)/14));

	for i=1:nfutday
		push!(df1D,[dayf+Day(i) vrate1D]);
		push!(df2D,[dayf+Day(i) vrate2D]);
	end

	# Write CSV's
	CSV.write(fname[1:end-4]*"_2D.csv",df2D);
	CSV.write(fname[1:end-4]*"_1D.csv",df1D);
	
end
