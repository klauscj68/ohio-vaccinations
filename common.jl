# Routines for defining the model parameters needed to run the age
# stratified vaccination model
using CSV, DataFrames, StaticArrays, Dates

# data
"""
Structure for encoding the independent model parameters in the Bubar
age-stratified model. Includes two inner constructors for defining
either by an array or by a dictionary.
"""
struct data
	# Disease progression
	d_E::Float64
	d_I::Float64

	# Epidemic transmissibility
	β::SVector{9,Float64}
	r0::Float64
	IFR::SVector{9,Float64}

	# Population
	N::SVector{9,Float64}

	# Vaccine 
	α::Float64
	ω::Float64
	vtot::Float64
	vh::Float64
	csv_vac::String
	distr_vac::SVector{9,Float64}	
	vrate::Float64

	# Contact matrix
	C::SMatrix{9,9,Float64}

	# Initial conditions
	csv_odh::String
	I0::SVector{9,Float64}
	E0::SVector{9,Float64}
	R0::SVector{9,Float64}
	D0::SVector{9,Float64}
	Sv0::SVector{9,Float64}

	# ODE solver
	tspan::Vector{Float64}
	δt::Float64

	function data(datamat::Vector{Union{Float64,String}})
		d_E,d_I = datamat[1:2];
		β = SVector{9,Float64}(datamat[3:11]); r0 = datamat[12]; 
			IFR = SVector{9,Float64}(datamat[13:20]);
		N = SVector{9,Float64}(datamat[21:29]);
		α,ω,vtot,vh = datamat[30:33];  
		csv_vac = datamat[34]; distr_vac = SVector{9,Float64}(datamat[35:43]);
			 vrate = datamat[44];
		C = SMatrix{9,9,Float64}(reshape(datamat[45:125],(9,9)));
		ram = reshape(datamat[126:170],(9,5)); 
			I0 = SVector{9,Float64}(ram[:,1]);
			E0 = SVector{9,Float64}(ram[:,2]); 
			R0 = SVector{9,Float64}(ram[:,3]); 
			D0 = SVector{9,Float64}(ram[:,4]); 
			Sv0 = SVector{9,Float64}(ram[:,5]);
		tspan = datamat[171:172]; δt = datamat[173];
		csv_odh = datamat[174];

		return new(d_E,d_I,
			   β,r0,IFR,
			   N,
			   α,ω,vtot,vh,vrate,
			   csv_vac,distr_vac,vrate,
			   C,
			   csv_odh,I0,E0,R0,D0,Sv0,
			   tspan,δt)
	end

	function data(dict::Dict{Symbol,Any})
		# Promote to static arrays if needed
		for key in [:β,:IFR,:N,:distr_vac,:I0,:E0,:R0,:D0,:Sv0]
			if !(typeof(dict[key]) <: SVector)
				dict[key] = SVector{9,Float64}(dict[key]);
			end
		end

		if !(typeof(dict[:C]) <: SMatrix)
			dict[:C] = SMatrix{9,9,Float64}(dict[:C]);
		end	

		return new(dict[:d_E],dict[:d_I],
			   dict[:β],dict[:r0],dict[:IFR],
			   dict[:N],
			   dict[:α],dict[:ω],dict[:vtot],dict[:vh],dict[:vrate],
			   dict[:csv_vac],dict[:distr_vac],distr[:vrate],
			   dict[:C],
			   dict[:csv_odh],dict[:I0],dict[:E0],dict[:R0],dict[:D0],dict[:Sv0],
			   dict[:tspan],dict[:δt])
	end
end

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
	mydata[:csv_vac] = "";

	#   If not csv specify unnormalized age distr
	mydata[:distr_vac] = [8872*.5,  # 0-9
			      8872*.5,   # 10-19
			      115479.,    # 20-29
			      163733.,   # 30-39
			      175522.,   # 40-49
			      204994.,   # 50-59
			      157069+317842.,   # 60-69
			      307675+220906.,   # 70-79
			      325176.];    # 80+

	#   Constant vaccination rate if not from csv
	mydata[:vrate] = 24063.47; # 1997268 vax started Dec 15th - 
				   #                     March 8 (83 days)
	
	# Contact matrix
	C = CSV.read("contact_matrix_usa.csv",DataFrame);
	mydata[:C] = convert(Matrix{Float64},C[:,2:end]);
	
	# Initial conditions
	#  CSV to read initial conditions from, leave empty if not a csv
	mydata[:csv_odh] = "";

	#  Aggregate unvaccinated population
	mydata[:I0] = [2421.693,
		       2503.307,
		       6748.,
		       6342.,
		       6070.,
		       6259.,
		       5049.,
		       3054.,
		       2477.];
	
	#  Aggregate initial exposed apartment
	mydata[:E0] = mydata[:d_E]*[525.1508,
				    542.8492,
				    1538.,
				    1395.,
				    1364.,
				    1333.,
				    1149.,
				    693.,
				    585.];

	#  Aggregate deceased
	mydata[:D0] = [.983428,
		       1.016572,
		       10.,
		       61.,
		       136.,
		       510.,
		       1317.,
		       2569.,
		       5115.];

	#  Aggregate recovered
	mydata[:R0] = [34183.48,
		       35335.52,
		       108015.,
		       86467.,
		       82821.,
		       88559.,
		       68793.,
		       41240.,
		       32412.] - mydata[:D0];

	#  Aggregate vaccinated at time 0
	mydata[:Sv0] = zeros((9,));

	# ODE solver params
	#  Time span for simulation
	mydata[:tspan] = [0., 135.];

	#  Runga kutta time step
	mydata[:δt] = .1;

	return mydata

end # Initial conditions for depmat
function datamat(sheet::data)
# Restore a dictionary from the sheet, as dictionaries are mutable
	mydata = Dict{Symbol,Any}();
	for key in fieldnames(data)
		mydata[key] = getfield(sheet,key);
	end

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

	return mydata
end

#%% depmat
"""
Define a dictionary using the inputs of datamat() and auxmat() which
updates the values of parameters depending on others
"""
function depmat(sheet::data,auxmat::Dict{Symbol,Float64})
	mydata = Dict{Symbol,Union{SVector{9,Float64},SMatrix{9,9,Float64}}}();
	
	# Check if initial conditions are loaded directly from ODH csv
	if !isempty(sheet.csv_odh)
		dfodh = CSV.read(sheet.csv_odh,DataFrame);
		# Convert days to time values: 0 is first day of csv
		time 
		
		
		
		# Unwilling to vaccinate
		λ = sheet.vh/(1-sheet.vh);
		mydata[:Eu0] = λ*sheet.E0;
		mydata[:Iu0] = λ*sheet.I0;
		mydata[:Ru0] = λ*sheet.R0;
		mydata[:Du0] = λ*sheet.D0;
	
		# Initial susceptible pool
		mydata[:S0] = sheet.N - (sheet.E0+mydata[:Eu0] + 
					 sheet.I0+mydata[:Iu0] + 
					 sheet.R0+mydata[:Ru0] + 
					 sheet.D0+mydata[:Du0]);
end
