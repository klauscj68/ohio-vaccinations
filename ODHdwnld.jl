# Script to download and parse the ODH raw data into format that can be 
# used by software suite
using CSV, DataFrames, HTTP, Dates

#%% parseODH
"""
Download the ODH data and parse it for Bayesian fitting
"""
function parseODH(date::String)
	file1 = "ODHraw_"*date*".csv";
	file2 = "ODH_"*date*".csv";

	# Download the ODH data
	Traw = CSV.File(HTTP.get("https://coronavirus.ohio.gov/static/dashboards/COVIDDeathData_CountyOfResidence.csv").body);
	CSV.write(file1,Traw);
	Traw = CSV.read(Traw,DataFrame);

end
