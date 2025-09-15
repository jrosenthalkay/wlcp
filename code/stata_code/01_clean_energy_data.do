**
** clean energy data
**

** Clean WDI energy 
insheet using "$data/raw/WDI_energy_data_extract.csv" , clear names

ren countrycode iso3_o 

keep iso3_o seriesname yr* 

drop if mi(iso3_o) 

reshape long yr , i(iso3_o seriesname) j(year)

destring yr , force replace 

gen varname = "_nrg_o_coal_rents_share" if seriesname =="Coal rents (% of GDP)"
replace varname = "_nrg_o_energy_use_kg_oil_pc" if seriesname =="Energy use (kg of oil equivalent per capita)"
replace varname = "_nrg_o_fossil_cons_perc" if seriesname == "Fossil fuel energy consumption (% of total)"
replace varname = "_nrg_o_gdp_per_unit_use_nom" if seriesname == "GDP per unit of energy use (PPP $ per kg of oil equivalent)"
replace varname = "_nrg_o_gdp_per_unit_use_real" if seriesname == "GDP per unit of energy use (constant 2021 PPP $ per kg of oil equivalent)"
replace varname = "_nrg_o_gas_rent_share" if seriesname == "Natural gas rents (% of GDP)"
replace varname = "_nrg_o_oil_rent_share" if seriesname == "Oil rents (% of GDP)"
replace varname = "_nrg_o_renewable_cons_perc" if seriesname == "Renewable energy consumption (% of total final energy consumption)"

drop seriesname 

reshape wide yr , i(iso3_o year) j(varname) string

foreach var of varlist yr* { 
	local mystr = subinstr("`var'","yr_","",.)
	ren `var' `mystr'
}

compress
save "$data/int/wdi_nrg_o" , replace


** repeat, exchange o-names with d- (destination)
insheet using "$data/raw/WDI_energy_data_extract.csv" , clear names

ren countrycode iso3_d

keep iso3_d seriesname yr* 

drop if mi(iso3_d) 

reshape long yr , i(iso3_d seriesname) j(year)

destring yr , force replace 

gen varname = "_nrg_d_coal_rents_share" if seriesname =="Coal rents (% of GDP)"
replace varname = "_nrg_d_energy_use_kg_oil_pc" if seriesname =="Energy use (kg of oil equivalent per capita)"
replace varname = "_nrg_d_fossil_cons_perc" if seriesname == "Fossil fuel energy consumption (% of total)"
replace varname = "_nrg_d_gdp_per_unit_use_nom" if seriesname == "GDP per unit of energy use (PPP $ per kg of oil equivalent)"
replace varname = "_nrg_d_gdp_per_unit_use_real" if seriesname == "GDP per unit of energy use (constant 2021 PPP $ per kg of oil equivalent)"
replace varname = "_nrg_d_gas_rent_share" if seriesname == "Natural gas rents (% of GDP)"
replace varname = "_nrg_d_oil_rent_share" if seriesname == "Oil rents (% of GDP)"
replace varname = "_nrg_d_renewable_cons_perc" if seriesname == "Renewable energy consumption (% of total final energy consumption)"

drop seriesname 

reshape wide yr , i(iso3_d year) j(varname) string

foreach var of varlist yr* { 
	local mystr = subinstr("`var'","yr_","",.)
	ren `var' `mystr'
}

compress
save "$data/int/wdi_nrg_d" , replace
