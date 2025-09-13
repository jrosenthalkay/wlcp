*
* Clean WDI GDP
*

insheet using "$data/raw/WDI_main.csv" , clear

ren v3 country_name 
ren v4 iso3

forv vv = 5/68 { 
	
	local y = 1955+`vv'
	ren v`vv' y`y'
}

drop if _n == 1
keep if inlist(v1,"Population, total","GDP per capita (current US$)","Net migration")

replace v1 = "total_pop" if v1 == "Population, total"
replace v1 = "gdp_per_capita" if v1 == "GDP per capita (current US$)"
replace v1 = "net_migration" if v1 == "Net migration"

drop v2

reshape long y , i(iso3 v1) j(year)

destring y , force replace
ren y y_

reshape wide y_ , i(iso3 year) j(v1) string

* merge on temperature
preserve
insheet using "C:\Users\l1aps02\FRB SF Dropbox\Aleisha Sawyer\WLCC_replication\data\int\country_avg_temp_timeseries_ghs_popweight.csv" , clear
destring weighted_temp , replace force 
destring weighted_lat , replace force 
destring weighted_lon , replace force 
drop if year < 1960 

collapse weighted_* , by(iso3 year)
tempfile x 
save `x'
restore 

merge 1:1 iso3 year using `x' , keep(1 3) nogen

encode iso3 , gen(iso3n)
xtset iso3n year

**
** interpolate GDP/cap for productivity estimation **
**

** first, interpolate population
gen lpop = ln(y_total_pop)
gen mipop = mi(y_total_pop)

* Get list of all country codes
levelsof iso3n, local(countries)

* Run regression for each country
foreach c of local countries {
	cap drop __poptemp__
	
    * Run regression for this country
    cap qui reg lpop year if iso3n==`c'
    
    * Predict values (including for missing observations)
    cap qui predict __poptemp__ if iso3n==`c', xb
    
    * Fill in missing values only
    cap qui replace y_total_pop = exp(__poptemp__) if missing(lpop) & iso3n==`c'
    
    * Drop prediction
    cap drop __poptemp__
}

** second, interpolate lgdp
gen lgdp = ln(y_gdp_per_cap)
gen migdp = mi(y_gdp_per_cap)

* Get list of all country codes
levelsof iso3n, local(countries)

* Run regression for each country
foreach c of local countries {
	cap drop __gdptemp__
	
    * Run regression for this country
    cap qui reg lgdp year if iso3n==`c'
    
    * Predict values (including for missing observations)
    cap qui predict __gdptemp__ if iso3n==`c', xb
    
    * Fill in missing values only
    cap qui replace y_gdp_per_cap = exp(__gdptemp__) if migdp & iso3n==`c'
    
    * Drop prediction
    cap drop __gdptemp__
}

compress 
save "$data/int/global_gdp_panel" , replace
