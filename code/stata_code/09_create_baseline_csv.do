// create full baseline for estimation 

* carbon intensity
insheet using "$data/raw/co2-intensity.csv" , clear 

ren AnnualCO carbon_intensity
ren code iso3 

collapse carbon_intensity , by(iso3 year)

tempfile c02_intensity
save `c02_intensity'

* carbon baseline 
insheet using "$data/raw/WDI_carbon.csv" , clear name

keep if strpos(seriesname,"excluding")>0

reshape long yr , i(countrycode countryname) j(year)

ren countrycode iso3
ren yr co2_em_mt
destring co2_em_mt , replace force 
keep if inrange(year,2000,2016)

collapse co2_em_mt , by(iso3)

tempfile co2em
save `co2em'

* energy trade shares
insheet using "$data/model_input/nrg_trade_share.csv" , clear

tempfile nrg_trade 
save `nrg_trade'

* supply elasticities 
insheet using "$data/model_input/nu_coal.csv" , clear 

ren miflag coal_miflag 
tempfile nuc
save `nuc'

* supply elasticities 
insheet using "$data/model_input/nu_fossil.csv" , clear 

ren miflag fossil_miflag
tempfile nuf
save `nuf'

* energy mix
insheet using "$data/raw/owid-energy-data.csv" , clear

gen ghg_total = greenhouse_gas_emissions
gen e_x = oil_production + gas_production
gen coal_share = coal_share_energy
gen fossil_share = fossil_share_energy - coal_share_energy
gen renewable_share = low_carbon_share_energy 

gen tot = coal_share + fossil_share + renewable_share

replace coal_share = coal_share/tot
replace fossil_share = fossil_share/tot
replace renewable_share = renewable_share/tot

ren iso_code iso3 
keep country iso3 year *_share gdp population electricity* e_x primary_energy_consumption ghg_total

merge m:1 iso3 year using `c02_intensity' , nogen keep(1 3)

drop if mi(iso3)
drop if year < 1985

merge 1:1 iso3 year using "$data/int/global_gdp_panel"  , keepusing(lgdp) keep(1 3)

ren lgdp lgdp_pc 
replace lgdp = ln(gdp/population) if mi(lgdp_pc)

gen lpop = ln(pop)
gen l_elec_demand = ln(electricity_demand)
gen l_elec_gen = ln(electricity_generation)

* extrapolate 
foreach var in coal_share fossil_share renewable_share carbon_intensity {
	nl (`var' = 1/(1+exp(-({b0} + {d1}*l_elec_demand + {d2}*l_elec_gen + {b1}*lgdp_pc + {b2}*lpop + {b3}*year)))) ///
	if !mi(`var') & !mi(lgdp_pc) & !mi(year) & !mi(l_elec_gen) & !mi(l_elec_demand)
predict `var'_hat 
}

egen hat_tot = rowtotal(*_hat)

foreach var in coal_share fossil_share renewable_share carbon_intensity {
	replace `var'_hat = `var'_hat / hat_tot 
}

* extrapolate - 2
foreach var in coal_share fossil_share renewable_share carbon_intensity {
	nl (`var' = 1/(1+exp(-({b0} + {b1}*lgdp_pc + {b2}*lpop + {b3}*year)))) ///
	if !mi(`var') & !mi(lgdp_pc) & !mi(year) 
predict `var'_hat2
}

egen hat_tot2 = rowtotal(*_hat2)

foreach var in coal_share fossil_share renewable_share carbon_intensity {
	replace `var'_hat2 = `var'_hat2 / hat_tot2 
}


collapse *_share *_hat *_hat2 e_x primary_energy_consumption ghg_total carbon_intensity if inrange(year,2000,2016) , by(iso3)

foreach var in coal_share fossil_share renewable_share carbon_intensity {

gen `var'_miflag = mi(`var')
replace `var' = `var'_hat if mi(`var')
replace `var' = `var'_hat2 if mi(`var')

}

drop *_hat* 

tempfile nrg_mix 
save `nrg_mix'

* rents 
use "$data/int/wdi_nrg_o" , clear

gen fossil_rent_share = (nrg_o_oil_rent + nrg_o_gas_rent)/100
gen coal_rent_share = (nrg_o_coal_rent)/100

collapse fossil_rent_share coal_rent_share if inrange(year,2000,2016) , by(iso3)
ren iso3 iso3 

tempfile rents
save `rents'

*
use "$data/int/global_gdp_panel"  , clear 

ren y_gdp_per_cap gdp_per_capita
ren y_total_pop population 
ren weighted_temp temp 
ren weighted_lat lat 
ren weighted_lon lon 

collapse temp gdp_per_capita population lat lon if inrange(year,2000,2016) , by(iso3)

merge 1:1 iso3 using `nrg_trade' , nogen keep(1 3)
merge 1:1 iso3 using `nuc' , nogen keep(1 3)
merge 1:1 iso3 using `nuf' , nogen keep(1 3)
merge 1:1 iso3 using `rents' , nogen keep(1 3)
merge 1:1 iso3 using `nrg_mix' , nogen keep(1 3)
merge 1:1 iso3 using `co2em' , nogen keep(1 3)

* drops 
drop inv_nu_c inv_nu_f

* new variables
gen py = gdp*(1-fossil_rent_share-coal_rent_share)
gen extraction_cost_share = fossil_rent_share * inv_nu_f_trunc + coal_rent_share * inv_nu_c_trunc
gen consumption_share = 1-extraction_cost_share

* missing data?
egen missing_flag = rowmiss(_all)

** extrapolate trade **
gen lgdp = ln(gdp_per_capita)
gen abslat = abs(lat)
reg nrg_export_share lgdp temp

nl (nrg_export_share = 1/(1+exp(-({b0} + {b1}*temp + {b2}*lgdp + {b3}*fossil_share + {b4}*extraction_cost + {b5}*abslat )))) ///
	if !mi(nrg_export_share) & !mi(lgdp) & !mi(fossil_share) & !mi(extraction_cost)
predict nrg_export_share_hat 

nl (nrg_import_share = 1/(1+exp(-({b0} + {b1}*temp + {b2}*lgdp + {b3}*fossil_share + {b4}*extraction_cost + {b5}*abslat )))) ///
	if !mi(nrg_export_share) & !mi(lgdp) & !mi(fossil_share) & !mi(extraction_cost)
predict nrg_import_share_hat 

* gen miflags 
gen nrg_import_share_miflag = mi(nrg_import_share)
replace nrg_import_share = nrg_import_share_hat if nrg_import_share_hat < nrg_export_share_hat & nrg_import_share_miflag

gen nrg_export_share_miflag = mi(nrg_export_share)
replace nrg_export_share = nrg_export_share_hat if nrg_import_share_hat > nrg_export_share_hat & nrg_export_share_miflag

replace nrg_export_share = 0 if mi(nrg_export_share) & !mi(nrg_import_share)
replace nrg_import_share = 0 if mi(nrg_import_share) & !mi(nrg_export_share)

drop nrg*hat 

* missing data - check 2
drop missing_flag
egen missing_flag = rowmiss(_all)

tab missing_flag

tab iso3 if missing_flag // all small

drop if missing_flag 

// qf 
gen qf = 371 / (41.868* 1e6) *3600 *1e9 /// TOE to KJ KJ to kWH, kwH to twh

drop ghg_total

outsheet * using "$data/model_input/baseline_csv_suff_stat.csv" , comma replace
