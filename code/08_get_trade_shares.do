* rents 
use "$data/int/wdi_nrg_o" , clear

gen fossil_rent_share = (nrg_o_oil_rent + nrg_o_gas_rent)/100
gen coal_rent_share = (nrg_o_coal_rent)/100

collapse fossil_rent_share if inrange(year,2000,2016) , by(iso3)
ren iso3 iso3 

tempfile rents
save `rents'

use "$data/int/global_gdp_panel"   , clear

keep iso3 year y_*
gen gdp = y_gdp_per_cap + y_total_pop

keep if inrange(year,2000,2016)

collapse gdp , by(iso3)

merge 1:1 iso3 using `rents' , keep(1 3) nogen

gen gdp_adj = gdp*(1-fossil_rent_share)

tempfile gdp 
save `gdp'

* Start with your initial data processing
use "$data/int/trade_collapse_noNRG", clear
bys iso3_d year: egen total_trade = total(trade)
gen trade_shares = trade/total_trade
sort iso3_d year iso3_o

collapse (mean) trade_shares (mean) trade , by(iso3_o iso3_d)

tempfile trade_s 
save `trade_s'
* Get unique list of all countries that appear as either origin or destination
preserve
keep iso3_o
duplicates drop
tempfile origins
save `origins'
restore

preserve
keep iso3_d
duplicates drop
rename iso3_d iso3_o
tempfile destinations
save `destinations'
restore

* Combine the unique country lists
preserve
use `origins', clear
append using `destinations'
duplicates drop
rename iso3_o iso3
tempfile allcountries
save `allcountries'
restore

* Create all possible country pairs
use `allcountries', clear
rename iso3 iso3_o
cross using `allcountries'
rename iso3 iso3_d

* Merge with the original data
merge 1:1 iso3_o iso3_d using `trade_s', keep(1 3)

* merge on importer gdp 
ren iso3_d iso3
merge m:1 iso3 using `gdp' , nogen keep(1 3)
ren iso3 iso3_d

bys iso3_d : egen tot_trade = total(trade)

* normalize to US 
su tot_trade if iso3_d == "USA"
local t1 = `r(mean)'
su gdp_adj if iso3_d == "USA"
local t2 = `r(mean)'
gen adj = `t1'/`t2'

gen gdp_adj2 = gdp_adj * adj

replace trade = gdp_adj2 - tot_trade if mi(trade) & iso3_o == iso3_d

drop trade_shares 
drop tot_trade
bys iso3_d : egen tot_trade = total(trade)
gen trade_shares = trade/tot_trade 
replace trade_shares = 0 if mi(trade_shares)

compress 
outsheet * using "$data/model_input/trade_shares.csv" , comma replace
