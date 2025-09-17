
* load trade data
use "$data/int/trade_collapse_noNRG" , clear

* set data
bys iso3_o iso3_d (year) : gen yearN = 17-_N
replace yearN = yearN+1 if yearN>0

expand yearN 

bys iso3_o iso3_d (year) : replace year = 1999 + _n if mi(year)

* true zeros
replace trade = 0 if mi(trade) & iso3_o != iso3_d 

collapse (first) trade , by(iso3_o iso3_d year)

** get temperature data **
preserve 

use "$data/int/global_gdp_panel"  , clear 
 
keep iso3 year weighted_temp lgdp y_total_pop
ren lgdp lgdp_pc_o
ren iso3 iso3_o 
ren weighted_temp temp_o 
ren y_total_pop pop_o 
tempfile gdp_o
save `gdp_o'

restore 

** grab country level data **
preserve 

use "$data/int/global_gdp_panel"  , clear 
 
keep iso3 year weighted_temp lgdp
ren lgdp lgdp_pc_d
ren iso3 iso3_d 
ren weighted_temp temp_d 
tempfile gdp_d 
save `gdp_d'

restore 

* merge on country temperature 
merge m:1 iso3_o year using `gdp_o' , keep(1 3) nogen
merge m:1 iso3_d year using `gdp_d' , keep(1 3) nogen

merge m:1 iso3_o year using "$data/int/wdi_nrg_o" , keep(1 3) nogen
merge m:1 iso3_d year using "$data/int/wdi_nrg_d" , keep(1 3) nogen

* vars
egen double oy = group(iso3_o year)
egen double dy = group(iso3_d year)

gen jself_trade = trade if iso3_o == iso3_d 
bys iso3_d year : egen  self_trade = mean(jself_trade)

gen l_trade_rat = ln( trade / self_trade )

gen trade_rat = trade / self_trade

* set panel variables variables
egen double panel_id = group(iso3_o iso3_d)

* temperature variables
bys iso3_o : egen mtemp = mean(temp_o)

gen temp2_o = temp_o * temp_o
gen temp2_d = temp_d * temp_d

gen tdif = temp_o - temp_d 
gen tdif2 = temp2_o - temp2_d 

compress
save "$data/int/trade_estimation_panel" , replace