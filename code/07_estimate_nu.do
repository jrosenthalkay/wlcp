**
** estimating fossil supply elasticities
**

**
** prepare int'l energy price data
**

* intl oil price
insheet using "$data/raw/crude-oil-prices.csv" , clear

tempfile p 
save `p'

* intl gas price
insheet using "$data/raw/natural-gas-prices.csv" , clear

collapse gasprice , by(year)

tempfile pg 
save `pg'

* intl coal price 
insheet using "$data/raw/coal-prices.csv" , clear 

* plot coal price time series
local colors `" "navy" "maroon" "forest_green" "dkorange" "purple" "teal" "cranberry" "olive_teal" "midblue" "gold" "'

levelsof code , local(clevs) 
local plotsring ""
local i = 0
foreach cc of local clevs { 
	local i = `i' + 1
	local colorindex = mod(`i', `:word count `colors'') + 1
    local colorname : word `colorindex' of `colors'
	local plotstring `"`plotstring' (line coalprice year if code == "`cc'" , lcol(`colorname') lpat(solid)) "'
}

graph twoway `plotstring' ///
	, xtit(Year) ytit("Price of coal (USD/tonne)") ///
	legend(on pos(6) rows(1) order(1 "AUS" 2 "COL" 3 "IDN" 4 "JPN" 5 "USA" 6 "ZAF")) ///
	tit("Coal prices across countries" , pos(11))
graph export "$figs/coalprice.png" , width(3200) replace
graph display
collapse coalprice , by(year)

tempfile pc 
save `pc'


**
** prepare estimation panel
**

* country data
use "$data/int/global_gdp_panel"  , clear 
 
keep iso3 year weighted_temp lgdp lpop
ren lgdp lgdp_pc
ren iso3 iso3
ren weighted_temp temp
tempfile gdp
save `gdp'

* start with nat'l accounts energy data panel
use "$data/int/wdi_nrg_o" , clear

* merge on int'l energy prices
ren iso3 iso3 
merge m:1 year using `p' , keep(1 3) nogen
merge m:1 year using `pg' , keep(1 3) nogen
merge m:1 year using `pc' , keep(1 3) nogen

* merge on GDP panel 
merge m:1 iso3 year using `gdp' , keep(1 3) nogen

* create temperature variables
bys year : egen Ncountry = total( !mi(temp) )
bys year : egen gtemp = mean(temp)
replace gtemp = (Ncountry)/(Ncountry-1) * gtemp - 1/Ncountry * temp 

bys iso3 : egen m_oilshare = mean(nrg_o_oil_ren)
bys iso3 : egen m_gasshare = mean(nrg_o_gas_ren)

gen totshare = (m_gasshare+m_oilshare)
gen fossilprice = m_gasshare/totshare * gasprice + m_oilshare/totshare * oilpricecrude

* variables
gen lprice_oil = ln( oilpricecrude )
gen lprice_gas = ln(gasprice)
gen lprice_coal = ln(coalprice)

gen lrent_oil = ln(nrg_o_oil_ren)
gen lrent_gas = ln(nrg_o_gas_ren)
gen lrent_pool = ln( nrg_o_gas_ren + nrg_o_oil_ren)
gen lfossilprice = ln(fossilprice)
gen lrent_coal = ln(nrg_o_coal_rents_share)

gen lgdp_tot = lgdp_pc + lpop

//
// fossil elasticities
//

encode iso3 , gen(iso3_n)
tsset iso3_n year 

* country specific estimates
gen nu = . 
gen nu_se = .

qui levelsof iso3 , local(clevs) 
foreach cc of local clevs {
	
	di in red "`cc'"
	count if !mi(D.lrent_pool) & !mi(D.lfossilprice) & !mi(D.temp) & !mi(D.lgdp_pc) & iso3 == "`cc'"
	
	if `r(N)'>3 {
		qui reg D.lrent_pool D.lfossilprice D.lgdp_tot if iso3 == "`cc'" , r
		qui replace nu = _b[D.lfossilprice] if iso3 == "`cc'"
		qui replace nu_se = _se[D.lfossilprice] if iso3 == "`cc'"
	}
	
}

* pooled estimate
reghdfe D.lrent_pool D.lfossilprice D.lgdp_tot , absorb(iso3_n) cluster(iso3_n year)
global nu = _b[D.lfossilprice]
global nu_se = _se[D.lfossilprice]

cap drop w0
gen w0 = (1/nu_se^2) / ((1/${nu_se}^2) + (1/nu_se^2))

cap drop nu_eb 
gen nu_eb = w0*nu + (1-w0)*${nu}

cap drop se_eb
gen se_eb = nu_se * sqrt(1 - w0)

cap drop nu_trunc
gen nu_trunc = nu_eb + se_eb * normalden( (1 - nu_eb) / se_eb ) / (1 - normal( (1 - nu_eb) / se_eb ) ) - 1 

replace nu = nu - 1
* replace missings with country averages
gen miflag = mi(nu_trunc)
replace nu_trunc = $nu if mi(nu_trunc)

graph twoway (hist nu if inrange(nu,-5,10), start(-5) w(0.25) lcol(%0) fcol(gs8%50)) ///
	(hist nu_trunc if !miflag , start(0) w(0.25) lcol(black) fcol(%0) ) ///
	, xtit("1/{&nu}{sub:i} (Oil-gas supply elasticity)") ytit("Density") legend(on pos(6) rows(1) order(1 "OLS estimate" 2 "Empirical Bayes estimate")) ///
	name(geb1, replace)
	
*
preserve 
collapse nu nu_trunc (max) miflag , by(iso3)

ren nu inv_nu_f
ren nu_trunc inv_nu_f_trunc
outsheet * using "$data/model_input/nu_fossil.csv" , comma replace 

drop if miflag 

sort inv_nu_f_trunc 
gen rank = _n 

gen labflag = ""
replace labflag = "Saudi Arabia" if iso3=="SAU"
replace labflag = "USA" if iso3=="USA"
replace labflag = "Russia" if iso3=="RUS"
replace labflag = "China" if iso3=="CHN"

local Nmax = _N + 1
graph twoway (scatter inv_nu_f_trunc rank, msym(O) mcol(ebblue) ) ///
	(scatter inv_nu_f_trunc rank if !mi(labflag) , msym(Oh) mcol(black) mlab(labflag) mlabcolor(black) mlabangle(30) ) ///
	, xtit("Supply elasticity rank") ytit("1/{&nu}{sup:f}{sub:i}") ///
	legend(off pos(6)) xsc(ra(0 `Nmax')) xlab(0(20)`Nmax') tit("Fossil supply elasticities", pos(11)) name(gf, replace) 
	
restore

//
//
// coal elasticities  
//
//
tsset iso3_n year 

* country specific estimates
cap drop nu_c 
cap drop nu_c_se
gen nu_c = . 
gen nu_c_se = .

qui levelsof iso3 , local(clevs) 
foreach cc of local clevs {
	
	di in red "`cc'"
	count if !mi(D.lrent_coal) & !mi(D.lrent_coal) & !mi(D.temp) & !mi(D.lgdp_pc) & iso3 == "`cc'"
	
	if `r(N)'>3 {
		qui reg D.lrent_coal D.lprice_coal D.lgdp_tot if iso3 == "`cc'" , r
		qui replace nu_c = _b[D.lprice_coal] if iso3 == "`cc'"
		qui replace nu_c_se = _se[D.lprice_coal] if iso3 == "`cc'"
	}
	
}

* pooled estimate
reghdfe D.lrent_coal D.lprice_coal D.lgdp_tot , absorb(iso3_n) cluster(iso3_n year)
global nu_c = _b[D.lprice_coal]
global nu_c_se = _se[D.lprice_coal]

cap drop w0
gen w0 = (1/nu_c_se^2) / ((1/${nu_c_se}^2) + (1/nu_c_se^2))

cap drop nu_c_eb 
gen nu_c_eb = w0*nu_c + (1-w0)*${nu_c}

cap drop se_eb
gen se_eb = nu_c_se * sqrt(1 - w0)

cap drop nu_c_trunc
gen nu_c_trunc = nu_c_eb + se_eb * normalden( (1- nu_c_eb ) / se_eb ) / (1 - normal(  (1- nu_c_eb) / se_eb ) ) - 1

replace nu_c = nu_c - 1

* replace missings with country averages
cap drop miflag
gen miflag = mi(nu_c_trunc)
replace nu_c_trunc = $nu_c if mi(nu_c_trunc)

graph twoway (hist nu_c if inrange(nu_c,-4,5), start(-4) w(0.25) lcol(%0) fcol(gs8%50)) ///
	(hist nu_c_trunc if !miflag , start(0) w(0.25) lcol(black) fcol(%0) ) ///
	, xtit("1/{&nu}{sub:i}{sup:c} (Coal supply elasticity)") ytit("Density") legend(on pos(6) rows(1) order(1 "OLS estimate" 2 "Empirical Bayes estimate")) ///
	name(geb2, replace)
	
preserve 
collapse nu_c nu_c_trunc (max) miflag , by(iso3)

ren nu_c inv_nu_c 
ren nu_c_trunc inv_nu_c_trunc 

outsheet * using "$data/model_input/nu_coal.csv" , comma replace 

drop if miflag 

sort inv_nu_c_trunc 
gen rank = _n 

gen labflag = ""
replace labflag = "Saudi Arabia" if iso3=="SAU"
replace labflag = "USA" if iso3=="USA"
replace labflag = "India" if iso3=="IND"
replace labflag = "Russia" if iso3=="RUS"
replace labflag = "China" if iso3=="CHN"

local Nmax = _N + 1
graph twoway (scatter inv_nu_c_trunc rank, msym(O) mcol(ebblue) ) ///
	(scatter inv_nu_c_trunc rank if !mi(labflag) , msym(Oh) mcol(black) mlab(labflag) mlabcolor(black) mlabangle(-30) ) ///
	, xtit("Supply elasticity rank") ytit("1/{&nu}{sup:c}{sub:i}") ///
	legend(off pos(6)) xsc(ra(0 `Nmax')) xlab(0(20)`Nmax') tit("Coal supply elasticities" , pos(11)) name(gc, replace)
	
restore

graph combine gf gc , xsiz(20) ysiz(10)
graph export "$figs/supply_elas_fig.png" , width(3200) replace 


graph combine geb1 geb2 , xcommon ycommon xsiz(20) ysiz(12) tit("Empirical Bayes shrinkage estimates" , pos(11))
graph export "$figs/eb_shrink.png" , width(3200) replace 
