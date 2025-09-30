**
** damage estimation
**

use "$data/int/trade_estimation_panel" , clear

** Estimation ** 

* energy controls
global nrg_control_o "nrg_o_oil_rent_share nrg_o_renewable_cons_perc"
global nrg_control_d "nrg_d_oil_rent_share nrg_d_renewable_cons_perc"

**
** Reg 1 poisson regression, o/d diff temps **
**
ppmlhdfe trade_rat temp_o temp2_o temp_d temp2_d ///
	c.( lgdp_pc_o $nrg_control_o )##c.( lgdp_pc_o $nrg_control_o ) ///
	c.( lgdp_pc_d $nrg_control_d )##c.( lgdp_pc_d $nrg_control_d ) ///
	, absorb(panel_id year) cluster(iso3_o iso3_d year)
eststo y_nrg_poisson_both , noesample 
	
* compute precision weighted avg for slope of damage function, delta method SE
nlcom -2 *( ( _se[temp2_o]^(-2) / (_se[temp2_o]^(-2) + _se[temp2_d]^(-2)) ) * _b[temp2_o] ///
	+ 2*( _se[temp2_d]^(-2) / (_se[temp2_o]^(-2) + _se[temp2_d]^(-2)) ) * _b[temp2_d] ) 
estadd scalar gam = r(b)[1,1] / 4 
local gse : display %5.4fc sqrt(r(V)[1,1]) / 4
estadd local gam_se = "$ ( `gse' ) $"

* compute precision weighted average of optimal temperature
nlcom _b[temp_o] / ( -2 * _b[temp2_o] )
local tstar1_se = sqrt(r(V)[1,1])

nlcom _b[temp_d] / ( -2 * _b[temp2_d] )
local tstar2_se = sqrt(r(V)[1,1])

local ow = `tstar1_se'^(-2) / ( `tstar1_se'^(-2) + `tstar2_se'^(-2) )

nlcom `ow' * _b[temp_o] / ( -2 * _b[temp2_o] ) + (1-`ow') * _b[temp_d] / ( -2 * _b[temp2_d] )
estadd scalar tstar = r(b)[1,1]
local tse : display %5.4fc sqrt(r(V)[1,1])
estadd local tstar_se = "$ ( `tse' ) $"

* table formatting -- fixed effect checkmarks 
estadd local pair_fe "\checkmark"
estadd local origin_fp "\checkmark"

**
** Reg 2: temperature differences as regressors
**

* estimate poisson model
ppmlhdfe trade_rat tdif tdif2 ///
	c.( lgdp_pc_o $nrg_control_o )##c.( lgdp_pc_o $nrg_control_o ) ///
	c.( lgdp_pc_d $nrg_control_d )##c.( lgdp_pc_d $nrg_control_d ) ///
	, absorb(panel_id year) cluster(iso3_o iso3_d year)
eststo y_nrg_poisson_one , noesample 

* compute gamma, tstar and their SE
nlcom -2 * _b[tdif2]
estadd scalar gam = r(b)[1,1] / 4 
local gse : display %5.4fc sqrt(r(V)[1,1]) / 4
estadd local gam_se = "$ ( `gse' ) $"

nlcom _b[tdif] / ( -2 * _b[tdif2] )
estadd scalar tstar = r(b)[1,1]
local tse : display %5.4fc sqrt(r(V)[1,1])
estadd local tstar_se = "$ ( `tse' ) $"

estadd local pair_fe "\checkmark"
estadd local origin_fp "\checkmark"

**
** Reg 3: repeat Reg 1, but OLS, not poisson
**
reghdfe l_trade_rat temp_o temp2_o temp_d temp2_d ///
	c.( lgdp_pc_o $nrg_control_o )##c.( lgdp_pc_o $nrg_control_o ) ///
	c.( lgdp_pc_d $nrg_control_d )##c.( lgdp_pc_d $nrg_control_d ) ///
	, absorb(panel_id year) cluster(iso3_o iso3_d year)
eststo y_nrg_ols_both , noesample 
	
*
nlcom -2 *( ( _se[temp2_o]^(-2) / (_se[temp2_o]^(-2) + _se[temp2_d]^(-2)) ) * _b[temp2_o] ///
	+ 2*( _se[temp2_d]^(-2) / (_se[temp2_o]^(-2) + _se[temp2_d]^(-2)) ) * _b[temp2_d] ) 
estadd scalar gam = r(b)[1,1] / 4 
local gse : display %5.4fc sqrt(r(V)[1,1]) / 4
estadd local gam_se = "$ ( `gse' ) $"

nlcom _b[temp_o] / ( -2 * _b[temp2_o] )
local tstar1_se = sqrt(r(V)[1,1])

nlcom _b[temp_d] / ( -2 * _b[temp2_d] )
local tstar2_se = sqrt(r(V)[1,1])

local ow = `tstar1_se'^(-2) / ( `tstar1_se'^(-2) + `tstar2_se'^(-2) )

nlcom `ow' * _b[temp_o] / ( -2 * _b[temp2_o] ) + (1-`ow') * _b[temp_d] / ( -2 * _b[temp2_d] )
estadd scalar tstar = r(b)[1,1]
local tse : display %5.4fc sqrt(r(V)[1,1])
estadd local tstar_se = "$ ( `tse' ) $"

estadd local pair_fe "\checkmark"
estadd local origin_fp "\checkmark"

**
** Reg 4: Reg 2, but OLS not poisson
** 
reghdfe l_trade_rat tdif tdif2 ///
	c.( lgdp_pc_o $nrg_control_o )##c.( lgdp_pc_o $nrg_control_o ) ///
	c.( lgdp_pc_d $nrg_control_d )##c.( lgdp_pc_d $nrg_control_d ) ///
	, absorb(panel_id year) cluster(iso3_o iso3_d year)
eststo y_nrg_ols_one , noesample 

nlcom -2 * _b[tdif2]
estadd scalar gam = r(b)[1,1] / 4 
local gse : display %5.4fc sqrt(r(V)[1,1]) / 4
estadd local gam_se = "$ ( `gse' ) $"

nlcom _b[tdif] / ( -2 * _b[tdif2] )
estadd scalar tstar = r(b)[1,1]
local tse : display %5.4fc sqrt(r(V)[1,1])
estadd local tstar_se = "$ ( `tse' ) $"

estadd local pair_fe "\checkmark"
estadd local origin_fp "\checkmark"

**
** create table 
** 
label var temp_o "Exporter Temperature (C)"
label var temp2_o "Exporter Temperature$ ^2 $ "
label var temp_d "Importer Temperature (C)"
label var temp2_d "Importer Temperature$ ^2 $ "
label var tdif "Exporter-Importer Temperature (C) difference"
label var tdif2 "Exporter-Importer Temperature$ ^2 $ difference "

esttab y_nrg_ols_both y_nrg_ols_one y_nrg_poisson_both y_nrg_poisson_one using "$tabs/grav_results.tex", ///
    b(%4.3fc) se(%4.3fc) star(* 0.10 ** 0.05 *** 0.01) ///
    label booktabs replace ///
	keep(temp_o temp2_o temp_d temp2_d tdif tdif2 ) numbers nomtit nonotes ///
    mgroups("OLS" "Poisson", pattern(1 0 1 0) span prefix(\multicolumn{@span}{c}{) suffix(}) erepeat(\cmidrule(lr){@span})) ///
    stats(tstar tstar_se gam gam_se pair_fe origin_fp r2 r2_p N ///
		, fmt(%4.3fc %12.0s %4.3fc %12.0s %12.0s %12.0s %4.3fc %4.3fc %9.0fc) ///
		labels("$ T^* $" "~" "$ \gamma $" "~" "Importer-Exporter pair FE" "Origin GDP/cap control" "Origin energy controls" "R$^2$" "Pseudo-R$^2$" "Observations"))
		
**
** Figure 
**

** use model 2 -- poisson with temp difference.
ppmlhdfe trade_rat tdif tdif2 ///
	c.( lgdp_pc_o $nrg_control_o )##c.( lgdp_pc_o $nrg_control_o ) ///
	c.( lgdp_pc_d $nrg_control_d )##c.( lgdp_pc_d $nrg_control_d ) ///
	, absorb(panel_id year) cluster(iso3_o iso3_d year)
	
cap drop b 
cap drop bU 
cap drop bL
cap drop tgrid 
gen tgrid = -6 + _n/4 if _n < 161 
predictnl b = 100*(exp( ( _b[tdif] + 2*_b[tdif2]*tgrid )/(4) ) - 1 ), ci(bU bL)

cap drop labflag
gen labflag = ""
replace labflag = "USA" if tgrid == 14.25
replace labflag = "Russia" if tgrid == 5.25
replace labflag = "India" if tgrid == 26

cap drop pop_w
gen pop_w = floor(pop_o/1e6)

graph twoway (rarea bU bL tgrid if tgrid < 31 , lcol(%0) fcol(ebblue%40) ) ///
	(line b tgrid if tgrid < 31, lwid(thick) lpat(solid) lcol(ebblue) ) ///
	(scatter b tgrid if !mi(labflag), msym(Oh) msiz(medlarge) mlab(labflag) mcol(black) mlabcolor(black) mlabangle(0)) ///
	(hist mtemp if iso3_o == iso3_d & year==2016  ,  start(-3) w(2)  lcol(%0) fcol(gs8%90) yaxis(2)) ///
	(hist mtemp [fw=pop_w] if iso3_o == iso3_d & year==2016  , start(-3) w(2) lcol(black%80) fcol(%0) yaxis(2)) ///
	, xtit("Country baseline temperature ({sup:o}C)") ytit("Marginal damages (%)") ///
	ysc(ra(0 0.8) axis(2)) ysc(ra(-60 60) axis(1)) ylab(-60(20)60 , axis(1) ) ylab(0(0.1)0.8 , axis(2)) ///
	yline(0 , lpat(solid) lcol(gs12)) ///
	legend(on pos(6) rows(2) order(2 "Marginal effects" 1 "95% CI" 4 "Baseline temperature distribution" 5 "Population weighted"))
graph export "$figs/me_fig.png" , width(3200) replace
	
