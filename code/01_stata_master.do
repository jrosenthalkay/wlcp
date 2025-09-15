/*

stata master for "The winners and losers of climate policies: a sufficient statistics approach" (economic journal)

authors: thomas bourany & jordan rosenthal-kay

stata/R code author: jordan rosenthal-kay

created 16 june 2025
last updated: 26 august 2025

*/

* set up relative paths
// >>>>>>>> SET YOUR ROOT FOLDER HERE <<<<<<<<
global main "C:/Users/yourname/path/to/wlcp"
// ^ Change this path to the root of the replication package on your machine

global code "${main}/code/stata_code"
global data "${main}/data"
global figs "${main}/output/figures"
global tabs "${main}/output/tables"

* packages required
foreach package in ftools reghdfe ppmlhdfe estout {
	cap which `package' // check package exists
	if _rc { // if package does not exist, install from ssc
		ssc install `package'
	}
}

** data cleaning and construction **

do "${code}/01_clean_energy_data"

do "${code}/02_clean_wdi.do"

do "${code}/03_clean_trade_data.do"

do "${code}/04_build_estimation_panel.do"

** damage function and supply elasticity estimation **

do "${code}/05_damage_estimation.do"

do "${code}/06_estimate_nu.do"

** input to counterfactuals **

do "${code}/07_get_trade_shares.do"

do "${code}/08_get_energy_import_export_shares.do"

do "${code}/09_create_baseline.csv"