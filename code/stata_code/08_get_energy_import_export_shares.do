
*
* get imports and exports of energy
*

* energy exports
insheet using "$data/raw/ITPD_E_R01.csv" , clear 

keep if broad_sector == "Mining_Energy"
drop if strpos(industry_descr,"ining")>0

collapse (rawsum) trade , by(exporter_iso3 year)

ren exporter_iso3 iso3 

ren trade nrg_exports
compress 
save "$data/int/nrg_export_totals" , replace

* nrg imports
insheet using "$data/raw/ITPD_E_R01.csv" , clear 

keep if broad_sector == "Mining_Energy"
drop if strpos(industry_descr,"ining")>0

collapse (rawsum) trade , by(importer_iso3 year)

ren importer_iso3 iso3 
ren trade nrg_imports 
compress 
save "$data/int/nrg_import_totals" , replace


* get ITPD-E GDP
insheet using "$data/raw/ITPD_E_R01.csv" , clear 

collapse (rawsum) trade , by(year importer_iso3)

ren trade trade_tot 
ren importer_iso3 iso3 

merge 1:1 iso3 year using "$data/int/nrg_export_totals"
merge 1:1 iso3 year using "$data/int/nrg_import_totals" , nogen

* merge on gdp
preserve 

use "$data/int/global_gdp_panel"  , clear 

gen gdp = y_gdp_per_capita * y_total_pop 
keep iso3 year gdp
tempfile gdp 
save `gdp'
display
restore 

merge 1:1 iso3 year using `gdp' , keep(1 3) nogen

* compute import / export shares. apply a clamp
gen net_imports = nrg_imports - nrg_exports

gen net_import_share1 = net_imports / trade_tot 
gen net_import_share2 = net_imports*1e6 / gdp 

egen net_import_share_a = rowmin(net_import_share1 net_import_share2)
egen net_import_share_b = rowmax(net_import_share1 net_import_share2)
gen net_import_share3 = net_import_share_a*(net_import_share1>0) + net_import_share_b*(net_import_share1<=0)

gen nrg_export_share = -net_import_share3*(net_import_share3<=0) if !mi(net_import_share3)
gen nrg_import_share = net_import_share3*(net_import_share3>0) if !mi(net_import_share3)


*
preserve 

collapse nrg_export_share nrg_import_share , by(iso3)

outsheet * using "$data/model_input/nrg_trade_share.csv"

restore


