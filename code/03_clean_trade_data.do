**
** set up trade data
**

insheet using "$data/raw/ITPD_E_R01.csv" , clear 

* remove trade in energy (we calibrate this sector differently)
drop if broad_sector == "Mining_Energy"

collapse (rawsum) trade , by(exporter_iso3 importer_iso3 year)

ren exporter_iso3 iso3_o
ren importer_iso3 iso3_d

duplicates report iso3_o iso3_d year 

compress 
save "$data/int/trade_collapse_noNRG" , replace
