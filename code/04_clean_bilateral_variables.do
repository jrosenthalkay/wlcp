*
* clean bilateral trade variables
*

* ingest data
insheet using "$data/raw/release_2.1_1990_1999.csv" , clear

* pick vector of bilateral cost shifters
* nb: some of these will be constructed
*global bivarlist "colony_of_origin_ever colony_of_destination_ever common_colonizer common_legal_origin contiguity common_language same_country same_region  member_eu_joint"

* clean variables, generate i=j variables
*replace region_o = "south_east_asia" if region_o == "suth_east_asia"
*replace region_d = "south_east_asia" if region_d == "suth_east_asia"

*replace iso3_o = "SRB" if iso3_o == "YUG"
*replace iso3_d = "SRB" if iso3_d == "YUG"

*gen same_country = iso3_o == iso3_d
*gen same_region = region_o == region_d

* subset and collapse data 
*keep iso3_o iso3_d $bivarlist distance

*collapse (max) $bivarlist distance , by(iso3_o iso3_d)

* clean bilateral variables
*foreach var of varlist $bivarlist {
*replace `var' = 0 if iso3_o == iso3_d & "`var'" != "same_country"
*}
*replace distance = 1 if iso3_o == iso3_d // internal distance shouldn't matter here as will control for self-flow 

* save
*compress
save "$data/int/bilateral_costs" , replace 
