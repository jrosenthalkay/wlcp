
### Readme file for "WLCP/code"

### How to run the code

0. `aa_clean_berkeley_temp.R` requires you to set the root. 
   - This file creates country temperature data by averaging pre-processed BEST gridded temp data using GHS population weights. This output is an input to the stata processing and estimation steps. Its output (`data/int/country_avg_temp_timeseries_ghs_popweight.csv`) is provided as part of the replication package, as this step is computationally costly: it requires handling a large populaiton raster and pre-processed data. See the README in the `data/raw` for information on these inputs.

1. `00_stata_master.do` runs all other data processing and estimation using relative paths. It requires you to set the root folder. Edit the `global main` line near the top to point to the root folder where this repository lives. All paths are defined relative to `main`.
2. `compute_suff_stats.jl` uses the output from the data processing and estimation steps to compute our welfare formulas. It also requires that you manually set the root folder name.
3. Finally `results_figures_tables_wrapper.R` creates all maps and bar graphs in the paper, using the sufficient statistics output. This file also requires you manually set the root.

**Data processing and cleaning files**

These handle initial data processing steps:

| File                            | Task                                                         | Input                                                        | Output                                                 |
| ------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------ |
| aa_clean_berkeley_temp.R        | Aggregates Berkeley Earth temperature data to the country level, using GHS-POP to form weights | /data/raw/GHS_POP_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif<br />data/raw/berktemp_1750-2019.csv | data/int/country_avg_temp_timeseries_ghs_popweight.csv |
| 01_clean_energy_data.do         | Cleans WDI data and formats to match trade data.             | data/raw/WDI_energy_data_extract.csv                         | data/int/wdi_nrg_o.dta <br />data/int/wdi_nrg_d.dta    |
| 02_clean_wdi.do                 | This file extracts and cleans GDP and population data from the World Development Indicators. For calibration, we extrapolate population and GDP for missing observations by fitting a gdp growth or population growth time trend for each country. | data/raw/WDI_main.csv<br />data/int/country_avg_temp_timeseries_ghs_popweight.csv | data/int/global_gdp_panel.dta                          |
| 03_clean_trade_data.do          | Aggregates the ITPD_E data to the country-pair-year level, omitting energy trade. Note: this takes a while to run. | data/raw/ITPD_E_R01.csv                                      | data/int/trade_collapse_noNRG.dta                      |
| 04_clean_bilateral_variables.do | Cleans the bilateral gravity variables                       | data/raw/release_2.1_1990_1999.csv                           | data/int/bilateral_costs.dta                           |
| 05_build_estimation_panel.do    | Builds the panel for damage function estimation              | data/int/trade_collapse_noNRG.dta<br />data/int/bilateral_costs.dta<br />data/int/global_gdp_panel.dta<br />data/int/wdi_nrg_o.dta<br />data/int/wdi_nrg_d.dta | data/int/trade_estimation_panel.dta                    |



**Estimation/calibration input code**

These files handle estimation of the damage function and of the energy supply elasticities.

| File                                  | Task                                                         | Input                                                        | Output                                                       |
| ------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 06_damage_estimation.do               | Estimates the damage function on productivity using the import penetration regression in a variety of ways; saves regression table, plots marginal effects figure. | data/int/trade_estimation_panel                              | output/tables/grav_results.tex<br />output/figures/me_fig.png |
| 07_estimate_nu.do                     | Estimates energy supply elasticities                         | data/raw/coal-prices.csv<br />data/raw/crude-oil-prices.csv<br />natural-gas-prices.csv<br />data/int/global_gdp_panel.dta<br />data/int/wdi_nrg_o.dta | output/figures/coalprice.png<br />data/model_input/nu_fossil.csv<br />data/model_input/nu_coal.csv<br />output/figures/supply_elas_fig.png<br />output/figures/eb_shrink.png |
| 08_get_trade_shares.do                | Creates input data for the sufficient statistics: Nat'l accounts data + trade shares | data/int/wdi_nrg_o.dta<br />data/int/global_gdp_panel.dta<br />data/int/trade_collapse_noNRG.dta | data/model_input/trade_shares.csv                            |
| 09_get_energy_import_export_shares.do | Creates a matrix of energy imports and exports               | data/raw/ITPD_E_R01.csv<br />data/int/global_gdp_panel.dta   | data/int/nrg_export_totals.dta<br />data/int/nrg_import_totals.dta<br />data/model_input/nrg_trade_share.csv |
| 10_create_baseline_csv.do             | Creates the main CSV to calibrate country level parameters in the model. Energy shares are extrapolated for countries where missing using NLLS on for sigmoid functions of country-level covariates. | data/raw/WDI_carbon.csv<br />data/int/nrg_trade_share.csv<br />data/model_input/nu_coal.csv<br />data/model_input/nu_fossil.csv<br /><br />data/raw/owid-energy-data.csv<br />data/int/wdi_nrg_o.dta<br />data/int/global_gdp_panel | data/model_input/baseline_csv_suff_stat.csv                  |



**Model and analysis code**

The sufficient statistics and some supplementary analysis is performed in julia. R code is used for visualizing the results of the counterfactuals.

| File                             | Task                                                         | Input                                                        | Output                                                       |
| -------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| compute_suff_stats.jl            | Computes the sufficient statistics for the main analysis.    | data/model_input/baseline_csv_suff_stat.csv<br />data/model_input/trade_shares.csv | output/figures/data_tradeshare_matrix_v5_25.png<br />output/figures/data_incomeshare_matrix_v5_25.png<br />output/model_output/*.csv |
| results_figures_tables_wrapper.R | Wrapper to run the contents of code/results_figures/ which makes the main maps, bar graphs, and Table 3 (Appendix) of the paper | results_figures/figs_e?_*.R<br />results_figures/make_policy_table.R | output/figures/welfare_*.png<br />output/tables/carbon_tax_table.tex |

*NB: results_figures_tables_wrapper.R calls the following:*

| File                                | Task                                   | Input                     | Output                             |
| ----------------------------------- | -------------------------------------- | ------------------------- | ---------------------------------- |
| results_figures/figs_e?_*.R         | Makes map and bar plots of the results | output/model_output/*.csv | output/figures/welfare_*.png       |
| results_figures/make_policy_table.R | creates Appendix Table 3               | output/model_output/*.csv | output/tables/carbon_tax_table.tex |

