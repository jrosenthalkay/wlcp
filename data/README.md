
### Readme file for "WLCC/data"

Files in `raw` are unprocessed inputs (all web downloads), while `int` files are intermediate products used in data processing. Files in `model_input` are the data used to compute the sufficient statistics.



Below is a description of each file in the folder: 

**`raw`**

- *berktemp_1750-2019.csv* -- 1 degree x 1 degree average annual temperature (TAVG) from [Berkeley Earth](https://berkeleyearth.org/data/). This data was transformed from a netCDF (.nc) file to a (.csv) outside of this project's workflow. 
  
  - It is omitted from the replication folder on GitHub. Please email the authors for a copy of this file if necessary; it is 110MB. It is only used as an input to construct historical country temperatures (which are provided in `data/int`) in the file `00_clean_berkeley_temp.R`. 
  
- *co2-intensity.csv* -- Carbon intensity data from Our World in Data (OWID): ["Carbon intensity: CO₂ emissions per dollar of GDP"](https://ourworldindata.org/grapher/co2-intensity?time=latest)

> Global Carbon Budget (2024); Bolt and van Zanden – Maddison Project Database 2023 – with major processing by Our World in Data. “Annual CO₂ emissions per GDP (kg per international-$)” [dataset]. Global Carbon Project, “Global Carbon Budget”; Bolt and van Zanden, “Maddison Project Database 2023” [original data].
Source: Global Carbon Budget (2024), Bolt and van Zanden – Maddison Project Database 2023 – with major processing by Our World In Data

- *coal-prices.csv* -- international coal prices from OWID: ["Coal prices"](https://ourworldindata.org/grapher/coal-prices?v=1&csvType=full&useColumnShortNames=false) downloaded 9 march 2025. OWID cites,
  
  >Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2024) – with major processing by Our World in Data. “Coal” [dataset]. Energy Institute, “Statistical Review of World Energy” [original data].
  >Source: Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2024) – with major processing by Our World In Data
  
- *crude-oil-prices.csv* -- Oil prices from Our World In Data (OWID): ["Crude oil prices"](https://ourworldindata.org/grapher/crude-oil-prices?v=1&csvType=full&useColumnShortNames=false), downloaded 8 march 2025. OWID encourages the citation for their data sources:
  
  > Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2024) – with major processing by Our World in Data. “Oil price - Crude prices since 1861” [dataset]. Energy Institute, “Statistical Review of World Energy” [original data].
  > Source: Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2024) – with major processing by Our World In Data
  
- *GHS_POP_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif* -- This is the GHS Population Raster, available from the [Global Human Settlement Layer](https://human-settlement.emergency.copernicus.eu/download.php) for 2015, global, at 30 arcsec resolution. 
  
  - Also not included due to its size (~380MB). Please download if you would like to replicate the country level temperature construction in `00_clean_berkeley_temp.R`. 
  
- *ITPD_E_R01.csv* -- This is release 1 of the ITPD estimation data from the [US Gravity portal](https://www.usitc.gov/data/gravity/itpde.htm). 
  
  - It is 4GB so we do not include it in the replication package.
  
- *natural-gas-prices.csv* -- again from OWID: ["Natural gas prices"](https://ourworldindata.org/grapher/natural-gas-prices?v=1&csvType=full&useColumnShortNames=false) downloaded 8 march 2025. OWID cites,
  
  > Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2024) – with major processing by Our World in Data. “Gas price” [dataset]. Energy Institute, “Statistical Review of World Energy” [original data].
  > Source: Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2024) – with major processing by Our World In Data
  
- *owid-energy-data.csv* -- data on nations' energy mix from OWID.

- *WDI_carbon.csv* -- carbon emissions data from the WDI; in particular the variable `Carbon dioxide (CO2) emissions (total) excluding LULUCF (Mt CO2e)` .

- *WDI_energy_data_extract.csv* -- This is a simple extract of the following variables from the [World Development Indicators](https://databank.worldbank.org/source/world-development-indicators).
- *WDI_main.csv* -- standard WDI data (GDP, mortality, etc)



**`int`**

- *wdi_nrg_d.dta* and *wdi_nrg_o.dta* -- these are cleaned energy data files to be merged onto the trade panel.
- *trade_collapsed_noNRG.dta* -- this is the aggregated total volume of trade net energy trade from the ITPD data.
- *country_avg_temp_timeseries_ghs_popweight.csv* -- this is the panel of temperature at the county-year level
- *global_gdp_panel.csv* country level panel of GDP and population; for estimation and calibration
- *bilateral_costs.csv* -- collapsed and cleaned standard gravity variables for projecting bilateral trade costs
- *trade_estimation_panel.dta* -- cleaned bilateral trade panel with energy controls, GDP and population, and bilateral trade costs.



**`model_input`**

- *nu_fossil.csv* -- $\nu^f$ estimates
- *nu_coal.csv* -- $\nu^c$ estimates
- *trade_shares.csv* -- the goods trade shares matrix, $\mathbf{T}$. 
- *nrg_trade_share.csv* -- trade shares in energy.
- *baseline_csv_suff_stat.csv* -- main model input: country level energy mix, gdp, population, etc augmented with supply elasticities, energy rents, and so on.


**`model_output`**

contains .csvs created for the different policy experiments, produced by the julia script, `code/02_compute_suff_stats.jl`. These files are used to produce the figures and tables handled by `code/03_results_figures_tables_wrapper.R`. 
