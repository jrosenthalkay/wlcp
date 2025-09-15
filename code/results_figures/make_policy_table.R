# enviro
library(countrycode)
library(tidyr)
library(dplyr)

# load population (can be used for weighted means)
popcsv <- read.csv(file.path(dbox,'data/model_input/baseline_csv_suff_stat.csv'))
pop_threshold <- 0  

# function that shortens name strings for printing
clean_country_names <- function(vec) {
  vec <- gsub("\\s*\\(.*?\\)", "", vec)  # extract main string
  vec <- gsub(".*Macao.*", "Macao", vec)  # fix Macao
  vec <- gsub(".*Myanmar.*", "Myanmar", vec)  # fix Myanmar
  vec <- trimws(vec)  # trim whitespace
  return(vec)
}

# Function to add policy rows
add_policy_rows <- function(welfare_csv, country_csv, policy_name, output_file, adj=TRUE) {
  
  # Read data
  welfare_data <- read.csv(file.path(path, welfare_csv))
  if(adj){
    welfare_data <- read.csv(file.path(path, welfare_csv))  %>%
      mutate(across(everything(), ~ .x * 100))
    country_data <- read.csv(file.path(path, country_csv)) %>%
      mutate(across(-iso3, ~ .x * 100))
  } else {
    welfare_data <- read.csv(file.path(path, welfare_csv))
    country_data <- read.csv(file.path(path, country_csv))
  }
  
  # join population and compute rents
  country_data <- country_data %>%
    left_join(popcsv, by = "iso3") %>%
    mutate(
      d_rents = dW_e + dW_π,
      d_tot = dW_p,
      d_damages = dW_D
    ) %>%
    filter(population >= pop_threshold)
  
  # World-level stats
  dW_wld <- welfare_data$dW_wld_u
  dϵ_eq <- welfare_data$dϵ_eq
  dqf <- welfare_data$dqf
  d_damages_world <- welfare_data$dW_Dy_wld_u 
  d_nrg_world <- welfare_data$dW_e_wld_u + welfare_data$dW_π_wld_u 
  d_trade_world <- welfare_data$dW_p_wld_u 
  
  # Top 3 winners and losers by welfare
  top_winners <- country_data %>% arrange(desc(dW)) %>% slice(1:5)
  top_losers <- country_data %>% arrange(dW) %>% slice(1:5)
  
  # Averages (simple means)
  avg_winner_welfare <- mean(top_winners$dW, na.rm = TRUE)
  avg_loser_welfare <- mean(top_losers$dW, na.rm = TRUE)
  avg_winner_rents <- mean(top_winners$d_rents, na.rm = TRUE)
  avg_loser_rents <- mean(top_losers$d_rents, na.rm = TRUE)
  avg_winner_tot <- mean(top_winners$d_tot, na.rm = TRUE)
  avg_loser_tot <- mean(top_losers$d_tot, na.rm = TRUE)
  avg_winner_damages <- mean(top_winners$d_damages, na.rm = TRUE)
  avg_loser_damages <- mean(top_losers$d_damages, na.rm = TRUE)
  
  # Format numbers
  dW_wld_fmt <- sprintf("%.2f", dW_wld)
  dϵ_eq_fmt <- sprintf("%.2f", dϵ_eq)
  dqf_fmt <- sprintf("%.2f", dqf)
  d_damages_world_fmt <- sprintf("%.2f", d_damages_world)
  d_nrg_world_fmt <- sprintf("%.2f",d_nrg_world)
  d_trade_world_fmt <- sprintf("%.2f",d_trade_world)
  
  avg_winner_fmt <- sprintf("%.2f", avg_winner_welfare)
  avg_loser_fmt <- sprintf("%.2f", avg_loser_welfare)
  avg_winner_rents_fmt <- sprintf("%.2f", avg_winner_rents)
  avg_loser_rents_fmt <- sprintf("%.2f", avg_loser_rents)
  avg_winner_tot_fmt <- sprintf("%.2f", avg_winner_tot)
  avg_loser_tot_fmt <- sprintf("%.2f", avg_loser_tot)
  avg_winner_damages_fmt <- sprintf("%.2f", avg_winner_damages)
  avg_loser_damages_fmt <- sprintf("%.2f", avg_loser_damages)
  
  # Country names
  winner_names_vec <- countrycode(top_winners$iso3, "iso3c", "country.name")
  winner_names_vec <- clean_country_names(winner_names_vec)
  winner_names <- paste(winner_names_vec, collapse = ", ")
  
  loser_names_vec <- countrycode(top_losers$iso3, "iso3c", "country.name")
  loser_names_vec <- clean_country_names(loser_names_vec)
  loser_names <- paste(loser_names_vec, collapse = ", ")
  
  # Add change in epsilon and qf to policy name if not BAU
  if (policy_name != "Business-as-usual") {
    policy_name <- paste0("\\shortstack[l]{",policy_name, "\\\\ $\\Delta\\epsilon$: ", dϵ_eq_fmt, "\\%, $\\Delta q^f$: ", dqf_fmt, "\\% }")
  }
  
  # Write rows
  cat("\\multicolumn{5}{l}{\\emph{", policy_name, "}} \\\\ \n", file = output_file, append = TRUE, sep = "")
  cat("World & ", dW_wld_fmt, " & ", d_damages_world_fmt, " & ", d_nrg_world_fmt, " & ", d_trade_world_fmt, " \\\\ \n", file = output_file, append = TRUE, sep = "")
  cat("Biggest winners: ", winner_names, " & ", avg_winner_fmt, " & ", avg_winner_damages_fmt, " & ", avg_winner_rents_fmt, " & ", avg_winner_tot_fmt, " \\\\ \n", file = output_file, append = TRUE, sep = "")
  cat("Biggest losers: ", loser_names, " & ", avg_loser_fmt, " & ", avg_loser_damages_fmt, " & ", avg_loser_rents_fmt, " & ", avg_loser_tot_fmt, " \\\\ \\\\[-1.0em] \n", file = output_file, append = TRUE, sep = "")
}

# Create table structure
output_file <- file.path(output_path, "tables/carbon_tax_table.tex")
cat("\\begin{tabular}{lcccc} \n", file = output_file)
cat("~ & \\shortstack{\\% change \\\\ welfare} & \\shortstack{\\% change \\\\ climate damages} & \\shortstack{\\% change \\\\ energy effects} & \\shortstack{\\% change \\\\ trade effects} \\\\ \n", file = output_file, append = TRUE)

# Add policies
cat("\\hline \\\\[-1.0em] \n", file = output_file, append = TRUE)
add_policy_rows("welfare_climate_l_12pc_dXagg.csv", "welfare_climate_l_12pc.csv", "Business-as-usual", output_file, adj=FALSE)

cat("\\hline \\\\[-1.0em] \n", file = output_file, append = TRUE)
add_policy_rows("welfare_carbontaxWLD_l_50_dXagg.csv", "welfare_carbontaxWLD_l_50.csv", "Global carbon tax, \\$50 USD/tonne", output_file)

cat("\\hline \\\\[-1.0em] \n", file = output_file, append = TRUE)
add_policy_rows("welfare_carbontaxCHN_l_50_dXagg.csv", "welfare_carbontaxCHN_l_50.csv", "China carbon tax, \\$50 USD/tonne", output_file)

cat("\\hline \\\\[-1.0em] \n", file = output_file, append = TRUE)
add_policy_rows("welfare_carbontaxUS_l_50_dXagg.csv", "welfare_carbontaxUS_l_50.csv", "USA carbon tax, \\$50 USD/tonne", output_file)

cat("\\hline \\\\[-1.0em] \n", file = output_file, append = TRUE)
add_policy_rows("welfare_carbontax_carbontariff_EU_l_50_dXagg.csv", "welfare_carbontax_carbontariff_EU_l_50.csv", "European Union climate club", output_file)

cat("\\hline \\\\[-1.0em] \n", file = output_file, append = TRUE)
add_policy_rows("welfare_carbontax_carbontariff_As_l_50_dXagg.csv", "welfare_carbontax_carbontariff_As_l_50.csv", "ASEAN climate club", output_file)

# Close table
cat("\\hline \n", file = output_file, append = TRUE)
cat("\\end{tabular} \n", file = output_file, append = TRUE)

# Print message
cat("LaTeX table written to:", output_file, "\n")

# Print to console to see the output
readLines(output_file)
