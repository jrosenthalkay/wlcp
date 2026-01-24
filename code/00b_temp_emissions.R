
# libraries 
library(haven)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyverse)
library(countrycode)
library(ggpattern)

# --- paths -----
# >>>>>>>> SET YOUR ROOT FOLDER HERE <<<<<<<<
root <- '/Users/jordanrosenthalkay/Dropbox/GitHub/wlcp'
# ^ Change this path to the root of the replication package on your machine

output_path <- file.path(root,'output')
data <- file.path(root,'data/raw')
savepath <- file.path(output_path,'figures')

# --- load data -----
emissions <- file.path(data,'annual-co-emissions-by-region.csv') %>% read.csv()
temp_anomaly <- file.path(data,'temperature-anomaly.csv') %>% read.csv()

# --- collapse and merge data -----
emissions <- emissions %>% filter(Entity=="World") %>%
  group_by(Year) %>% 
  summarise(
    tot_emissions = sum(Annual.CO..emissions)/1e9 # convert to gt
  ) %>%
  ungroup() %>%
  mutate(cum_emissions = cumsum(tot_emissions))

temp_anomaly <- temp_anomaly %>% filter(Entity == "World") %>% group_by(Year) %>%
  summarise( 
    avg_anomaly = Global.average.temperature.anomaly.relative.to.1861.1890,
    ub_anomaly = Upper.bound.of.the.annual.temperature.anomaly..95..confidence.interval.,
    lb_anomaly = Lower.bound.of.the.annual.temperature.anomaly..95..confidence.interval.
    )%>%
  ungroup()

df <- emissions %>% left_join(temp_anomaly)

# --- plot -----

p <- df %>% filter(Year>=1900) %>%
  ggplot(aes(x = cum_emissions, y = avg_anomaly)) +
  # 95% band from the data
  geom_ribbon(aes(ymin = lb_anomaly, ymax = ub_anomaly),
              fill = "grey85", alpha = 0.25, show.legend = FALSE) +
  # observed series
  geom_line(aes(color = "Temperature anomaly with 95% CI"), linewidth = 0.8, na.rm = TRUE) +
  # linear fit
  geom_smooth(aes(color = "Linear fit"), linetype='dashed',
              method = "lm", se = FALSE, linewidth = 0.9, na.rm = TRUE) +
  labs(
    x = "Cumulative anthropogenic emissions (GtCO2) since 1750",
    y = "Global temperature rise from pre-industrial level (°C)",
    color = "",
    title = "Temperature rise vs. cumulative emissions since 1900"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_color_manual(
    values = c("Temperature anomaly with 95% CI" = "black", "Linear fit" = "steelblue"),
    breaks = c("Temperature anomaly with 95% CI", "Linear fit")
  ) + scale_x_continuous(labels = scales::label_comma()) +
  geom_point(
    data = df %>% filter(Year == 1949 | Year == 2019),
    aes(x = cum_emissions, y = avg_anomaly),
    color = "#E69F00", size = 3
  ) +
  ggrepel::geom_label_repel(
    data = df %>% filter(Year == 1949),
    aes(x = cum_emissions, y = avg_anomaly, label = "1949"),
    color = "#E69F00",
    box.padding = 0.4,
    min.segment.length = 0,
    seed = 123
  ) +
  ggrepel::geom_label_repel(
    data = df %>% filter(Year == 2019),
    aes(x = cum_emissions, y = avg_anomaly, label = "2019"),
    color = "#E69F00",
    box.padding = 0.4,
    min.segment.length = 0,
    seed = 123
  )

p

savename <- 'linear_climatesystem.png'
ggsave(file.path(savepath,savename), plot = p, width = 8, height = 6, dpi = 300)


# t-test
m <- lm(log(avg_anomaly) ~ I(log(cum_emissions/1000)), 
        data = df %>% filter(Year >= 1949))

est <- coef(m)[2]
se  <- sqrt(vcov(m)[2,2])
t_stat <- (est - 1) / se
df_resid <- df.residual(m)
p_val <- 2 * pt(-abs(t_stat), df_resid)

c(t_stat = t_stat, df = df_resid, p_value = p_val)

